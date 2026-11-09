#!/usr/bin/env python3
"""Sanity-check MSA-Pairformer LLR against the missense/UniProt linker on chr22.

Purpose: before committing to a full protein->variant coalescing pipeline for
``msa_pairformer_llr_chunks0-9.parquet``, quantify

  1. Coverage: what fraction of MSA-Pairformer's chr22 UniProt accessions are
     present in the ``linker_missense_enst_transcript_aa_uniprot`` chr22 shard.
  2. Fan-out: how many linker rows exist per
     ``(uniprot_id, aa_pos, aa_ref, aa_alt)`` after joining.
  3. Isoform disagreement: for how many ``(uniprot_id, aa_pos)`` tuples do
     different isoforms in the linker disagree on ``aa_ref`` (which would mean
     joining on bare ``uniprot_id`` is unsafe)?
  4. Canonical vs. MANE vs. any: for a coalesced-to-variant table matching the
     AlphaMissense schema (llr, llr_mane_max, llr_canon_max, llr_any_max), how
     many chr22 SNVs would each policy retain?

The output is a human-readable report only; no files are produced. Run once,
share the numbers with the score's author, decide the final isoform policy,
then use ``preprocess_msa_pairformer_llr.py`` to build the actual table.

Run from ``merge_v2/preprocess/``::

    python msa_pairformer_chr22_sanity_check.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent  # merge_v2/

SCORE_PATH = PROJECT_DIR / "data" / "raw_data" / "scores" / "msa_pairformer_llr_chunks0-9.parquet"
LINKER_CHR22 = (
    PROJECT_DIR / "data" / "raw_data" / "linker"
    / "linker_missense_enst_transcript_aa_uniprot_by_chrom_tsv"
    / "linker_missense_enst_transcript_aa_uniprot_22.tsv.bgz"
)


def q(v: str) -> str:
    return "'" + v.replace("'", "''") + "'"


def main() -> int:
    for path, label in ((SCORE_PATH, "MSA-Pairformer score"), (LINKER_CHR22, "chr22 linker shard")):
        if not path.exists():
            sys.exit(f"ERROR: {label} not found: {path}")

    con = duckdb.connect()
    con.execute("SET preserve_insertion_order = false")

    # Materialize the chr22 linker into a temp table (DuckDB reads .bgz via gzip
    # transparently; the shard is ~55 MB compressed).
    print("Loading chr22 linker shard ...", flush=True)
    con.execute(
        "CREATE TEMP TABLE linker22 AS "
        f"SELECT * FROM read_csv({q(str(LINKER_CHR22))}, delim='\t', header=true, "
        "compression='gzip', nullstr=['NA',''], "
        "types={'chrom':'VARCHAR','pos':'BIGINT','ref':'VARCHAR','alt':'VARCHAR',"
        "'aa_pos':'INTEGER','aa_ref':'VARCHAR','aa_alt':'VARCHAR',"
        "'gene_symbol':'VARCHAR','enst':'VARCHAR','ensg':'VARCHAR','ensp':'VARCHAR',"
        "'uniprot_id':'VARCHAR','uniprot_isoform':'VARCHAR',"
        "'mane_select':'BOOLEAN','canonical':'BOOLEAN','transcript_mane_select':'BOOLEAN',"
        "'locus':'VARCHAR','alleles':'VARCHAR'})"
    )

    n_linker22 = con.execute("SELECT count(*) FROM linker22").fetchone()[0]
    print(f"  chr22 linker rows: {n_linker22:,}")

    # -------------------------------------------------------------------------
    # 1. Coverage: chr22 uniprot_ids in each source
    # -------------------------------------------------------------------------
    print("\n=== 1. Coverage of chr22 UniProt accessions ===")
    n_uid_linker22 = con.execute(
        "SELECT count(DISTINCT uniprot_id) FROM linker22"
    ).fetchone()[0]
    print(f"  distinct uniprot_id on chr22 (linker): {n_uid_linker22:,}")

    con.execute(
        "CREATE TEMP TABLE score_uids AS "
        f"SELECT DISTINCT uniprot_id FROM read_parquet({q(str(SCORE_PATH))})"
    )
    n_uid_score = con.execute("SELECT count(*) FROM score_uids").fetchone()[0]
    print(f"  distinct uniprot_id (MSA-Pairformer, all chroms): {n_uid_score:,}")

    # chr22 overlap
    con.execute(
        "CREATE TEMP TABLE chr22_uids AS "
        "SELECT DISTINCT l.uniprot_id "
        "FROM linker22 l SEMI JOIN score_uids s USING (uniprot_id)"
    )
    n_chr22_overlap = con.execute("SELECT count(*) FROM chr22_uids").fetchone()[0]
    n_chr22_missing = con.execute(
        "SELECT count(DISTINCT uniprot_id) FROM linker22 "
        "WHERE uniprot_id NOT IN (SELECT uniprot_id FROM score_uids)"
    ).fetchone()[0]
    print(
        f"  chr22 uniprot_ids present in BOTH: {n_chr22_overlap:,} / "
        f"{n_uid_linker22:,} linker-chr22 ({100*n_chr22_overlap/max(n_uid_linker22,1):.1f}%)"
    )
    print(f"  chr22 uniprot_ids in linker but NOT in MSA-Pairformer: {n_chr22_missing:,}")

    # -------------------------------------------------------------------------
    # 2. Isoform structure of the linker (chr22)
    # -------------------------------------------------------------------------
    print("\n=== 2. Isoform structure of the chr22 linker ===")
    row = con.execute(
        "SELECT "
        "  count(*) AS n_rows, "
        "  count(DISTINCT (uniprot_id, uniprot_isoform)) AS n_uid_isoforms, "
        "  count(DISTINCT uniprot_id) AS n_uid, "
        "  count(*) FILTER (WHERE uniprot_isoform = uniprot_id) AS n_canonical_iso, "
        "  count(*) FILTER (WHERE canonical) AS n_canonical_transcript, "
        "  count(*) FILTER (WHERE mane_select) AS n_mane, "
        "  count(*) FILTER (WHERE canonical AND mane_select) AS n_canon_and_mane "
        "FROM linker22"
    ).fetchone()
    (n_rows, n_uid_isof, n_uid, n_can_iso, n_can_tx, n_mane, n_can_and_mane) = row
    print(f"  rows: {n_rows:,}")
    print(f"  distinct (uniprot_id, uniprot_isoform): {n_uid_isof:,}")
    print(f"  distinct uniprot_id (bare accession): {n_uid:,}")
    print(f"  rows where uniprot_isoform == uniprot_id (canonical isoform tag): {n_can_iso:,} ({100*n_can_iso/n_rows:.1f}%)")
    print(f"  rows where transcript.canonical == TRUE: {n_can_tx:,} ({100*n_can_tx/n_rows:.1f}%)")
    print(f"  rows where mane_select == TRUE: {n_mane:,} ({100*n_mane/n_rows:.1f}%)")
    print(f"  rows where canonical AND mane_select: {n_can_and_mane:,}")

    # Isoforms per uniprot_id
    print("\n  isoforms per uniprot_id (chr22):")
    for row in con.execute(
        "SELECT n_isoforms, count(*) AS n_uids FROM ("
        "  SELECT uniprot_id, count(DISTINCT uniprot_isoform) AS n_isoforms "
        "  FROM linker22 GROUP BY uniprot_id"
        ") GROUP BY n_isoforms ORDER BY n_isoforms"
    ).fetchall():
        print(f"    {row[0]} isoform(s): {row[1]:,} uniprot_ids")

    # -------------------------------------------------------------------------
    # 3. Do isoforms disagree on aa_ref at the same aa_pos?
    # -------------------------------------------------------------------------
    print("\n=== 3. Isoform disagreement on aa_ref ===")
    # For each (uniprot_id, aa_pos), count distinct aa_ref values across isoforms.
    row = con.execute(
        "WITH uap AS ("
        "  SELECT uniprot_id, aa_pos, count(DISTINCT aa_ref) AS n_refs "
        "  FROM linker22 GROUP BY 1,2"
        ") "
        "SELECT count(*) AS n_uap, "
        "       count(*) FILTER (WHERE n_refs > 1) AS n_disagree "
        "FROM uap"
    ).fetchone()
    n_uap, n_disagree = row
    print(
        f"  distinct (uniprot_id, aa_pos): {n_uap:,}; "
        f"with disagreeing aa_ref across isoforms: {n_disagree:,} "
        f"({100*n_disagree/max(n_uap,1):.3f}%)"
    )

    # -------------------------------------------------------------------------
    # 4. Join with MSA-Pairformer (score is at bare-uniprot_id grain)
    # -------------------------------------------------------------------------
    print("\n=== 4. Joining MSA-Pairformer chr22 subset -> chr22 linker ===")
    # Filter MSA-Pairformer to just the chr22-relevant uniprot_ids to keep the join small.
    con.execute(
        "CREATE TEMP TABLE score22 AS "
        "SELECT s.uniprot_id, s.position AS aa_pos, s.aa_ref, s.aa_alt, s.llr "
        f"FROM read_parquet({q(str(SCORE_PATH))}) AS s "
        "SEMI JOIN chr22_uids u USING (uniprot_id)"
    )
    n_score22 = con.execute("SELECT count(*) FROM score22").fetchone()[0]
    print(f"  MSA-Pairformer rows for chr22 uniprot_ids: {n_score22:,}")

    # Rows in the score whose (uniprot_id, aa_pos, aa_ref, aa_alt) has AT LEAST
    # ONE matching linker row (join on bare uniprot_id, without any isoform filter).
    con.execute(
        "CREATE TEMP TABLE join_all AS "
        "SELECT s.uniprot_id, s.aa_pos, s.aa_ref, s.aa_alt, s.llr, "
        "       l.chrom, l.pos, l.ref, l.alt, l.enst, l.uniprot_isoform, "
        "       l.canonical, l.mane_select "
        "FROM score22 s JOIN linker22 l "
        "  ON s.uniprot_id = l.uniprot_id AND s.aa_pos = l.aa_pos "
        " AND s.aa_ref = l.aa_ref AND s.aa_alt = l.aa_alt"
    )
    n_join_all = con.execute("SELECT count(*) FROM join_all").fetchone()[0]
    n_score22_matched = con.execute(
        "SELECT count(DISTINCT (uniprot_id, aa_pos, aa_ref, aa_alt)) FROM join_all"
    ).fetchone()[0]
    print(
        f"  join-all rows (transcript+isoform fan-out): {n_join_all:,}; "
        f"unique (uniprot_id, aa_pos, aa_ref, aa_alt) matched: {n_score22_matched:,} "
        f"(match rate: {100*n_score22_matched/max(n_score22,1):.2f}%)"
    )
    n_score22_unmatched = n_score22 - n_score22_matched
    print(f"  score rows with NO linker match at all: {n_score22_unmatched:,}")

    # Fan-out: rows-per-score-key
    print("\n  fan-out distribution (linker rows per score key):")
    for row in con.execute(
        "SELECT n_rows, count(*) AS n_keys FROM ("
        "  SELECT count(*) AS n_rows FROM join_all "
        "  GROUP BY uniprot_id, aa_pos, aa_ref, aa_alt"
        ") GROUP BY n_rows ORDER BY n_rows LIMIT 20"
    ).fetchall():
        print(f"    {row[0]} linker row(s): {row[1]:,} score keys")

    # Distinct SNV per score key: how many (chrom,pos,ref,alt) does each score
    # key map to (usually 1 codon => up to 3 SNVs per amino-acid substitution)?
    print("\n  distinct SNVs per score key (usually 1 codon = 1 SNV per aa change):")
    for row in con.execute(
        "SELECT n_snv, count(*) AS n_keys FROM ("
        "  SELECT count(DISTINCT (chrom, pos, ref, alt)) AS n_snv FROM join_all "
        "  GROUP BY uniprot_id, aa_pos, aa_ref, aa_alt"
        ") GROUP BY n_snv ORDER BY n_snv LIMIT 10"
    ).fetchall():
        print(f"    {row[0]} SNV(s): {row[1]:,} score keys")

    # -------------------------------------------------------------------------
    # 5. Match rate under different isoform policies
    # -------------------------------------------------------------------------
    print("\n=== 5. Coalesced-to-SNV match rate under alternate isoform policies ===")
    # Policy A: canonical-only linker (uniprot_isoform IS NULL OR ends with -1
    # OR uniprot_isoform == uniprot_id). Compute SNVs retained under each and
    # how many the score column would be populated for.
    for policy_label, policy_sql in [
        ("any isoform (join_all)", "TRUE"),
        ("canonical isoform (uniprot_isoform IS NULL OR uniprot_isoform = uniprot_id)",
         "l.uniprot_isoform IS NULL OR l.uniprot_isoform = l.uniprot_id"),
        ("transcript.canonical = TRUE", "l.canonical"),
        ("mane_select = TRUE", "l.mane_select"),
    ]:
        row = con.execute(
            "SELECT "
            "  count(DISTINCT (l.chrom, l.pos, l.ref, l.alt)) AS n_snv, "
            "  count(DISTINCT (s.uniprot_id, s.aa_pos, s.aa_ref, s.aa_alt)) AS n_keys "
            "FROM score22 s JOIN linker22 l "
            "  ON s.uniprot_id = l.uniprot_id AND s.aa_pos = l.aa_pos "
            " AND s.aa_ref = l.aa_ref AND s.aa_alt = l.aa_alt "
            f"WHERE {policy_sql}"
        ).fetchone()
        print(f"  {policy_label}: {row[0]:,} distinct SNVs, {row[1]:,} matched score keys")

    # -------------------------------------------------------------------------
    # 6. AA sequence agreement check: pick a few uniprot_ids and reconstruct
    # aa_ref sequence from each side
    # -------------------------------------------------------------------------
    print("\n=== 6. aa_ref sequence agreement (spot check on 5 random chr22 uniprot_ids) ===")
    sample_ids = con.execute(
        "SELECT uniprot_id FROM chr22_uids USING SAMPLE 5"
    ).fetchall()
    for (uid,) in sample_ids:
        # aa_ref sequence in score, per-position (must be same at any given position within a uid)
        score_seq_rows = con.execute(
            "SELECT aa_pos, count(DISTINCT aa_ref) FROM score22 "
            f"WHERE uniprot_id = {q(uid)} GROUP BY aa_pos"
        ).fetchall()
        score_max_ref = max((r[1] for r in score_seq_rows), default=0)

        # aa_ref sequence in linker canonical isoform
        linker_seq_rows = con.execute(
            "SELECT aa_pos, count(DISTINCT aa_ref) FROM linker22 "
            "WHERE (uniprot_isoform IS NULL OR uniprot_isoform = uniprot_id) "
            f"  AND uniprot_id = {q(uid)} "
            "GROUP BY aa_pos"
        ).fetchall()
        linker_max_ref = max((r[1] for r in linker_seq_rows), default=0)

        # Overlap comparison: at positions present in both, do they agree on aa_ref?
        row = con.execute(
            "WITH s AS ("
            "  SELECT DISTINCT aa_pos, aa_ref FROM score22 "
            f"  WHERE uniprot_id = {q(uid)}"
            "), l AS ("
            "  SELECT DISTINCT aa_pos, aa_ref FROM linker22 "
            "  WHERE (uniprot_isoform IS NULL OR uniprot_isoform = uniprot_id) "
            f"    AND uniprot_id = {q(uid)}"
            ") "
            "SELECT count(*) AS n_pos_both, "
            "       count(*) FILTER (WHERE s.aa_ref = l.aa_ref) AS n_agree, "
            "       count(*) FILTER (WHERE s.aa_ref <> l.aa_ref) AS n_disagree "
            "FROM s FULL OUTER JOIN l USING (aa_pos) "
            "WHERE s.aa_ref IS NOT NULL AND l.aa_ref IS NOT NULL"
        ).fetchone()
        n_pos_both, n_agree, n_disagree = row
        # Also compute score-only and linker-only positions
        n_score_only = con.execute(
            f"SELECT count(DISTINCT aa_pos) FROM score22 WHERE uniprot_id = {q(uid)} "
            "AND aa_pos NOT IN (SELECT aa_pos FROM linker22 "
            "  WHERE (uniprot_isoform IS NULL OR uniprot_isoform = uniprot_id) "
            f"    AND uniprot_id = {q(uid)})"
        ).fetchone()[0]
        n_linker_only = con.execute(
            "SELECT count(DISTINCT aa_pos) FROM linker22 "
            "WHERE (uniprot_isoform IS NULL OR uniprot_isoform = uniprot_id) "
            f"  AND uniprot_id = {q(uid)} "
            f"  AND aa_pos NOT IN (SELECT aa_pos FROM score22 WHERE uniprot_id = {q(uid)})"
        ).fetchone()[0]
        print(
            f"  {uid}: overlap_pos={n_pos_both}, agree={n_agree}, disagree={n_disagree}, "
            f"score-only={n_score_only}, canonical-linker-only={n_linker_only}, "
            f"score_multiref_at_pos_max={score_max_ref}, canonical_linker_multiref_max={linker_max_ref}"
        )

    con.close()
    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
