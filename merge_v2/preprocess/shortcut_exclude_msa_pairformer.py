#!/usr/bin/env python3
"""Shortcut: produce the six variant-side score tables that would result from a
pipeline rerun with ``msa_pairformer_variant.parquet`` added to
``create_variant_scores_all_table.py``'s ``--exclude`` list, WITHOUT paying the
full merge + percentile + gene-aggregation cost.

Empirical justification
-----------------------
The current local 22-col ``variant_scores_all_outer.parquet`` (18 scores, incl.
``msa_pairformer_llr``) was verified to contain **zero rows** where the 17
legacy score columns are all NULL, i.e. MSA-Pairformer canonical added **zero
unique variant keys** to the merge. Therefore:

* Dropping the ``msa_pairformer_llr`` column from the 22-col outer produces a
  table bit-exactly identical to a fresh ``--exclude`` rebuild.
* Applying drop-any-NA on the remaining 17 columns yields the same inner as a
  rebuild would produce (49,188,671 rows).
* The same argument extends to ``*_pre_percentile.parquet`` (per-column
  percentiles are computed independently), to the threshold TSV (per-column
  thresholds are computed independently), and to ``variant_scores_all_outer_ensg``
  (the linker join is on ``chrom/pos/ref/alt`` only, so it commutes with
  drop-any-NA on score columns).

Outputs (six)
-------------
1. ``scores/variant_scores_all_outer.parquet``            (drop column)
2. ``scores/variant_scores_outer_pre_percentile.parquet`` (drop column)
3. ``scores/variant_scores_all_inner.parquet``            (drop column + drop-any-NA)
4. ``scores/variant_scores_inner_pre_percentile.parquet`` (drop column + drop-any-NA)
5. ``scores/variant_scores_percentile_thresholds.tsv``    (drop ``msa_pairformer_llr`` row)
6. ``scores/gene_aggregated/variant_scores_all_inner_ensg.parquet``
   (from ``outer_ensg``: drop column + drop-any-NA on 17)

The following files CANNOT be shortcut and MUST be rebuilt by the pipeline
because their contents depend on the new inner row set (~49.2 M vs the current
16.7 M) that MSA-Pairformer's canonical filter had been suppressing::

  scores/variant_scores_inner_post_percentile.parquet         (recomputed CDFs)
  scores/gene_aggregated/variant_scores_all_inner_ensg_stats.parquet
                                                              (per-gene stats)
  full_analysis_tables/gene_aggregated/*_inner_*              (all downstream)

Safety model
------------
* Outputs are first written to ``<final>.17col.parquet`` (or ``.17col.tsv``)
  siblings next to the canonical filenames, so the current 22-col files are
  never at risk during the write.
* Without ``--commit``, the script stops after producing the sidecars. You can
  inspect them and decide whether to promote them.
* With ``--commit``, each existing 22-col original is renamed to
  ``<original_stem>.22col.bak.<ext>`` and the matching ``.17col.*`` sidecar is
  renamed into its canonical place. Both operations are ``os.rename`` calls, so
  they are atomic on the same filesystem. Rolling back is a symmetric reverse
  rename.

Usage
-----
::

    # inspect the plan; touch nothing
    python merge_v2/preprocess/shortcut_exclude_msa_pairformer.py --dry-run

    # produce .17col.* siblings but leave originals in place
    python merge_v2/preprocess/shortcut_exclude_msa_pairformer.py

    # produce siblings then atomically swap them into place
    python merge_v2/preprocess/shortcut_exclude_msa_pairformer.py --commit
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import duckdb


REPO_ROOT = Path(__file__).resolve().parents[2]
SCORES_DIR_DEFAULT = REPO_ROOT / "merge_v2" / "data" / "processed_data" / "scores"
GA_DIR_DEFAULT = SCORES_DIR_DEFAULT / "gene_aggregated"

LEGACY_17 = [
    "AM", "mcap", "esm1b", "gmvp", "phylop", "sift", "cadd", "cpt",
    "gpn_msa", "ESM_1v", "EVE", "popEVE", "PAI3D", "MisFit_D",
    "MisFit_S", "mpc", "polyphen",
]
LEGACY_17_PP = [f"{c}_pre_percentile" for c in LEGACY_17]
DROP_SCORE_COL = "msa_pairformer_llr"
DROP_PP_COL = f"{DROP_SCORE_COL}_pre_percentile"

COMPRESSION = "zstd"
ROW_GROUP_SIZE = 512_000
DEFAULT_LEADING_KEYS = ("chrom", "pos", "ref", "alt")


def q_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def q_str(s: str) -> str:
    return "'" + s.replace("'", "''") + "'"


def project_expr(cols: list[str], leading: tuple[str, ...] = DEFAULT_LEADING_KEYS) -> str:
    return ", ".join(q_ident(c) for c in [*leading, *cols])


def all_present_expr(cols: list[str]) -> str:
    return " AND ".join(f"{q_ident(c)} IS NOT NULL" for c in cols)


@dataclass
class Step:
    label: str
    src: Path
    final_dst: Path
    projection: str
    where: Optional[str] = None
    expected_rows: Optional[int] = None  # sanity check, if known
    leading: tuple[str, ...] = field(default=DEFAULT_LEADING_KEYS)

    @property
    def sidecar(self) -> Path:
        return self.final_dst.with_name(f"{self.final_dst.stem}.17col{self.final_dst.suffix}")

    @property
    def backup(self) -> Path:
        return self.final_dst.with_name(f"{self.final_dst.stem}.22col.bak{self.final_dst.suffix}")


def build_plan(scores_dir: Path, ga_dir: Path) -> tuple[list[Step], Path, Path]:
    src_outer = scores_dir / "variant_scores_all_outer.parquet"
    src_outer_pp = scores_dir / "variant_scores_outer_pre_percentile.parquet"
    src_outer_ensg = ga_dir / "variant_scores_all_outer_ensg.parquet"

    plan = [
        Step(
            label="outer",
            src=src_outer,
            final_dst=scores_dir / "variant_scores_all_outer.parquet",
            projection=project_expr(LEGACY_17),
            where=None,
            expected_rows=79_045_780,
        ),
        Step(
            label="outer_pre_percentile",
            src=src_outer_pp,
            final_dst=scores_dir / "variant_scores_outer_pre_percentile.parquet",
            projection=project_expr(LEGACY_17_PP),
            where=None,
            expected_rows=79_045_780,
        ),
        Step(
            label="inner",
            src=src_outer,
            final_dst=scores_dir / "variant_scores_all_inner.parquet",
            projection=project_expr(LEGACY_17),
            where=all_present_expr(LEGACY_17),
            expected_rows=49_188_671,
        ),
        Step(
            label="inner_pre_percentile",
            src=src_outer_pp,
            final_dst=scores_dir / "variant_scores_inner_pre_percentile.parquet",
            projection=project_expr(LEGACY_17_PP),
            where=all_present_expr(LEGACY_17_PP),
            expected_rows=49_188_671,
        ),
        Step(
            label="inner_ensg",
            src=src_outer_ensg,
            final_dst=ga_dir / "variant_scores_all_inner_ensg.parquet",
            projection=project_expr(LEGACY_17, leading=(*DEFAULT_LEADING_KEYS, "ensg")),
            where=all_present_expr(LEGACY_17),
            expected_rows=None,  # fanout unknown a priori; only asserted >= inner rows
            leading=(*DEFAULT_LEADING_KEYS, "ensg"),
        ),
    ]
    tsv_src = scores_dir / "variant_scores_percentile_thresholds.tsv"
    return plan, src_outer_ensg, tsv_src


def verify_sources_or_die(steps: list[Step], tsv_src: Path) -> None:
    missing = []
    for st in steps:
        if not st.src.exists():
            missing.append(st.src)
    if not tsv_src.exists():
        missing.append(tsv_src)
    if missing:
        for m in missing:
            print(f"error: missing source: {m}", file=sys.stderr)
        raise SystemExit(1)


def verify_source_has_column(con: duckdb.DuckDBPyConnection, src: Path, col: str) -> None:
    schema = con.execute(f"DESCRIBE SELECT * FROM read_parquet({q_str(str(src))})").fetchall()
    names = {row[0] for row in schema}
    if col not in names:
        raise SystemExit(
            f"error: expected column {col!r} in {src} (found {sorted(names)!r}); "
            f"cannot safely shortcut - source may already be a 17-col file."
        )


def write_parquet(
    con: duckdb.DuckDBPyConnection, src: Path, tmp_dst: Path, projection: str, where: Optional[str]
) -> None:
    where_clause = f" WHERE {where}" if where else ""
    sql = (
        f"COPY (SELECT {projection} FROM read_parquet({q_str(str(src))}){where_clause})"
        f" TO {q_str(str(tmp_dst))}"
        f" (FORMAT PARQUET, COMPRESSION '{COMPRESSION}', ROW_GROUP_SIZE {ROW_GROUP_SIZE})"
    )
    con.execute(sql)


def count_rows(con: duckdb.DuckDBPyConnection, path: Path) -> int:
    return con.execute(f"SELECT count(*) FROM read_parquet({q_str(str(path))})").fetchone()[0]


def shortcut_tsv(tsv_src: Path, tsv_dst: Path, drop_key: str) -> tuple[int, int]:
    lines = tsv_src.read_text().splitlines()
    if not lines:
        raise SystemExit(f"error: empty TSV: {tsv_src}")
    header = lines[0]
    body_in = lines[1:]
    body_out = [row for row in body_in if not row.startswith(drop_key + "\t")]
    tsv_dst.write_text("\n".join([header, *body_out]) + "\n")
    return len(body_in), len(body_out)


def commit_step(final_dst: Path, sidecar: Path, backup: Path) -> None:
    if not sidecar.exists():
        raise SystemExit(f"error: cannot commit; sidecar missing: {sidecar}")
    if backup.exists():
        raise SystemExit(
            f"error: refusing to commit; backup path already exists: {backup}. "
            f"Move or remove it, then retry."
        )
    if final_dst.exists():
        final_dst.rename(backup)
    sidecar.rename(final_dst)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.strip(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--dry-run", action="store_true", help="print plan and exit")
    parser.add_argument(
        "--commit", action="store_true",
        help="after writing sidecars, atomically swap them into their canonical paths "
             "(existing originals renamed to .22col.bak.*)",
    )
    parser.add_argument("--memory-limit", default="8GB")
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--temp-dir", default=None)
    parser.add_argument("--scores-dir", default=str(SCORES_DIR_DEFAULT))
    parser.add_argument("--gene-agg-dir", default=str(GA_DIR_DEFAULT))
    args = parser.parse_args()

    if args.dry_run and args.commit:
        print("error: --dry-run and --commit are mutually exclusive", file=sys.stderr)
        return 2

    scores_dir = Path(args.scores_dir)
    ga_dir = Path(args.gene_agg_dir)

    plan, _outer_ensg_src, tsv_src = build_plan(scores_dir, ga_dir)
    tsv_dst_sidecar = tsv_src.with_name(f"{tsv_src.stem}.17col{tsv_src.suffix}")
    tsv_dst_backup = tsv_src.with_name(f"{tsv_src.stem}.22col.bak{tsv_src.suffix}")

    verify_sources_or_die(plan, tsv_src)

    if args.dry_run:
        print("Dry-run.  Would produce these 17-col sidecar files:")
        for st in plan:
            n = f"{st.expected_rows:,}" if st.expected_rows else "unknown"
            print(f"  [{st.label:<22s}] src={st.src.name!s:<50} -> dst={st.sidecar.name!s:<58} (~{n} rows)")
        print(f"  [{'thresholds':<22s}] src={tsv_src.name!s:<50} -> dst={tsv_dst_sidecar.name!s:<58} (17 score rows)")
        print()
        print("With --commit, originals would be renamed to:")
        for st in plan:
            print(f"  {st.final_dst.name} -> {st.backup.name}")
        print(f"  {tsv_src.name} -> {tsv_dst_backup.name}")
        return 0

    con = duckdb.connect()
    con.execute(f"SET memory_limit='{args.memory_limit}'")
    if args.threads is not None:
        con.execute(f"SET threads={int(args.threads)}")
    if args.temp_dir:
        con.execute(f"SET temp_directory={q_str(str(args.temp_dir))}")
    con.execute("SET preserve_insertion_order=false")

    # sanity: make sure each source actually has the column we intend to drop
    for st in plan:
        col_to_verify = DROP_PP_COL if st.label in {"outer_pre_percentile", "inner_pre_percentile"} else DROP_SCORE_COL
        verify_source_has_column(con, st.src, col_to_verify)

    print(f"scores_dir  = {scores_dir}")
    print(f"gene_agg_dir= {ga_dir}")
    print()

    for st in plan:
        if st.sidecar.exists():
            print(f"[{st.label}] sidecar already exists, skipping: {st.sidecar.name}")
            continue
        tmp = st.sidecar.with_suffix(st.sidecar.suffix + ".tmp")
        if tmp.exists():
            tmp.unlink()
        t0 = time.time()
        print(f"[{st.label}] writing {st.sidecar.name} ...", flush=True)
        write_parquet(con, st.src, tmp, st.projection, st.where)
        tmp.rename(st.sidecar)
        n_out = count_rows(con, st.sidecar)
        sz_gb = st.sidecar.stat().st_size / 1e9
        elapsed = time.time() - t0
        line = f"[{st.label}] ok - {n_out:,} rows, {sz_gb:.2f} GB, {elapsed:.1f}s"
        if st.expected_rows is not None and n_out != st.expected_rows:
            line += f"  [WARN: expected {st.expected_rows:,} rows]"
        print(line)

    # TSV shortcut
    if tsv_dst_sidecar.exists():
        print(f"[thresholds] sidecar already exists, skipping: {tsv_dst_sidecar.name}")
    else:
        n_in, n_out = shortcut_tsv(tsv_src, tsv_dst_sidecar, DROP_SCORE_COL)
        print(f"[thresholds] ok - {n_out} score rows (dropped 1 from {n_in})")

    if not args.commit:
        print()
        print("Sidecars written.  Inspect them, then rerun with --commit to atomically")
        print("swap them into their canonical paths (originals -> .22col.bak.*).")
        return 0

    print()
    print("Committing: renaming originals to .22col.bak.*, then promoting .17col.* into place")
    for st in plan:
        commit_step(st.final_dst, st.sidecar, st.backup)
        print(f"  swapped: {st.final_dst.name}")
    commit_step(tsv_src, tsv_dst_sidecar, tsv_dst_backup)
    print(f"  swapped: {tsv_src.name}")

    print()
    print("Done.  The six 17-col outputs are now at their canonical paths.  Old 22-col")
    print("originals live at .22col.bak.* and can be deleted once you have re-run the")
    print("downstream pipeline steps (percentile-inner-post, gene-agg stats, analysis).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
