#!/usr/bin/env python3
"""Coalesce MSA-Pairformer LLR (protein grain) into a variant-grain score parquet.

The raw score at

    merge_v2/data/raw_data/scores/msa_pairformer_llr_chunks0-9.parquet

is keyed on ``(uniprot_id, position, aa_ref, aa_alt)`` (~193.5M rows) -- one row
per amino-acid substitution per UniProt accession. To feed it into the standard
``create_variant_scores_all_table.py`` merge, it must be projected down to
``(chrom, pos, ref, alt)`` grain via the missense/UniProt linker under

    merge_v2/data/raw_data/linker/linker_missense_enst_transcript_aa_uniprot_by_chrom_tsv/

which supplies ``(chrom, pos, ref, alt, aa_pos, aa_ref, aa_alt, uniprot_id,
uniprot_isoform, mane_select, ...)`` per missense variant + transcript + isoform.
Because a single ``(uniprot_id, aa_pos, aa_ref, aa_alt)`` score row can match
multiple isoforms of the same UniProt accession -- and ~11% of ``(uniprot_id,
aa_pos)`` tuples in the linker have disagreeing ``aa_ref`` across isoforms
(measured on chr22; see the sanity check next to this script) -- the fold-down
must pick an explicit isoform policy. **This script uses the canonical UniProt
isoform policy**: rows are kept only where the linker's ``uniprot_isoform`` is
NULL or equal to ``uniprot_id`` (i.e. the linker's designation of "canonical
by default"; matches the sanity check's convention). Historical iterations of
this preprocessor also emitted MANE-select and any-isoform columns; those were
retired because the canonical UniProt isoform is the endorsed reduction for
downstream analyses. If the other policies are ever needed the earlier revision
is in the git history.

The output schema matches the "minimal" variant-score shape (parity with
``scores_prior.parquet``: no ``locus.*``/``alleles`` legacy columns), so it
drops straight into ``create_variant_scores_all_table.py`` once registered in
``data_config/{input_data_locations,score_input_data}.json``::

    chrom               VARCHAR
    pos                 BIGINT
    ref                 VARCHAR
    alt                 VARCHAR
    msa_pairformer_llr  FLOAT   -- canonical UniProt isoform, max over any residual duplicates

Engine / memory strategy
------------------------
The score is scanned once per chromosome via a per-chromosome loop that
INNER-JOINs one linker shard at a time (~500K-3M rows compressed per shard)
against the full score (193.5M rows), applies the canonical-isoform WHERE
filter, and streams the GROUP-BY-reduced result to a per-chrom Parquet
checkpoint under ``.ckpt_<output stem>/``:

    for chrom in chr1..chr22, chrX, chrY:
        linker_shard  =  linker_missense_enst_transcript_aa_uniprot_<n>.tsv.bgz
        checkpoint    =  .ckpt_<output stem>/<chrom>.parquet
        COPY (
            SELECT l.chrom, l.pos, l.ref, l.alt,
                   CAST(max(s.llr) AS FLOAT) AS msa_pairformer_llr
            FROM <linker shard>  l
            INNER JOIN <score>   s
              ON s.uniprot_id = l.uniprot_id AND s.position = l.aa_pos
             AND s.aa_ref     = l.aa_ref     AND s.aa_alt   = l.aa_alt
            WHERE l.uniprot_isoform IS NULL OR l.uniprot_isoform = l.uniprot_id
            GROUP BY l.chrom, l.pos, l.ref, l.alt
        ) TO <checkpoint> (FORMAT PARQUET, ...)

DuckDB builds the hash table on the small linker shard side and streams the
score through it, spilling to ``--temp-dir`` under ``--memory-limit``. Peak
memory is bounded to a single shard + one running aggregation state per
``(chrom, pos, ref, alt)`` group; a crash loses at most one chromosome's worth
of work (the finished checkpoints are picked up on restart).

After every shard has a checkpoint, a single ``COPY (SELECT * FROM
read_parquet([...])) TO <output>`` consolidates them in canonical chromosome
order. Checkpoints are removed on a successful, non-``--keep-intermediates``
run.

.bgz linker reads
-----------------
Per ``merge_v2/.cursor/rules/duckdb-file-reading.mdc``, ``.bgz`` files require
``compression='gzip'`` and ``nullstr=['NA','']`` (the linker uses ``NA`` as a
BOOLEAN null sentinel). An explicit ``types={...}`` map is pinned on every
``read_csv`` call to skip the ~50 s/shard type-inference scan; the map matches
the schema baked into ``msa_pairformer_chr22_sanity_check.py`` so any drift
here or there fails loudly on the first shard.

Downstream wiring
-----------------
Once this script has produced ``msa_pairformer_variant.parquet``, upload it to
GCS (matching how the other single-file scores in ``raw_data/scores/`` are
distributed) and add the URI to::

    merge_v2/data_config/input_data_locations.json      # score_tables.variant_level
    merge_v2/data_config/score_input_data.json          # score_fields projection

At that point the standard pipeline picks it up like any other variant score.

Run from ``merge_v2/preprocess/``::

    python preprocess_msa_pairformer_llr.py
    python preprocess_msa_pairformer_llr.py --dry-run
    python preprocess_msa_pairformer_llr.py --overwrite
    python preprocess_msa_pairformer_llr.py --memory-limit 20GB --threads 4
    python preprocess_msa_pairformer_llr.py --input <score>.parquet --output <path>

Prerequisites (both already fetched by the standard downloaders):

* ``merge_v2/data/raw_data/scores/msa_pairformer_llr_chunks0-9.parquet``
* ``merge_v2/data/raw_data/linker/linker_missense_enst_transcript_aa_uniprot_by_chrom_tsv/
  linker_missense_enst_transcript_aa_uniprot_{1..22,X,Y}.tsv.bgz``
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import duckdb

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: merge_v2/preprocess/...)
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent  # merge_v2/
RAW_SCORES_DIR = PROJECT_DIR / "data" / "raw_data" / "scores"
RAW_LINKER_DIR = (
    PROJECT_DIR / "data" / "raw_data" / "linker"
    / "linker_missense_enst_transcript_aa_uniprot_by_chrom_tsv"
)
DEFAULT_INPUT = RAW_SCORES_DIR / "msa_pairformer_llr_chunks0-9.parquet"
DEFAULT_LINKER_DIR = RAW_LINKER_DIR
DEFAULT_OUTPUT = RAW_SCORES_DIR / "msa_pairformer_variant.parquet"

# Linker shard filenames like ``linker_missense_enst_transcript_aa_uniprot_22.tsv.bgz``.
LINKER_PREFIX = "linker_missense_enst_transcript_aa_uniprot_"
LINKER_SUFFIX = ".tsv.bgz"
_LINKER_NAME_RE = re.compile(
    rf"^{re.escape(LINKER_PREFIX)}([0-9]+|X|Y){re.escape(LINKER_SUFFIX)}$"
)

# Canonical chromosome order (chrY is included but effectively empty in the
# missense linker; visited so its residual rows -- if any -- still make it into
# the consolidated output and the schema stays uniform across checkpoints).
CHROMOSOMES = [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY"]

# Output schema (variant key + one canonical-isoform-filtered score column).
KEY_COLS: list[str] = ["chrom", "pos", "ref", "alt"]
SCORE_COLS: list[str] = ["msa_pairformer_llr"]

# CSV null sentinels for the .bgz linker: ``NA`` for BOOLEAN nulls, ``''`` for
# empty-string nulls. Matches the .cursor/rules recipe and the chr22 sanity
# check.
CSV_NULL_VALUES: list[str] = ["NA", ""]

# Full linker schema pinned by ``read_csv(..., types={...})``. Skipping type
# inference is a ~50 s/shard win; pinning every column (rather than just the
# ones we project) matches the chr22 sanity check so any schema drift trips
# the first shard's read rather than a silent cast surprise mid-loop.
LINKER_TYPES: dict[str, str] = {
    "chrom": "VARCHAR",
    "pos": "BIGINT",
    "ref": "VARCHAR",
    "alt": "VARCHAR",
    "aa_pos": "INTEGER",
    "aa_ref": "VARCHAR",
    "aa_alt": "VARCHAR",
    "gene_symbol": "VARCHAR",
    "enst": "VARCHAR",
    "ensg": "VARCHAR",
    "ensp": "VARCHAR",
    "uniprot_id": "VARCHAR",
    "uniprot_isoform": "VARCHAR",
    "mane_select": "BOOLEAN",
    "canonical": "BOOLEAN",
    "transcript_mane_select": "BOOLEAN",
    "locus": "VARCHAR",
    "alleles": "VARCHAR",
}

PLAN_FILE = "plan.json"
CKPT_PREFIX = ".ckpt_"


# ---------------------------------------------------------------------------
# SQL helpers.
# ---------------------------------------------------------------------------
def q(identifier: str) -> str:
    """Double-quote a SQL identifier (``ref``/``alt`` are DuckDB reserved words)."""
    return '"' + identifier.replace('"', '""') + '"'


def sql_str(value: str) -> str:
    """Single-quote a SQL string literal for DuckDB path interpolation."""
    return "'" + value.replace("'", "''") + "'"


def _linker_types_sql() -> str:
    """The ``types={...}`` map, formatted for a ``read_csv`` call site."""
    return (
        "{"
        + ", ".join(f"{sql_str(k)}: {sql_str(v)}" for k, v in LINKER_TYPES.items())
        + "}"
    )


def _linker_nullstr_sql() -> str:
    """The ``nullstr=[...]`` list, formatted for a ``read_csv`` call site."""
    return "[" + ", ".join(sql_str(v) for v in CSV_NULL_VALUES) + "]"


def linker_reader_sql(shard_path: Path) -> str:
    """``read_csv(...)`` expression for one .bgz linker shard.

    Uses the .bgz recipe from ``merge_v2/.cursor/rules/duckdb-file-reading.mdc``
    (``compression='gzip'``, ``nullstr=['NA','']``, explicit ``types={...}``).
    """
    return (
        "read_csv("
        f"{sql_str(str(shard_path))}, "
        "delim='\\t', header=true, "
        "compression='gzip', "
        f"nullstr={_linker_nullstr_sql()}, "
        f"types={_linker_types_sql()})"
    )


def score_reader_sql(score_path: Path) -> str:
    """``read_parquet(...)`` expression for the uniprot-grain score."""
    return f"read_parquet({sql_str(str(score_path))})"


def coalesce_select_sql(linker_reader: str, score_reader: str) -> str:
    """Per-chromosome ``INNER JOIN`` + canonical-isoform ``WHERE`` filter +
    ``GROUP BY`` reducing to variant grain.

    Emits one row per ``(chrom, pos, ref, alt)`` with a single
    ``max(s.llr)``-reduced ``msa_pairformer_llr`` column, taken only over
    linker rows where the UniProt isoform is either NULL (absent isoform tag
    treated as "canonical by default", matching the sanity check's convention)
    or equal to the UniProt id itself (i.e. the canonical isoform). Rows
    matching non-canonical isoforms are dropped by the ``WHERE`` filter and do
    not contribute to the aggregation.

    The ``max`` is a defensive aggregation: after the canonical-only filter,
    most ``(chrom, pos, ref, alt)`` groups collapse to a single linker row per
    UniProt id, but a group can still see multiple linker rows when the same
    variant maps to canonical isoforms of multiple UniProt entries (rare) or
    when the linker's own de-duplication kept a residual duplicate. Taking
    the max keeps the strongest LLR signal in either case.

    The output column is ``CAST(... AS FLOAT)`` regardless of the input's
    ``llr`` type (which happens to already be FLOAT), so every per-chrom
    checkpoint has an identical schema and the final ``UNION ALL``
    consolidation cannot hit a type mismatch on an empty-vs-populated
    chromosome (chrY, notably, has essentially zero linker rows).
    """
    return (
        "SELECT "
        "l.chrom AS chrom, "
        "l.pos AS pos, "
        "l.ref AS ref, "
        "l.alt AS alt, "
        f"CAST(max(s.llr) AS FLOAT) AS {q('msa_pairformer_llr')} "
        f"FROM {linker_reader} l "
        f"INNER JOIN {score_reader} s "
        "ON s.uniprot_id = l.uniprot_id "
        "AND s.position = l.aa_pos "
        "AND s.aa_ref = l.aa_ref "
        "AND s.aa_alt = l.aa_alt "
        "WHERE l.uniprot_isoform IS NULL OR l.uniprot_isoform = l.uniprot_id "
        "GROUP BY l.chrom, l.pos, l.ref, l.alt"
    )


def copy_sql(
    source_sql: str, out_path: Path, compression: str, row_group_size: int
) -> str:
    """``COPY (<select>) TO '<path>' (FORMAT PARQUET, ...)``."""
    return (
        f"COPY ({source_sql}) TO {sql_str(str(out_path))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(compression)}, "
        f"ROW_GROUP_SIZE {int(row_group_size)})"
    )


def configure(
    con: duckdb.DuckDBPyConnection,
    memory_limit: str | None,
    threads: int | None,
    temp_dir: Path,
) -> None:
    """Apply memory / spill / thread settings; ensure ``temp_dir`` exists."""
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET temp_directory = {sql_str(str(temp_dir))}")
    con.execute("SET preserve_insertion_order = false")
    if memory_limit:
        con.execute(f"SET memory_limit = {sql_str(memory_limit)}")
    if threads:
        con.execute(f"SET threads = {int(threads)}")


# ---------------------------------------------------------------------------
# Linker-shard discovery.
# ---------------------------------------------------------------------------
def discover_linker_shards(directory: Path) -> list[tuple[str, Path]]:
    """``[(chrom_label, path), ...]`` for the linker shards, in chromosome order.

    Files that do not match the ``linker_missense_enst_transcript_aa_uniprot_{N}.tsv.bgz``
    pattern (or that name a chromosome outside ``chr1..chr22, chrX, chrY``)
    are ignored, so extra unrelated files in the directory don't break the
    loop.
    """
    found: dict[str, Path] = {}
    for p in sorted(directory.glob(f"{LINKER_PREFIX}*{LINKER_SUFFIX}")):
        m = _LINKER_NAME_RE.match(p.name)
        if m is None:
            continue
        found["chr" + m.group(1)] = p
    return [(c, found[c]) for c in CHROMOSOMES if c in found]


# ---------------------------------------------------------------------------
# Checkpoint bookkeeping.
# ---------------------------------------------------------------------------
def ckpt_dir_for(output: Path) -> Path:
    return output.parent / f"{CKPT_PREFIX}{output.stem}"


def ckpt_path_for(ckpt_dir: Path, chrom_label: str) -> Path:
    return ckpt_dir / f"{chrom_label}.parquet"


def build_plan_dict(
    args: argparse.Namespace, linker_shards: list[tuple[str, Path]]
) -> dict:
    """Serialisable resume-validation plan; the checkpoint dir is rejected if
    this dictionary no longer matches on a subsequent run."""
    return {
        "input_score": str(args.input),
        "linker_dir": str(args.linker_dir),
        "output": str(args.output),
        "chromosomes": [c for c, _ in linker_shards],
        "linker_shards": {c: p.name for c, p in linker_shards},
        "output_columns": list(KEY_COLS + SCORE_COLS),
    }


def write_plan(ckpt_dir: Path, plan: dict) -> None:
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    (ckpt_dir / PLAN_FILE).write_text(json.dumps(plan, indent=2))


def load_plan(ckpt_dir: Path) -> dict | None:
    p = ckpt_dir / PLAN_FILE
    if not p.exists():
        return None
    return json.loads(p.read_text())


def validate_resume(ckpt_dir: Path, plan: dict) -> None:
    """Refuse to resume against a plan that no longer matches the inputs."""
    existing = load_plan(ckpt_dir)
    if existing is None:
        return
    # JSON round-trips tuples to lists; normalise before comparing.
    expected = json.loads(json.dumps(plan))
    if existing != expected:
        sys.exit(
            "ERROR: existing checkpoints in\n"
            f"  {ckpt_dir}\n"
            "were built with a different plan (input score path, linker dir, "
            "output path, or chromosome list changed). Re-run with "
            "--overwrite to start clean, or restore the original inputs."
        )


# ---------------------------------------------------------------------------
# Per-chromosome runner + consolidation.
# ---------------------------------------------------------------------------
def run_chromosome_batch(
    score_path: Path,
    linker_shard: Path,
    chrom_label: str,
    out_path: Path,
    temp_dir: Path,
    args: argparse.Namespace,
) -> int:
    """Stream one chromosome's linker shard INNER-JOINed with the score,
    reduced to variant grain, into a per-chrom Parquet checkpoint.

    Returns the checkpoint row count (read back from the Parquet footer, so
    exact and cheap).
    """
    shutil.rmtree(temp_dir, ignore_errors=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    linker_reader = linker_reader_sql(linker_shard)
    score_reader = score_reader_sql(score_path)
    select_sql = coalesce_select_sql(linker_reader, score_reader)

    con = duckdb.connect()
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        print(
            f"    streaming {chrom_label} (linker: {linker_shard.name}) -> "
            f"{out_path.name} ...",
            flush=True,
        )
        con.execute(copy_sql(select_sql, out_path, args.compression, args.row_group_size))
        n_rows = con.execute(
            f"SELECT count(*) FROM read_parquet({sql_str(str(out_path))})"
        ).fetchone()[0]
    finally:
        con.close()
    shutil.rmtree(temp_dir, ignore_errors=True)
    return int(n_rows)


def consolidate(
    ckpt_dir: Path,
    chrom_labels: list[str],
    out_path: Path,
    temp_dir: Path,
    args: argparse.Namespace,
) -> int:
    """Concat all per-chrom checkpoints into the final Parquet.

    Reads each ``<ckpt_dir>/<chrom_label>.parquet`` in canonical chromosome
    order via a single ``UNION ALL`` ``read_parquet(list)`` so DuckDB streams
    them sequentially rather than treating them as a partitioned dataset
    (avoids any glob-order surprises).
    """
    shutil.rmtree(temp_dir, ignore_errors=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    parts = [str(ckpt_path_for(ckpt_dir, c)) for c in chrom_labels]
    parts_sql = "[" + ", ".join(sql_str(p) for p in parts) + "]"
    select_sql = f"SELECT * FROM read_parquet({parts_sql})"

    con = duckdb.connect()
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        print(f"  streaming consolidation -> {out_path.name} ...", flush=True)
        con.execute(copy_sql(select_sql, out_path, args.compression, args.row_group_size))
        n_rows = con.execute(
            f"SELECT count(*) FROM read_parquet({sql_str(str(out_path))})"
        ).fetchone()[0]
    finally:
        con.close()
    shutil.rmtree(temp_dir, ignore_errors=True)
    return int(n_rows)


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT,
        help=f"UniProt-grain MSA-Pairformer LLR parquet (default: {DEFAULT_INPUT}).",
    )
    parser.add_argument(
        "--linker-dir", "--linker_dir", dest="linker_dir",
        type=Path, default=DEFAULT_LINKER_DIR,
        help=f"Directory of per-chromosome missense/UniProt linker shards "
             f"(default: {DEFAULT_LINKER_DIR}).",
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help=f"Output variant-grain Parquet path "
             f"(default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Discard any existing output and checkpoints and start clean.",
    )
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Keep the per-chromosome checkpoint directory after a successful "
             "run (for debugging).",
    )
    parser.add_argument(
        "--memory-limit", default=None,
        help="DuckDB memory limit, e.g. '12GB' (default: DuckDB's ~80%% of RAM).",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="DuckDB worker threads (default: DuckDB auto).",
    )
    parser.add_argument(
        "--temp-dir", type=Path, default=None,
        help="Directory for DuckDB spill files during each per-chromosome "
             "batch (default: <output dir>/.duckdb_spill).",
    )
    parser.add_argument(
        "--compression", default="zstd",
        help="Parquet compression codec (default: zstd; matches the pipeline).",
    )
    parser.add_argument(
        "--row-group-size", type=int, default=512_000,
        help="Parquet row group size (default: 512000; matches the pipeline).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the per-chromosome plan (and a representative SQL "
             "statement) and exit without running any joins.",
    )
    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(
            f"ERROR: MSA-Pairformer score not found: {args.input}\n"
            f"       Download it into {RAW_SCORES_DIR} first."
        )
    if not args.linker_dir.exists():
        sys.exit(
            f"ERROR: missense/UniProt linker directory not found: {args.linker_dir}\n"
            f"       Download the shards into it first (they ship as part of "
            f"the standard 'linker' download group)."
        )

    linker_shards = discover_linker_shards(args.linker_dir)
    if not linker_shards:
        sys.exit(
            f"ERROR: no {LINKER_PREFIX}*{LINKER_SUFFIX} files found under\n"
            f"  {args.linker_dir}\n"
            f"       Expected shards matching {LINKER_PREFIX}"
            "{1..22,X,Y}" + LINKER_SUFFIX + "."
        )

    out_dir = args.output.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or (out_dir / ".duckdb_spill")
    ckpt_dir = ckpt_dir_for(args.output)

    print(f"Input score:  {args.input}")
    print(f"Linker dir:   {args.linker_dir}")
    print(f"Output:       {args.output}")
    print(
        f"Linker shards: {len(linker_shards)} "
        f"({', '.join(c for c, _ in linker_shards)})"
    )
    print(f"Score columns to emit ({len(SCORE_COLS)}, FLOAT):")
    for c in SCORE_COLS:
        print(f"  {c}   (canonical UniProt isoform only; max over any residual duplicates)")
    print(f"Join type: INNER on (uniprot_id, aa_pos, aa_ref, aa_alt); "
          "canonical-isoform WHERE filter; "
          f"GROUP BY {KEY_COLS}")
    print(f"Checkpoints: {ckpt_dir}")
    print(f"Spill dir:   {temp_dir}\n")

    if args.dry_run:
        print("Per-chromosome plan:")
        for c, p in linker_shards:
            print(f"  {c:6s} -> {p.name}  ({ckpt_path_for(ckpt_dir, c).name})")
        print("\nRepresentative SQL (first chromosome):")
        first_chrom, first_path = linker_shards[0]
        sql = coalesce_select_sql(
            linker_reader_sql(first_path), score_reader_sql(args.input)
        )
        print(
            copy_sql(
                sql,
                ckpt_path_for(ckpt_dir, first_chrom),
                args.compression,
                args.row_group_size,
            )
            + ";"
        )
        return 0

    if args.overwrite:
        shutil.rmtree(ckpt_dir, ignore_errors=True)
        args.output.unlink(missing_ok=True)

    if args.output.exists():
        print("Output already exists; nothing to do (use --overwrite to rebuild).")
        return 0

    plan = build_plan_dict(args, linker_shards)
    validate_resume(ckpt_dir, plan)
    write_plan(ckpt_dir, plan)

    # Per-chromosome loop. Each iteration writes exactly one checkpoint; if a
    # checkpoint already exists we skip it (that's how we resume). Peak memory
    # per chrom = one shard's hash-build side + the running GROUP-BY state.
    total_rows = 0
    for i, (chrom_label, linker_path) in enumerate(linker_shards, start=1):
        ckpt_path = ckpt_path_for(ckpt_dir, chrom_label)
        header = f"[{i}/{len(linker_shards)}] {chrom_label}"
        if ckpt_path.exists():
            existing = duckdb.connect()
            try:
                n = existing.execute(
                    f"SELECT count(*) FROM read_parquet({sql_str(str(ckpt_path))})"
                ).fetchone()[0]
            finally:
                existing.close()
            total_rows += int(n)
            print(f"  {header}: checkpoint exists ({n:,} rows) -- skipping.")
            continue

        print(f"  {header}: joining {linker_path.name} ...")
        n = run_chromosome_batch(
            args.input, linker_path, chrom_label, ckpt_path, temp_dir, args
        )
        total_rows += n
        print(f"     checkpoint written: {n:,} rows -> {ckpt_path.name}")

    # Consolidation: single COPY reading all per-chrom checkpoints in
    # canonical chromosome order.
    chrom_labels = [c for c, _ in linker_shards]
    final_rows = consolidate(ckpt_dir, chrom_labels, args.output, temp_dir, args)

    # Sanity: the consolidated file should have the same row count as the sum
    # of the per-chromosome checkpoints.
    if final_rows != total_rows:
        print(
            f"WARNING: consolidated row count {final_rows:,} != "
            f"sum of per-chrom checkpoints {total_rows:,}",
            file=sys.stderr,
        )

    if not args.keep_intermediates:
        shutil.rmtree(ckpt_dir, ignore_errors=True)
    shutil.rmtree(temp_dir, ignore_errors=True)

    final_con = duckdb.connect()
    try:
        cols = [
            row[0]
            for row in final_con.execute(
                f"DESCRIBE SELECT * FROM read_parquet({sql_str(str(args.output))})"
            ).fetchall()
        ]
    finally:
        final_con.close()
    print(
        f"\nDone. Wrote {final_rows:,} rows x {len(cols)} columns to "
        f"{args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
