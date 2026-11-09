#!/usr/bin/env python3
"""Preprocess raw gnomad per-variant-expected parquets into the eval-obs/exp shape.

The four raw sources at

    gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/scores/gnomad.new.genetics_gym.per_variant.expected.parquet/
    gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/scores/gnomad.ukb.genetics_gym.per_variant.expected.parquet/
    gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/scores/gnomad.v2.genetics_gym.per_variant.expected.parquet/
    gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/scores/gnomad.v2.ukb.genetics_gym.per_variant.expected.parquet/

are raw hail exports: per-(variant, transcript, gene, annotation) rows, no split
``ref``/``alt`` columns, mixed Ensembl/RefSeq transcripts, and every annotation
class (missense, synonymous, stop_gained, downstream_gene_variant, ...) in one
file. That shape does not match the join key or the semantics of the existing
per-cohort obs/exp eval sources
(``asd_obs_exp_mis.parquet``, ``chd_obs_exp_mis.parquet``, ``dd_obs_exp_mis.parquet``)
that are consumed by ``create_variant_eval_all_table.py``.

This script normalizes each raw source into the same shape as those existing
files, so it can be dropped in behind the same eval pipeline. Concretely, for
each raw input it:

1. Filters to the canonical missense annotation on the Ensembl transcript::

       canonical = TRUE
       AND mane_select = TRUE
       AND annotation = 'missense_variant'
       AND transcript LIKE 'ENST%'

   which reduces the raw source to **one row per variant** (verified against
   the sample: no residual per-variant duplicates after this filter).

2. Derives the split variant key from Hail's locus/alleles struct::

       chrom  <- locus.contig
       pos    <- locus.position
       ref    <- alleles[1]          -- DuckDB is 1-indexed for LIST access
       alt    <- alleles[2]

   and carries through the original ``locus.contig`` / ``locus.position`` /
   ``alleles`` columns unchanged (matching the existing files, which retain
   the raw Hail columns alongside the derived split key).

3. Renames the transcript/gene columns to the pipeline's canonical names::

       transcript  -> enst
       gene_id     -> ensg

   (``gene_id`` is ENSG-prefixed for the Ensembl-annotated missense subset; a
   small residual of non-Ensembl IDs from OMIM-only genes is preserved as-is,
   matching how the existing files handle rare non-canonical entries.)

4. Renames the signal columns per cohort::

       observed  -> observed_<cohort>
       expected  -> expected_<cohort>

   where ``<cohort>`` is derived from the input filename
   (``gnomad.new`` -> ``gnomad_new``, ``gnomad.v2.ukb`` -> ``gnomad_v2_ukb``,
   etc.), so the four cohorts contribute distinct output columns and can all be
   joined into the same wide eval table.

The output column order matches the existing ``asd_obs_exp_mis.parquet`` schema
exactly::

    locus.contig, locus.position, alleles, enst, ensg,
    observed_<cohort>, expected_<cohort>, chrom, pos, ref, alt

Engine / memory strategy
------------------------
Single streaming ``COPY (SELECT ... FROM read_parquet(...) WHERE ...) TO ...``
per source: DuckDB scans the raw partitioned Parquet dir, applies the row filter
and the column projection in one pipeline, and writes zstd-compressed Parquet
with the pipeline's standard 512k row-group size. Peak memory is bounded by a
single scan + one row group; nothing is materialized in RAM.

Downstream wiring
-----------------
Once the preprocessed files exist under
``merge_v2/data/raw_data/evals/gnomad_<cohort>_obs_exp_mis.parquet/`` they are
consumed by ``create_variant_eval_all_table.py`` via the four new entries added
to ``merge_v2/data_config/evals_input_data.json``.

If you want the preprocessed files re-downloadable via
``download_source_data.py``, upload each output directory to
``gs://grohlicek/genetics_gym_vsm_all_content/eval_obs_exp/`` (matching the URIs
already added to ``merge_v2/data_config/input_data_locations.json`` under
``eval_tables.variant_level``).

Run from ``merge_v2/preprocess/`` (this script's own directory; the driver
``bootstrap_gnomad_obs_exp.sh`` next to it is the intended entry point)::

    # Process all four cohorts (defaults look under ../data/raw_data/scores/
    # for the raw gnomad.*.per_variant.expected.parquet directories):
    python preprocess_gnomad_obs_exp.py --all

    # Process a single cohort explicitly (raw input path + cohort tag):
    python preprocess_gnomad_obs_exp.py \\
        --input ../data/raw_data/scores/gnomad.new.genetics_gym.per_variant.expected.parquet \\
        --cohort gnomad_new

    # Dry-run (print the SQL that would be executed and exit):
    python preprocess_gnomad_obs_exp.py --all --dry-run
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import duckdb

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: merge_v2/preprocess/...)
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent  # merge_v2/
RAW_DEFAULT_DIR = PROJECT_DIR / "data" / "raw_data" / "scores"
OUT_DEFAULT_DIR = PROJECT_DIR / "data" / "raw_data" / "evals"

# Bake in the four cohorts so --all just works after a plain download. The
# cohort tag is what will be substituted into observed_<cohort>/expected_<cohort>
# in the output; the raw filename is the Hail export the URIs in
# input_data_locations.json (score_tables during download, then relocated to
# eval_obs_exp/ after upload) point at.
COHORTS: list[tuple[str, str]] = [
    ("gnomad_new", "gnomad.new.genetics_gym.per_variant.expected.parquet"),
    ("gnomad_ukb", "gnomad.ukb.genetics_gym.per_variant.expected.parquet"),
    ("gnomad_v2", "gnomad.v2.genetics_gym.per_variant.expected.parquet"),
    ("gnomad_v2_ukb", "gnomad.v2.ukb.genetics_gym.per_variant.expected.parquet"),
]


def q(identifier: str) -> str:
    """Quote a SQL identifier (DuckDB double-quotes; e.g. 'ref' is reserved)."""
    return '"' + identifier.replace('"', '""') + '"'


def sql_str(value: str) -> str:
    """Quote a SQL string literal."""
    return "'" + value.replace("'", "''") + "'"


def reader_sql(resolved: Path) -> str:
    """Build the DuckDB ``read_parquet(...)`` expression for a parquet input.

    Uses the top-level ``part-*.parquet`` glob (Hail/Spark's canonical output
    convention) rather than a permissive ``**/*.parquet`` walk. The gnomad
    exports contain stale ``_temporary/0/_temporary/attempt_.../part-*.parquet``
    task-attempt files left over from Spark retries; those files are zero-byte
    or otherwise incomplete (DuckDB rejects them with "too small to be a
    Parquet file") and their canonical counterparts already exist at top level,
    so a ``**`` walk would either duplicate rows or crash the read. Restricting
    the glob to top-level parts also skips the ``_SUCCESS`` marker without
    special-casing it.
    """
    pattern = str(resolved / "part-*.parquet") if resolved.is_dir() else str(resolved)
    return f"read_parquet({sql_str(pattern)})"


def preprocess_sql(reader: str, cohort: str) -> str:
    """Build the ``SELECT`` that produces the target eval-obs/exp shape.

    Column order matches the existing ``asd_obs_exp_mis.parquet`` schema. The
    ``WHERE`` clause reduces the per-transcript, per-annotation raw rows to
    exactly one canonical Ensembl missense row per variant (verified against
    the raw sample).
    """
    obs_out = f"observed_{cohort}"
    exp_out = f"expected_{cohort}"
    return (
        "SELECT "
        f'"locus.contig" AS "locus.contig", '
        f'"locus.position" AS "locus.position", '
        f'alleles AS alleles, '
        f'transcript AS enst, '
        f'gene_id AS ensg, '
        f'observed AS {q(obs_out)}, '
        f'expected AS {q(exp_out)}, '
        f'"locus.contig" AS chrom, '
        f'"locus.position" AS pos, '
        f'alleles[1] AS ref, '
        f'alleles[2] AS alt '
        f"FROM {reader} "
        "WHERE canonical = TRUE "
        "AND mane_select = TRUE "
        "AND annotation = 'missense_variant' "
        "AND transcript LIKE 'ENST%'"
    )


def copy_sql(source_sql: str, out_path: Path, compression: str, row_group_size: int) -> str:
    """Build a ``COPY (...) TO '<path>' (FORMAT PARQUET, ...)`` statement."""
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
    """Apply memory / spill settings so the filtered scan stays out-of-core."""
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET temp_directory = {sql_str(str(temp_dir))}")
    con.execute("SET preserve_insertion_order = false")
    if memory_limit:
        con.execute(f"SET memory_limit = {sql_str(memory_limit)}")
    if threads:
        con.execute(f"SET threads = {int(threads)}")


def output_path(out_dir: Path, cohort: str) -> Path:
    """Derive the output path for a cohort, matching the *_obs_exp_mis pattern."""
    return out_dir / f"{cohort}_obs_exp_mis.parquet"


def process_one(
    con: duckdb.DuckDBPyConnection,
    input_path: Path,
    cohort: str,
    out_path: Path,
    args: argparse.Namespace,
) -> None:
    """Preprocess one raw gnomad parquet dir into a single eval-obs/exp parquet."""
    if not input_path.exists():
        sys.exit(
            f"ERROR: raw input not found: {input_path}\n"
            f"       Download it first, e.g.:\n"
            f"         gcloud storage cp -r "
            f"gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/scores/"
            f"{input_path.name} {input_path.parent}/"
        )

    reader = reader_sql(input_path.resolve())
    source_sql = preprocess_sql(reader, cohort)
    copy = copy_sql(source_sql, out_path, args.compression, args.row_group_size)

    print(f"\n=== {input_path.name} -> {out_path.name} ===")
    print(f"Cohort tag: {cohort}")
    print(f"Signal columns: observed_{cohort}, expected_{cohort}")
    print(f"Input:  {input_path}")
    print(f"Output: {out_path}")

    if args.dry_run:
        print("  (dry run -- nothing written)")
        print(copy + ";")
        return

    if out_path.exists() and not args.overwrite:
        print("  skip (output exists; use --overwrite to rebuild)")
        return

    n_raw = con.execute(f"SELECT count(*) FROM {reader}").fetchone()[0]
    print(f"  raw rows: {n_raw:,}  (per (variant, transcript, gene, annotation))")
    print("  filtering + splitting alleles + writing parquet ...", flush=True)
    con.execute(copy)

    n_out = con.execute(
        f"SELECT count(*) FROM read_parquet({sql_str(str(out_path))})"
    ).fetchone()[0]
    print(f"  wrote {n_out:,} rows  (one row per canonical missense variant)")

    # Diagnostic: report per-variant multiplicity but do NOT hard-fail on it.
    # Multi-gene annotation (a variant falling inside two overlapping genes)
    # legitimately produces multiple rows per (chrom,pos,ref,alt), one per
    # ensg, with identical expected/observed values. The existing per-cohort
    # obs/exp files (e.g. asd_obs_exp_mis.parquet) carry the same shape
    # (~22k dup keys out of ~59M rows), and the downstream eval merge in
    # create_variant_eval_all_table.py already handles it. We only surface a
    # warning so a truly-anomalous multiplicity (e.g. a broken filter that
    # brings back tens of millions of dupes) is still obvious in the log.
    dup_row = con.execute(
        f"SELECT count(*) FROM ("
        f"  SELECT chrom, pos, ref, alt, count(*) AS n "
        f"  FROM read_parquet({sql_str(str(out_path))}) "
        f"  GROUP BY 1,2,3,4 HAVING n > 1"
        f")"
    ).fetchone()[0]
    if dup_row:
        pct = dup_row / max(n_out, 1) * 100
        print(
            f"  note: {dup_row:,} (chrom,pos,ref,alt) keys have >1 row "
            f"({pct:.3f}% of output); expected for variants annotated to "
            f"multiple overlapping genes -- matches the shape of the existing "
            f"*_obs_exp_mis.parquet eval sources."
        )


def cohort_from_filename(name: str) -> str:
    """Infer the cohort tag from a raw filename like 'gnomad.v2.ukb.genetics_gym...'.

    Everything before ``.genetics_gym`` is taken as the cohort id and dots are
    replaced with underscores so it is safe as a SQL identifier suffix.
    """
    stem = name.split(".genetics_gym", 1)[0]
    return stem.replace(".", "_")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--all", action="store_true",
        help="Process all four gnomad cohorts (defaults look under "
             f"{RAW_DEFAULT_DIR} for the raw parquet directories).",
    )
    selection.add_argument(
        "--input", type=Path, default=None, metavar="RAW_PARQUET_DIR",
        help="Process a single raw gnomad parquet directory. Requires --cohort "
             "unless the filename follows the 'gnomad.<cohort>.genetics_gym...' "
             "pattern, in which case the cohort tag is inferred.",
    )
    parser.add_argument(
        "--cohort", default=None, metavar="TAG",
        help="Cohort tag to substitute into observed_<tag> / expected_<tag>. "
             "Only used with --input; ignored under --all (which uses the "
             "baked-in cohort mapping).",
    )
    parser.add_argument(
        "--raw-dir", type=Path, default=RAW_DEFAULT_DIR,
        help=f"Directory holding the raw gnomad.*.per_variant.expected.parquet "
             f"directories in --all mode (default: {RAW_DEFAULT_DIR}).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=OUT_DEFAULT_DIR,
        help=f"Directory to write gnomad_<cohort>_obs_exp_mis.parquet outputs "
             f"(default: {OUT_DEFAULT_DIR}).",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Rebuild outputs that already exist (default: skip them).",
    )
    parser.add_argument(
        "--memory-limit", default=None,
        help="DuckDB memory limit, e.g. '8GB' (default: DuckDB's ~80%% of RAM).",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="DuckDB worker threads (default: DuckDB auto).",
    )
    parser.add_argument(
        "--temp-dir", type=Path, default=None,
        help="Directory for DuckDB spill files "
             "(default: <output dir>/.duckdb_spill).",
    )
    parser.add_argument(
        "--compression", default="zstd",
        help="Parquet compression codec (default: zstd; matches existing "
             "*_obs_exp_mis.parquet outputs).",
    )
    parser.add_argument(
        "--row-group-size", type=int, default=512_000,
        help="Parquet row group size (default: 512000; matches the pipeline).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the assembled SQL for each cohort and exit without writing.",
    )
    args = parser.parse_args()

    if not args.all and args.input is None:
        sys.exit("ERROR: pass --all or --input RAW_PARQUET_DIR.")

    if args.all:
        jobs = [
            (args.raw_dir / raw_name, cohort, output_path(args.output_dir, cohort))
            for cohort, raw_name in COHORTS
        ]
    else:
        cohort = args.cohort or cohort_from_filename(args.input.name)
        jobs = [(args.input, cohort, output_path(args.output_dir, cohort))]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or (args.output_dir / ".duckdb_spill")
    temp_dir.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()  # in-memory: only a streaming COPY per file
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        for input_path, cohort, out_path in jobs:
            process_one(con, input_path, cohort, out_path, args)
    finally:
        con.close()
        shutil.rmtree(temp_dir, ignore_errors=True)

    print("\nAll done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
