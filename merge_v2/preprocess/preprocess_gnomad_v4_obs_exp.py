#!/usr/bin/env python3
"""Preprocess gnomad_v4_obs_exp_per_variant.parquet into the eval-obs/exp shape.

The raw source at

    gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/obs_exp/gnomad_v4_obs_exp_per_variant.parquet

is a Spark-style partitioned Hail export with a different shape than the
existing per-cohort obs/exp eval sources consumed by
``create_variant_eval_all_table.py``.  Specifically:

- **No transcript-level columns** (``canonical``, ``mane_select``,
  ``annotation``, ``transcript``).  The data is per-(variant, gene), not
  per-transcript, so the ``WHERE`` filters from ``preprocess_gnomad_obs_exp.py``
  don't apply.
- **Dual gene_id encoding** — every gene appears twice (once as ENSG, once as a
  numeric Entrez/NCBI ID), a Hail VEP artefact.  ~51% of rows in the first
  part file carry non-ENSG IDs.
- **Exact duplicate rows** within the ENSG subset (~22 k groups of dupes in the
  first part file alone), likely a residual transcript-level fan-out that wasn't
  collapsed upstream.
- **Generic column names** (``observed``/``expected``) without a cohort suffix.
- **Hail-only keys** (``locus.contig``, ``locus.position``, ``alleles[]``) — no
  split ``chrom``/``pos``/``ref``/``alt``.

This script normalizes the raw source into the same shape as the other
``*_obs_exp_mis.parquet`` eval files:

1. **INNER JOINs with the variant→gene linker** (``linker_all.parquet``) to
   restrict to missense variants only.  Without this, 77% of the raw variant
   keys are non-missense and would bloat downstream tables.
2. Filters to ``gene_id LIKE 'ENSG%'`` (drops the Entrez duplicates).
3. Deduplicates exact-duplicate rows (``SELECT DISTINCT``).
4. Derives ``chrom``, ``pos``, ``ref``, ``alt`` from the Hail columns.
5. Renames ``gene_id`` → ``ensg``, ``observed`` → ``observed_gnomad_v4``,
   ``expected`` → ``expected_gnomad_v4``.
6. Sets ``enst`` to ``NULL`` (column absent in this source but present in
   the canonical schema; the downstream merge ignores it).

Output schema matches the existing obs/exp evals::

    locus.contig, locus.position, alleles, enst, ensg,
    observed_gnomad_v4, expected_gnomad_v4, chrom, pos, ref, alt

Engine / memory strategy
------------------------
Single streaming ``COPY (SELECT DISTINCT ... FROM read_parquet(...) WHERE ...)
TO ...``: DuckDB scans the 4 340 Spark part files, applies the ENSG filter +
``DISTINCT``, and writes zstd-compressed Parquet with the pipeline's standard
512 k row-group size.  Peak memory is bounded by one scan + the dedup hash
table; nothing else is materialised in RAM.

Usage::

    # From merge_v2/preprocess/
    python preprocess_gnomad_v4_obs_exp.py
    python preprocess_gnomad_v4_obs_exp.py --dry-run
    python preprocess_gnomad_v4_obs_exp.py --overwrite --memory-limit 8GB
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import duckdb

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent  # merge_v2/
DEFAULT_INPUT = PROJECT_DIR / "data" / "raw_data" / "evals" / "gnomad_v4_obs_exp_per_variant.parquet"
DEFAULT_OUTPUT = PROJECT_DIR / "data" / "raw_data" / "evals" / "gnomad_v4_obs_exp_mis.parquet"
DEFAULT_LINKER = PROJECT_DIR / "data" / "raw_data" / "linker" / "linker_all.parquet"

COHORT = "gnomad_v4"


def q(identifier: str) -> str:
    return '"' + identifier.replace('"', '""') + '"'


def sql_str(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def reader_sql(resolved: Path) -> str:
    """``read_parquet(...)`` for a Spark-partitioned parquet directory.

    Uses ``part-*.parquet`` to skip ``_SUCCESS`` and any stale
    ``_temporary/`` task-attempt artefacts.
    """
    pattern = str(resolved / "part-*.parquet") if resolved.is_dir() else str(resolved)
    return f"read_parquet({sql_str(pattern)})"


def preprocess_sql(reader: str, linker_path: Path) -> str:
    obs_out = f"observed_{COHORT}"
    exp_out = f"expected_{COHORT}"
    linker = f"read_parquet({sql_str(str(linker_path))})"
    return (
        "SELECT DISTINCT "
        'r."locus.contig" AS "locus.contig", '
        'r."locus.position" AS "locus.position", '
        "r.alleles AS alleles, "
        "NULL::VARCHAR AS enst, "
        "r.gene_id AS ensg, "
        f"r.observed AS {q(obs_out)}, "
        f"r.expected AS {q(exp_out)}, "
        'r."locus.contig" AS chrom, '
        'r."locus.position" AS pos, '
        "r.alleles[1] AS ref, "
        "r.alleles[2] AS alt "
        f"FROM {reader} AS r "
        f"INNER JOIN (SELECT DISTINCT chrom, pos, ref, alt FROM {linker}) AS lk "
        'ON r."locus.contig" = lk.chrom '
        'AND r."locus.position" = lk.pos '
        "AND r.alleles[1] = lk.ref "
        "AND r.alleles[2] = lk.alt "
        "WHERE r.gene_id LIKE 'ENSG%'"
    )


def copy_sql(source_sql: str, out_path: Path, compression: str, row_group_size: int) -> str:
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
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET temp_directory = {sql_str(str(temp_dir))}")
    con.execute("SET preserve_insertion_order = false")
    if memory_limit:
        con.execute(f"SET memory_limit = {sql_str(memory_limit)}")
    if threads:
        con.execute(f"SET threads = {int(threads)}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT,
        help=f"Raw Spark-partitioned parquet directory (default: {DEFAULT_INPUT}).",
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help=f"Output parquet path (default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--linker", type=Path, default=DEFAULT_LINKER,
        help=f"Variant→gene linker parquet (default: {DEFAULT_LINKER}). "
             "Used to INNER JOIN, restricting output to missense variants.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--memory-limit", default=None)
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--temp-dir", type=Path, default=None)
    parser.add_argument("--compression", default="zstd")
    parser.add_argument("--row-group-size", type=int, default=512_000)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    input_path = args.input.resolve()
    out_path = args.output.resolve()

    linker_path = args.linker.resolve()

    if not input_path.exists():
        sys.exit(
            f"ERROR: raw input not found: {input_path}\n"
            f"       Download it first, e.g.:\n"
            f"         gcloud storage cp -r "
            f"gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/obs_exp/"
            f"{input_path.name} {input_path.parent}/"
        )
    if not linker_path.exists():
        sys.exit(f"ERROR: linker not found: {linker_path}")

    reader = reader_sql(input_path)
    source_sql = preprocess_sql(reader, linker_path)
    copy = copy_sql(source_sql, out_path, args.compression, args.row_group_size)

    print(f"=== gnomad_v4_obs_exp_per_variant -> {out_path.name} ===")
    print(f"Cohort tag: {COHORT}")
    print(f"Signal columns: observed_{COHORT}, expected_{COHORT}")
    print(f"Input:  {input_path}")
    print(f"Output: {out_path}")

    if args.dry_run:
        print("  (dry run -- nothing written)")
        print(copy + ";")
        return 0

    if out_path.exists() and not args.overwrite:
        print("  skip (output exists; use --overwrite to rebuild)")
        return 0

    temp_dir = args.temp_dir or (out_path.parent / ".duckdb_spill")
    con = duckdb.connect()
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)

        n_raw = con.execute(f"SELECT count(*) FROM {reader}").fetchone()[0]
        print(f"  raw rows: {n_raw:,}")

        print("  filtering ENSG + dedup + splitting alleles + writing parquet ...",
              flush=True)
        con.execute(copy)

        n_out = con.execute(
            f"SELECT count(*) FROM read_parquet({sql_str(str(out_path))})"
        ).fetchone()[0]
        print(f"  wrote {n_out:,} rows  ({n_raw - n_out:,} dropped by ENSG filter + dedup)")

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
                f"multiple overlapping genes."
            )
    finally:
        con.close()
        shutil.rmtree(temp_dir, ignore_errors=True)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
