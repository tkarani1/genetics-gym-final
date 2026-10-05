#!/usr/bin/env python3
"""Normalize the gnomAD v2.1.1 MPC-liftover export into variant-key shape.

The raw source (converted from the native Hail table
``gs://grohlicek/genetics_gym_vsm_all_content/from_ruchit/gnomad_v2.1.1_mpc_liftover_GRCh38.ht/``
by hand, since it wasn't already Hail-exported like this pipeline's other raw
score sources) is a Spark-partitioned Parquet directory at

    gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/scores/gnomad_v2_mpc_liftover_GRCh38.parquet

with schema (one row per lifted-over allele)::

    locus.contig: string           -- GRCh38 contig, e.g. 'chr1' (also carries
                                       ALT-contig and 'chrY' rows -- see below)
    locus.position: int32
    alleles: list<string>          -- always length 2: [ref, alt]
    transcript_grch37: string      -- MPC transcript context (not needed downstream)
    original_locus.contig: string  -- pre-liftover GRCh37 locus (not needed downstream)
    original_locus.position: int32
    original_alleles: list<string>
    ref_allele_mismatch: bool      -- liftover QC flag
    mpc_liftover: double           -- the score itself

This script:

1. Drops ``ref_allele_mismatch = True`` rows (bad/ambiguous liftover; 1,412 of
   67,996,604 rows in the current source).
2. Restricts ``locus.contig`` to the same chromosome set already used
   throughout this pipeline (``chr1``..``chr22``, ``chrX``) -- the source also
   carries a handful of ALT-contig rows (e.g. ``chr15_KI270850v1_alt``) and
   ``chrY`` rows that cannot join to anything else in the merge and would just
   be dead weight in the outer table.
3. Splits ``alleles`` into ``ref``/``alt``, renames ``locus.contig`` /
   ``locus.position`` to ``chrom``/``pos``, and keeps ``mpc_liftover``.
4. Drops the transcript/original-locus/QC columns -- not needed downstream
   once the QC filter above has been applied.

Output schema: ``chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR,
mpc_liftover FLOAT``.

Usage::

    python preprocess_mpc_liftover.py --input /path/to/raw/gnomad_v2_mpc_liftover_GRCh38.parquet
    python preprocess_mpc_liftover.py --input ... --dry-run
    python preprocess_mpc_liftover.py --input ... --overwrite
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import duckdb

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent  # merge_v2/
DEFAULT_OUTPUT = (
    PROJECT_DIR / "data" / "raw_data" / "scores" / "gnomad_v2_mpc_liftover_GRCh38_variant.parquet"
)

STANDARD_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def q(identifier: str) -> str:
    return '"' + identifier.replace('"', '""') + '"'


def sql_str(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def reader_sql(resolved: Path) -> str:
    """``read_parquet(...)`` for the Spark-partitioned raw directory."""
    pattern = str(resolved / "part-*.parquet") if resolved.is_dir() else str(resolved)
    return f"read_parquet({sql_str(pattern)})"


def preprocess_sql(reader: str) -> str:
    chroms = ", ".join(sql_str(c) for c in STANDARD_CHROMS)
    return (
        "SELECT "
        'r."locus.contig" AS chrom, '
        'r."locus.position" AS pos, '
        "r.alleles[1] AS ref, "
        "r.alleles[2] AS alt, "
        "TRY_CAST(r.mpc_liftover AS FLOAT) AS mpc_liftover "
        f"FROM {reader} AS r "
        "WHERE (r.ref_allele_mismatch IS NULL OR r.ref_allele_mismatch = FALSE) "
        f'AND r."locus.contig" IN ({chroms})'
    )


def copy_sql(source_sql: str, out_path: Path, compression: str, row_group_size: int) -> str:
    return (
        f"COPY ({source_sql}) TO {sql_str(str(out_path))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(compression)}, "
        f"ROW_GROUP_SIZE {int(row_group_size)})"
    )


def configure(con, memory_limit, threads, temp_dir: Path) -> None:
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET temp_directory = {sql_str(str(temp_dir))}")
    con.execute("SET preserve_insertion_order = false")
    if memory_limit:
        con.execute(f"SET memory_limit = {sql_str(memory_limit)}")
    if threads:
        con.execute(f"SET threads = {int(threads)}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", required=True, type=Path,
        help="Raw Spark-partitioned parquet directory (downloaded from GCS).",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
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

    if not input_path.exists():
        sys.exit(f"ERROR: raw input not found: {input_path}")

    reader = reader_sql(input_path)
    source_sql = preprocess_sql(reader)
    copy = copy_sql(source_sql, out_path, args.compression, args.row_group_size)

    print(f"=== gnomad_v2_mpc_liftover_GRCh38 (raw) -> {out_path.name} ===")
    print(f"Input:  {input_path}")
    print(f"Output: {out_path}")

    if args.dry_run:
        print("  (dry run -- nothing written)")
        print(copy + ";")
        return 0

    if out_path.exists() and not args.overwrite:
        print("  skip (output exists; use --overwrite to rebuild)")
        return 0

    out_path.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or (out_path.parent / ".duckdb_spill")
    con = duckdb.connect()
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)

        n_raw = con.execute(f"SELECT count(*) FROM {reader}").fetchone()[0]
        n_mismatch = con.execute(
            f"SELECT count(*) FROM {reader} WHERE ref_allele_mismatch = TRUE"
        ).fetchone()[0]
        n_nonstandard_chrom = con.execute(
            f'SELECT count(*) FROM {reader} WHERE "locus.contig" NOT IN '
            f"({', '.join(sql_str(c) for c in STANDARD_CHROMS)})"
        ).fetchone()[0]
        print(f"  raw rows: {n_raw:,}")
        print(f"  ref_allele_mismatch=True (dropped): {n_mismatch:,}")
        print(f"  non-standard chrom, e.g. ALT contigs/chrY (dropped): {n_nonstandard_chrom:,}")

        print("  filtering + splitting alleles + writing parquet ...", flush=True)
        con.execute(copy)

        n_out = con.execute(
            f"SELECT count(*) FROM read_parquet({sql_str(str(out_path))})"
        ).fetchone()[0]
        print(f"  wrote {n_out:,} rows ({n_raw - n_out:,} dropped total)")

        n_null_score = con.execute(
            f"SELECT count(*) FROM read_parquet({sql_str(str(out_path))}) WHERE mpc_liftover IS NULL"
        ).fetchone()[0]
        print(f"  null mpc_liftover: {n_null_score:,} ({n_null_score / max(n_out, 1) * 100:.3f}%)")

        dup_row = con.execute(
            f"SELECT count(*) FROM ("
            f"  SELECT chrom, pos, ref, alt, count(*) AS n "
            f"  FROM read_parquet({sql_str(str(out_path))}) "
            f"  GROUP BY 1,2,3,4 HAVING n > 1"
            f")"
        ).fetchone()[0]
        print(f"  duplicate (chrom,pos,ref,alt) keys: {dup_row:,} (handled by the merge step's own dedup)")
    finally:
        con.close()
        shutil.rmtree(temp_dir, ignore_errors=True)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
