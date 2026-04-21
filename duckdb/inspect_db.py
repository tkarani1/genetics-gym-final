#!/usr/bin/env python3
"""List all registered score and eval tables and report per-table statistics."""
from __future__ import annotations

import argparse
import sys

import duckdb


DEFAULT_DB_PATH = "scores.duckdb"


def _print_sample(con: duckdb.DuckDBPyConnection, quoted: str,
                   sample_rows: int | None) -> None:
    """Print randomly sampled rows from a table."""
    if sample_rows is None or sample_rows <= 0:
        return
    sample = con.execute(
        f"SELECT * FROM {quoted} USING SAMPLE {sample_rows}"
    ).fetchall()
    columns = [
        desc[0]
        for desc in con.execute(
            f"SELECT * FROM {quoted} LIMIT 0"
        ).description
    ]
    print(f"  {'  '.join(f'{c:>14}' for c in columns)}")
    for row in sample:
        print(f"  {'  '.join(f'{v!s:>14}' for v in row)}")
    print()


def _print_score_row(con: duckdb.DuckDBPyConnection, existing_tables: set,
                     table_name: str, source_column: str, source_path: str,
                     deduped: bool, fmt, sample_rows: int | None) -> None:
    if table_name not in existing_tables:
        print(f"{table_name:<25} {'MISSING TABLE':}")
        return
    quoted = f'"{table_name}"'
    stats = con.execute(f"""
        SELECT
            COUNT(*)                    AS total,
            COUNT(score)                AS scored,
            COUNT(*) - COUNT(score)     AS nulls,
            COUNT(DISTINCT "key")       AS unique_keys,
            MIN(score)                  AS min_score,
            MAX(score)                  AS max_score,
            AVG(score)                  AS mean_score,
            MEDIAN(score)               AS median_score
        FROM {quoted}
    """).fetchone()
    total, scored, nulls, unique_keys, mn, mx, mean, median = stats
    dedup_flag = "yes" if deduped else "no"
    print(
        f"{table_name:<25} {source_column:<25} {dedup_flag:>8} {total:>12,} "
        f"{scored:>12,} {nulls:>12,} {unique_keys:>12,} "
        f"{fmt(mn)} {fmt(mx)} {fmt(mean)} {fmt(median)}"
    )
    print(f"  path: {source_path}")
    if unique_keys < total:
        print(f"  *** {total - unique_keys:,} duplicate key(s) detected ***")
    _print_sample(con, quoted, sample_rows)


def _print_eval_row(con: duckdb.DuckDBPyConnection, existing_tables: set,
                    table_name: str, source_column: str, source_path: str,
                    deduped: bool, sample_rows: int | None) -> None:
    if table_name not in existing_tables:
        print(f"{table_name:<25} {'MISSING TABLE':}")
        return
    quoted = f'"{table_name}"'
    stats = con.execute(f"""
        SELECT
            COUNT(*)                                        AS total,
            COUNT(*) FILTER (WHERE is_pos = TRUE)           AS positives,
            COUNT(*) FILTER (WHERE is_pos = FALSE)          AS negatives,
            COUNT(*) FILTER (WHERE is_pos IS NULL)          AS nulls,
            COUNT(DISTINCT "key")                           AS unique_keys
        FROM {quoted}
    """).fetchone()
    total, positives, negatives, nulls, unique_keys = stats
    dedup_flag = "yes" if deduped else "no"
    print(
        f"{table_name:<25} {source_column:<25} {dedup_flag:>8} {total:>12,} "
        f"{positives:>12,} {negatives:>12,} {nulls:>12,} "
        f"{unique_keys:>12,}"
    )
    print(f"  path: {source_path}")
    if unique_keys < total:
        print(f"  *** {total - unique_keys:,} duplicate key(s) detected ***")
    _print_sample(con, quoted, sample_rows)


def _print_merged_row(con: duckdb.DuckDBPyConnection, existing_tables: set,
                      table_name: str, source_column: str, source_path: str,
                      sample_rows: int | None) -> None:
    if table_name not in existing_tables:
        print(f"{table_name:<25} {'MISSING TABLE':}")
        return
    quoted = f'"{table_name}"'
    row_count = con.execute(
        f"SELECT COUNT(*) FROM {quoted}"
    ).fetchone()[0]
    columns = [
        desc[0]
        for desc in con.execute(f"SELECT * FROM {quoted} LIMIT 0").description
    ]
    key_names = {"chrom", "pos", "ref", "alt", "ensg", "key"}
    data_cols = [c for c in columns if c not in key_names]
    print(
        f"{table_name:<25} {source_path:<15} {row_count:>12,} "
        f"{len(columns):>8} {len(data_cols):>12}"
    )
    print(f"  sources: {source_column}")
    print(f"  columns: {', '.join(data_cols)}")
    _print_sample(con, quoted, sample_rows)


def inspect_db(db_path: str, sample_rows: int | None = None) -> None:
    """Print metadata entries and summary statistics for every table.

    For each row in the *metadata* table, verifies the corresponding
    DuckDB table exists and reports row count, null count, unique key
    count, and min/max/mean/median of non-null scores.

    Parameters
    ----------
    sample_rows : int or None
        If set, print this many randomly sampled rows per table.
    """
    con = duckdb.connect(db_path, read_only=True)
    try:
        existing_tables = {
            r[0]
            for r in con.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema = 'main'"
            ).fetchall()
        }

        if "metadata" not in existing_tables:
            print("Database has no metadata table.", file=sys.stderr)
            return

        rows = con.execute(
            "SELECT source_column, source_path, table_name, table_type, deduped "
            "FROM metadata ORDER BY table_type, table_name"
        ).fetchall()

        if not rows:
            print("No tables registered.", file=sys.stderr)
            return

        def fmt(v: float | None) -> str:
            return f"{v:>12.6g}" if v is not None else f"{'N/A':>12}"

        score_rows = [r for r in rows if r[3] == "score"]
        eval_rows = [r for r in rows if r[3] == "eval"]
        merged_rows = [r for r in rows if r[3] == "merged_scores"]

        if score_rows:
            print("SCORE TABLES")
            print(f"{'Table':<25} {'Column':<25} {'Deduped':>8} {'Rows':>12} "
                  f"{'Scored':>12} {'Nulls':>12} {'Unique Keys':>12} "
                  f"{'Min':>12} {'Max':>12} {'Mean':>12} {'Median':>12}")
            print("-" * 170)
            for source_column, source_path, table_name, _, deduped in score_rows:
                _print_score_row(con, existing_tables, table_name,
                                 source_column, source_path, deduped, fmt,
                                 sample_rows)
            print()

        if eval_rows:
            print("EVAL TABLES")
            print(f"{'Table':<25} {'Column':<25} {'Deduped':>8} {'Rows':>12} "
                  f"{'Positive':>12} {'Negative':>12} {'Nulls':>12} "
                  f"{'Unique Keys':>12}")
            print("-" * 132)
            for source_column, source_path, table_name, _, deduped in eval_rows:
                _print_eval_row(con, existing_tables, table_name,
                                source_column, source_path, deduped,
                                sample_rows)
            print()

        if merged_rows:
            print("MERGED TABLES")
            print(f"{'Table':<25} {'Operation':<15} {'Rows':>12} "
                  f"{'Columns':>8} {'Data Cols':>12}")
            print("-" * 75)
            for source_column, source_path, table_name, _, _ in merged_rows:
                _print_merged_row(con, existing_tables, table_name,
                                  source_column, source_path,
                                  sample_rows)
            print()

    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inspect a DuckDB scores database: list tables and statistics.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--sample", type=int, default=None, metavar="N",
        help="Print N randomly sampled rows per table.",
    )
    args = parser.parse_args()
    inspect_db(args.db, sample_rows=args.sample)


if __name__ == "__main__":
    main()
