#!/usr/bin/env python3
"""Ingest a single score column from a parquet file into a DuckDB table."""
from __future__ import annotations

import argparse
import os
import sys

import duckdb


DEFAULT_DB_PATH = "scores.duckdb"


def ingest_score(
    db_path: str,
    score_name: str,
    score_path: str,
    table_name: str,
    score_type: str,
) -> None:
    """Read one score column from a parquet file and materialise it as a
    DuckDB table with canonical key columns, a hash key, and a dense rank.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file (must already be initialised).
    score_name : str
        Column name in the parquet file to use as the score value.
    score_path : str
        Path to the source parquet file.
    table_name : str
        Name for the new table inside the database.
    score_type : str
        ``"variant"`` — key columns are ``chrom/pos/ref/alt``.
        ``"gene"`` — key column is ``ensg``.
    """
    score_path = os.path.abspath(score_path)

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Database not found: {db_path}")
    if not os.path.isfile(score_path):
        raise FileNotFoundError(f"Parquet file not found: {score_path}")
    if score_type not in ("variant", "gene"):
        raise ValueError(f"score_type must be 'variant' or 'gene', got {score_type!r}")

    con = duckdb.connect(db_path)
    try:
        tables = {
            r[0]
            for r in con.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema = 'main'"
            ).fetchall()
        }
        if "metadata" not in tables:
            raise RuntimeError(
                "Database has no metadata table. Run initialize_db first."
            )
        if table_name in tables:
            raise ValueError(
                f"Table {table_name!r} already exists in the database."
            )

        parquet_columns = {
            r[0]
            for r in con.execute(
                f"SELECT name FROM parquet_schema('{score_path}')"
            ).fetchall()
        }
        if score_name not in parquet_columns:
            raise ValueError(
                f"Column {score_name!r} not found in parquet. "
                f"Available: {sorted(parquet_columns)}"
            )

        if score_type == "variant":
            missing = [k for k in ("chrom", "pos", "ref", "alt") if k not in parquet_columns]
            if missing:
                raise ValueError(
                    f"Variant key column(s) missing from parquet: {missing}"
                )
            key_select = (
                "chrom, "
                "pos::BIGINT AS pos, "
                "ref, "
                "alt, "
                "NULL::VARCHAR AS ensg"
            )
            hash_expr = (
                "hash(chrom || '|' || CAST(pos AS VARCHAR) "
                "|| '|' || ref || '|' || alt)"
            )
        else:
            if "ensg" not in parquet_columns:
                raise ValueError("Key column 'ensg' not found in parquet.")
            key_select = (
                "NULL::VARCHAR AS chrom, "
                "NULL::BIGINT AS pos, "
                "NULL::VARCHAR AS ref, "
                "NULL::VARCHAR AS alt, "
                "ensg"
            )
            hash_expr = "NULL::UBIGINT"

        quoted_table = f'"{table_name}"'
        quoted_score = f'"{score_name}"'

        create_sql = f"""
        CREATE TABLE {quoted_table} AS
        SELECT
            {key_select},
            {hash_expr} AS "key",
            {quoted_score}::DOUBLE AS score,
            NULL::DOUBLE AS temp_1,
            NULL::DOUBLE AS temp_2,
            CASE WHEN {quoted_score} IS NOT NULL
                 THEN DENSE_RANK() OVER (ORDER BY {quoted_score} NULLS LAST)
                 ELSE NULL
            END::INTEGER AS "rank"
        FROM read_parquet('{score_path}')
        ORDER BY score NULLS LAST;
        """
        con.execute(create_sql)

        row_count = con.execute(
            f"SELECT COUNT(*) FROM {quoted_table}"
        ).fetchone()[0]
        scored_count = con.execute(
            f"SELECT COUNT(score) FROM {quoted_table}"
        ).fetchone()[0]
        unique_keys = con.execute(
            f'SELECT COUNT(DISTINCT "key") FROM {quoted_table}'
        ).fetchone()[0]
        has_dupes = unique_keys < row_count

        con.execute(
            "INSERT INTO metadata (score_name, score_path, table_name, table_type, deduped) "
            "VALUES (?, ?, ?, 'score', ?)",
            [score_name, score_path, table_name, not has_dupes],
        )

        print(
            f"Ingested {score_name!r} → {table_name!r}: "
            f"{row_count} rows ({scored_count} scored, "
            f"{row_count - scored_count} null)",
            file=sys.stderr,
        )
        if has_dupes:
            print(
                f"  WARNING: {row_count - unique_keys:,} duplicate key(s) "
                f"detected. Run remove_duplicates before pairwise operations.",
                file=sys.stderr,
            )
    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ingest a score column from a parquet file into a DuckDB table.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--score_name", required=True,
        help="Name of the score column in the parquet file.",
    )
    parser.add_argument(
        "--score_path", required=True,
        help="Path to the source parquet file.",
    )
    parser.add_argument(
        "--table_name", required=True,
        help="Name for the DuckDB table.",
    )
    parser.add_argument(
        "--score_type", required=True, choices=["variant", "gene"],
        help="Key type: 'variant' (chrom/pos/ref/alt) or 'gene' (ensg).",
    )
    args = parser.parse_args()
    ingest_score(args.db, args.score_name, args.score_path,
                 args.table_name, args.score_type)


if __name__ == "__main__":
    main()
