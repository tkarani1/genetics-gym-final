#!/usr/bin/env python3
"""Ingest a boolean eval column from a parquet file into a DuckDB table."""
from __future__ import annotations

import argparse
import os
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event


def ingest_eval(
    db_path: str,
    eval_name: str,
    eval_path: str,
    table_name: str,
    analysis_level: str,
) -> None:
    """Read evaluation data from a parquet/TSV file and materialise it as a
    DuckDB table with canonical key columns, a hash key, and evaluation
    columns.

    Variant-level tables store a boolean ``is_pos`` column (read from the
    column named by *eval_name*) with ``n_case``/``n_ctrl`` set to NULL.

    Gene-level tables store integer ``n_case`` and ``n_ctrl`` columns (read
    from identically-named source columns) with ``is_pos`` set to NULL.
    For gene-level ingestion *eval_name* is a descriptive label only and
    does not need to match a source column.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file (must already be initialised).
    eval_name : str
        For variant-level: column name in the source file to use as ``is_pos``.
        For gene-level: descriptive label stored in metadata.
    eval_path : str
        Path to the source parquet (or TSV/CSV) file.
    table_name : str
        Name for the new table inside the database.
    analysis_level : str
        ``"variant"`` — key columns are ``chrom/pos/ref/alt``.
        ``"gene"`` — key column is ``ensg``.
    """
    eval_path = os.path.abspath(eval_path)

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Database not found: {db_path}")
    if not os.path.isfile(eval_path):
        raise FileNotFoundError(f"Source file not found: {eval_path}")
    if analysis_level not in ("variant", "gene"):
        raise ValueError(f"analysis_level must be 'variant' or 'gene', got {analysis_level!r}")

    lower = eval_path.lower()
    is_parquet = lower.endswith(".parquet")
    read_fn = (f"read_parquet('{eval_path}')" if is_parquet
               else f"read_csv_auto('{eval_path}')")

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

        if is_parquet:
            source_columns = {
                r[0] for r in con.execute(
                    f"SELECT name FROM parquet_schema('{eval_path}')"
                ).fetchall()
            }
        else:
            desc = con.execute(f"DESCRIBE SELECT * FROM {read_fn}").fetchall()
            source_columns = {r[0] for r in desc}

        if analysis_level == "variant":
            if eval_name not in source_columns:
                raise ValueError(
                    f"Column {eval_name!r} not found in source. "
                    f"Available: {sorted(source_columns)}"
                )
            missing = [k for k in ("chrom", "pos", "ref", "alt") if k not in source_columns]
            if missing:
                raise ValueError(
                    f"Variant key column(s) missing from source: {missing}"
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
            quoted_eval = f'"{eval_name}"'
            data_select = (
                f"{quoted_eval}::BOOLEAN AS is_pos, "
                "NULL::INTEGER AS n_case, "
                "NULL::INTEGER AS n_ctrl"
            )
        else:
            if "ensg" not in source_columns:
                raise ValueError("Key column 'ensg' not found in source.")

            has_bool = eval_name in source_columns
            has_counts = "n_case" in source_columns and "n_ctrl" in source_columns
            if not has_bool and not has_counts:
                raise ValueError(
                    f"Gene-level source must contain either "
                    f"{eval_name!r} (boolean) or n_case/n_ctrl (counts). "
                    f"Available: {sorted(source_columns)}"
                )

            key_select = (
                "NULL::VARCHAR AS chrom, "
                "NULL::BIGINT AS pos, "
                "NULL::VARCHAR AS ref, "
                "NULL::VARCHAR AS alt, "
                "ensg"
            )
            hash_expr = "hash(ensg)"

            quoted_eval = f'"{eval_name}"'
            is_pos_expr = f"{quoted_eval}::BOOLEAN AS is_pos" if has_bool else "NULL::BOOLEAN AS is_pos"
            n_case_expr = "n_case::INTEGER AS n_case" if has_counts else "NULL::INTEGER AS n_case"
            n_ctrl_expr = "n_ctrl::INTEGER AS n_ctrl" if has_counts else "NULL::INTEGER AS n_ctrl"
            data_select = f"{is_pos_expr}, {n_case_expr}, {n_ctrl_expr}"

        quoted_table = f'"{table_name}"'

        create_sql = f"""
        CREATE TABLE {quoted_table} AS
        SELECT
            {key_select},
            {hash_expr} AS "key",
            {data_select}
        FROM {read_fn};
        """
        con.execute(create_sql)

        row_count = con.execute(
            f"SELECT COUNT(*) FROM {quoted_table}"
        ).fetchone()[0]
        unique_keys = con.execute(
            f'SELECT COUNT(DISTINCT "key") FROM {quoted_table}'
        ).fetchone()[0]
        has_dupes = unique_keys < row_count

        con.execute(
            "INSERT INTO metadata (source_column, source_path, table_name, table_type, analysis_level, deduped) "
            "VALUES (?, ?, ?, 'eval', ?, ?)",
            [eval_name, eval_path, table_name, analysis_level, not has_dupes],
        )
        log_event(
            con, "ingest_eval", "create_table", table_name,
            f"source_column={eval_name}, source_path={eval_path}, "
            f"analysis_level={analysis_level}, rows={row_count}",
        )

        if analysis_level == "variant":
            pos_count = con.execute(
                f"SELECT COUNT(*) FROM {quoted_table} WHERE is_pos = TRUE"
            ).fetchone()[0]
            neg_count = con.execute(
                f"SELECT COUNT(*) FROM {quoted_table} WHERE is_pos = FALSE"
            ).fetchone()[0]
            null_count = row_count - pos_count - neg_count
            print(
                f"Ingested {eval_name!r} → {table_name!r}: "
                f"{row_count} rows ({pos_count} positive, "
                f"{neg_count} negative, {null_count} null)",
                file=sys.stderr,
            )
        else:
            parts = [f"{row_count} rows"]
            if has_bool:
                pos_count = con.execute(
                    f"SELECT COUNT(*) FROM {quoted_table} WHERE is_pos = TRUE"
                ).fetchone()[0]
                neg_count = con.execute(
                    f"SELECT COUNT(*) FROM {quoted_table} WHERE is_pos = FALSE"
                ).fetchone()[0]
                parts.append(f"{pos_count} positive, {neg_count} negative")
            if has_counts:
                total_case = con.execute(
                    f"SELECT SUM(n_case) FROM {quoted_table}"
                ).fetchone()[0] or 0
                total_ctrl = con.execute(
                    f"SELECT SUM(n_ctrl) FROM {quoted_table}"
                ).fetchone()[0] or 0
                parts.append(f"Σn_case={total_case:,}, Σn_ctrl={total_ctrl:,}")
            print(
                f"Ingested {eval_name!r} → {table_name!r}: "
                f"{' ('.join(parts)})",
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
        description="Ingest a boolean eval column from a parquet file into a DuckDB table.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--eval_name", required=True,
        help="Variant: column name to use as is_pos. Gene: descriptive label.",
    )
    parser.add_argument(
        "--eval_path", required=True,
        help="Path to the source file (parquet, TSV, or CSV).",
    )
    parser.add_argument(
        "--table_name", required=True,
        help="Name for the DuckDB table.",
    )
    parser.add_argument(
        "--analysis_level", required=True, choices=["variant", "gene"],
        help="Key type: 'variant' (chrom/pos/ref/alt) or 'gene' (ensg).",
    )
    args = parser.parse_args()
    ingest_eval(args.db, args.eval_name, args.eval_path,
                args.table_name, args.analysis_level)


if __name__ == "__main__":
    main()
