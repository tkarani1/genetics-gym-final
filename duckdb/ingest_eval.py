#!/usr/bin/env python3
"""Ingest evaluation data from a parquet/TSV file into a DuckDB table."""
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
    *,
    eval_column: str | None = None,
    case_column: str | None = None,
    ctrl_column: str | None = None,
    observed_column: str | None = None,
    expected_column: str | None = None,
) -> None:
    """Read evaluation data from a source file and materialise it as a DuckDB
    table with canonical key columns, a hash key, and evaluation columns.

    Every eval table has five data columns — ``is_pos`` (boolean),
    ``n_case`` (integer), ``n_ctrl`` (integer), ``observed`` (integer),
    and ``expected`` (double).  Which are populated vs NULL depends on the
    arguments provided:

    * If *eval_column* is given, the named source column is read as
      ``is_pos``; otherwise ``is_pos`` is NULL.
    * If *case_column* and *ctrl_column* are given, those source columns
      are read as ``n_case`` and ``n_ctrl``; otherwise both are NULL.
    * If *observed_column* and *expected_column* are given, those source
      columns are read as ``observed`` and ``expected``; otherwise both
      are NULL.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file (must already be initialised).
    eval_name : str
        Human-readable label for this evaluation, stored in metadata.
    eval_path : str
        Path to the source parquet (or TSV/CSV) file.
    table_name : str
        Name for the new table inside the database.
    analysis_level : str
        ``"variant"`` — key columns are ``chrom/pos/ref/alt``.
        ``"gene"`` — key column is ``ensg``.
    eval_column : str or None
        Source column to read as ``is_pos``.
    case_column : str or None
        Source column to read as ``n_case``.  Must be paired with *ctrl_column*.
    ctrl_column : str or None
        Source column to read as ``n_ctrl``.  Must be paired with *case_column*.
    observed_column : str or None
        Source column to read as ``observed``.  Must be paired with
        *expected_column*.
    expected_column : str or None
        Source column to read as ``expected``.  Must be paired with
        *observed_column*.
    """
    eval_path = os.path.abspath(eval_path)

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Database not found: {db_path}")
    if not os.path.isfile(eval_path) and not os.path.isdir(eval_path):
        raise FileNotFoundError(f"Source file not found: {eval_path}")
    if analysis_level not in ("variant", "gene"):
        raise ValueError(f"analysis_level must be 'variant' or 'gene', got {analysis_level!r}")

    has_bool = eval_column is not None
    has_counts = case_column is not None or ctrl_column is not None
    has_obs_exp = observed_column is not None or expected_column is not None
    if not has_bool and not has_counts and not has_obs_exp:
        raise ValueError(
            "At least one of --eval_column, --case_column/--ctrl_column, "
            "or --observed_column/--expected_column must be provided."
        )
    if (case_column is None) != (ctrl_column is None):
        raise ValueError(
            "--case_column and --ctrl_column must be provided together."
        )
    if (observed_column is None) != (expected_column is None):
        raise ValueError(
            "--observed_column and --expected_column must be provided together."
        )

    lower = eval_path.lower()
    is_parquet = lower.endswith(".parquet")
    if is_parquet and os.path.isdir(eval_path):
        read_fn = f"read_parquet('{eval_path}/*.parquet')"
    elif is_parquet:
        read_fn = f"read_parquet('{eval_path}')"
    else:
        read_fn = f"read_csv_auto('{eval_path}')"

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
            if os.path.isdir(eval_path):
                schema_target = next(
                    os.path.join(eval_path, f)
                    for f in os.listdir(eval_path)
                    if f.endswith(".parquet")
                )
            else:
                schema_target = eval_path
            source_columns = {
                r[0] for r in con.execute(
                    f"SELECT name FROM parquet_schema('{schema_target}')"
                ).fetchall()
            }
        else:
            desc = con.execute(f"DESCRIBE SELECT * FROM {read_fn}").fetchall()
            source_columns = {r[0] for r in desc}

        if has_bool and eval_column not in source_columns:
            raise ValueError(
                f"eval_column {eval_column!r} not found in source. "
                f"Available: {sorted(source_columns)}"
            )
        if has_counts:
            for col_arg, col_name in [("case_column", case_column),
                                       ("ctrl_column", ctrl_column)]:
                if col_name not in source_columns:
                    raise ValueError(
                        f"{col_arg} {col_name!r} not found in source. "
                        f"Available: {sorted(source_columns)}"
                    )
        if has_obs_exp:
            for col_arg, col_name in [("observed_column", observed_column),
                                       ("expected_column", expected_column)]:
                if col_name not in source_columns:
                    raise ValueError(
                        f"{col_arg} {col_name!r} not found in source. "
                        f"Available: {sorted(source_columns)}"
                    )

        def _col_or_null(col: str, cast: str | None = None) -> str:
            """Include a source column if present, otherwise NULL."""
            if col in source_columns:
                expr = f'"{col}"'
                if cast:
                    expr = f'{expr}::{cast}'
                return f'{expr} AS {col}'
            return f'NULL::{cast or "VARCHAR"} AS {col}'

        if analysis_level == "variant":
            missing = [k for k in ("chrom", "pos", "ref", "alt") if k not in source_columns]
            if missing:
                raise ValueError(
                    f"Variant key column(s) missing from source: {missing}"
                )
            hash_expr = (
                "hash(chrom || '|' || CAST(pos AS VARCHAR) "
                "|| '|' || ref || '|' || alt)"
            )
        else:
            if "ensg" not in source_columns:
                raise ValueError("Key column 'ensg' not found in source.")
            hash_expr = "hash(ensg)"

        key_select = ", ".join([
            _col_or_null("chrom"),
            _col_or_null("pos", "BIGINT"),
            _col_or_null("ref"),
            _col_or_null("alt"),
            _col_or_null("ensg"),
        ])

        is_pos_expr = (f'"{eval_column}"::BOOLEAN AS is_pos'
                       if has_bool else "NULL::BOOLEAN AS is_pos")
        n_case_expr = (f'"{case_column}"::INTEGER AS n_case'
                       if has_counts else "NULL::INTEGER AS n_case")
        n_ctrl_expr = (f'"{ctrl_column}"::INTEGER AS n_ctrl'
                       if has_counts else "NULL::INTEGER AS n_ctrl")
        observed_expr = (f'"{observed_column}"::INTEGER AS observed'
                         if has_obs_exp else "NULL::INTEGER AS observed")
        expected_expr = (f'"{expected_column}"::DOUBLE AS expected'
                         if has_obs_exp else "NULL::DOUBLE AS expected")
        data_select = (f"{is_pos_expr}, {n_case_expr}, {n_ctrl_expr}, "
                       f"{observed_expr}, {expected_expr}")

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
            "INSERT INTO metadata "
            "(source_column, source_path, table_name, table_type, "
            "analysis_level, deduped, eval_column, case_column, ctrl_column, "
            "observed_column, expected_column) "
            "VALUES (?, ?, ?, 'eval', ?, ?, ?, ?, ?, ?, ?)",
            [eval_name, eval_path, table_name, analysis_level,
             not has_dupes, eval_column, case_column, ctrl_column,
             observed_column, expected_column],
        )
        log_event(
            con, "ingest_eval", "create_table", table_name,
            f"eval_name={eval_name}, eval_column={eval_column}, "
            f"case_column={case_column}, ctrl_column={ctrl_column}, "
            f"observed_column={observed_column}, "
            f"expected_column={expected_column}, "
            f"source_path={eval_path}, "
            f"analysis_level={analysis_level}, rows={row_count}",
        )

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
        if has_obs_exp:
            sum_obs = con.execute(
                f"SELECT SUM(observed) FROM {quoted_table}"
            ).fetchone()[0] or 0
            sum_exp = con.execute(
                f"SELECT SUM(expected) FROM {quoted_table}"
            ).fetchone()[0] or 0.0
            parts.append(f"Σobserved={sum_obs:,}, Σexpected={sum_exp:,.4f}")
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
        description="Ingest evaluation data from a source file into a DuckDB table.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--eval_name", required=True,
        help="Human-readable label for this evaluation.",
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
    parser.add_argument(
        "--eval_column", default=None,
        help="Source column to read as is_pos (boolean).",
    )
    parser.add_argument(
        "--case_column", default=None,
        help="Source column to read as n_case (must pair with --ctrl_column).",
    )
    parser.add_argument(
        "--ctrl_column", default=None,
        help="Source column to read as n_ctrl (must pair with --case_column).",
    )
    parser.add_argument(
        "--observed_column", default=None,
        help="Source column to read as observed (must pair with --expected_column).",
    )
    parser.add_argument(
        "--expected_column", default=None,
        help="Source column to read as expected (must pair with --observed_column).",
    )
    args = parser.parse_args()
    ingest_eval(
        args.db, args.eval_name, args.eval_path,
        args.table_name, args.analysis_level,
        eval_column=args.eval_column,
        case_column=args.case_column,
        ctrl_column=args.ctrl_column,
        observed_column=args.observed_column,
        expected_column=args.expected_column,
    )


if __name__ == "__main__":
    main()
