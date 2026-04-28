#!/usr/bin/env python3
"""Compute pairwise percentiles between two score columns in a wide table.

Output is written to a Parquet file — no scratch columns are modified in
the database.
"""
from __future__ import annotations

import argparse
import os
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

KEY_COLS = ("chrom", "pos", "ref", "alt", "ensg", '"key"')


def pairwise_percentile(
    db_path: str,
    wide_table: str,
    anchor_column: str,
    target_column: str,
    output_path: str,
    *,
    memory_limit: str | None = None,
) -> None:
    """Compute intersection-based pairwise percentiles for two score columns
    and write the result to a Parquet file.

    For the rows where both *anchor_column* and *target_column* are non-null,
    ``CUME_DIST`` is computed for each.  Rows outside the intersection still
    appear in the output with ``NULL`` percentile values.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    wide_table : str
        Name of the wide score table (e.g. ``variant_scores``).
    anchor_column : str
        Score column to use as the anchor.
    target_column : str
        Score column to compute percentiles against the anchor.
    output_path : str
        Destination Parquet file.
    memory_limit : str or None
        Optional DuckDB memory limit (e.g. ``'8GB'``).
    """
    output_path = os.path.abspath(output_path)
    if os.path.exists(output_path):
        raise FileExistsError(f"Output file already exists: {output_path}")

    con = duckdb.connect(db_path, read_only=True)
    try:
        if memory_limit:
            con.execute(f"SET memory_limit = '{memory_limit}'")
        con.execute("SET preserve_insertion_order = false")

        exists = con.execute(
            "SELECT COUNT(*) FROM information_schema.tables "
            "WHERE table_name = ? AND table_schema = 'main'",
            [wide_table],
        ).fetchone()[0]
        if not exists:
            raise ValueError(f"Table {wide_table!r} does not exist.")

        actual_cols = {
            desc[0]
            for desc in con.execute(
                f'SELECT * FROM "{wide_table}" LIMIT 0'
            ).description
        }
        for col_name in (anchor_column, target_column):
            if col_name not in actual_cols:
                raise ValueError(
                    f"Column {col_name!r} not found in {wide_table!r}."
                )

        quoted = f'"{wide_table}"'
        key_select = ", ".join(KEY_COLS)
        both_not_null = (
            f'("{anchor_column}" IS NOT NULL AND "{target_column}" IS NOT NULL)'
        )
        pw_anchor = f"{anchor_column}_pairwise_{target_column}"
        pw_target = f"{target_column}_pairwise_{anchor_column}"

        query = f"""
            SELECT {key_select},
                "{anchor_column}",
                "{target_column}",
                CASE WHEN {both_not_null} THEN CUME_DIST() OVER (
                    PARTITION BY {both_not_null}
                    ORDER BY "{anchor_column}"
                ) END AS "{pw_anchor}",
                CASE WHEN {both_not_null} THEN CUME_DIST() OVER (
                    PARTITION BY {both_not_null}
                    ORDER BY "{target_column}"
                ) END AS "{pw_target}"
            FROM {quoted}
        """

        con.execute(
            f"COPY ({query}) TO '{output_path}' (FORMAT PARQUET)"
        )

        row_count = con.execute(
            f"SELECT COUNT(*) FROM read_parquet('{output_path}')"
        ).fetchone()[0]

        print(
            f"Pairwise percentiles {anchor_column!r} x {target_column!r} "
            f"→ {output_path}: {row_count:,} rows",
            file=sys.stderr,
        )
    finally:
        con.close()

    con_rw = duckdb.connect(db_path)
    try:
        log_event(
            con_rw, "pairwise_percentile", "export_parquet", wide_table,
            f"anchor={anchor_column}, target={target_column}, "
            f"output={output_path}, rows={row_count}",
        )
    finally:
        con_rw.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute pairwise percentiles between two score columns "
            "in a wide table and export to Parquet."
        ),
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--wide_table", required=True,
        help="Name of the wide score table (e.g. variant_scores).",
    )
    parser.add_argument(
        "--anchor_column", required=True,
        help="Score column to use as the anchor.",
    )
    parser.add_argument(
        "--target_column", required=True,
        help="Score column to compute percentiles against the anchor.",
    )
    parser.add_argument(
        "--output_path", required=True,
        help="Destination Parquet file.",
    )
    parser.add_argument(
        "--memory_limit", default=None,
        help="DuckDB memory limit (e.g. '8GB'). Default: DuckDB auto.",
    )
    args = parser.parse_args()
    pairwise_percentile(
        args.db, args.wide_table, args.anchor_column, args.target_column,
        args.output_path, memory_limit=args.memory_limit,
    )


if __name__ == "__main__":
    main()
