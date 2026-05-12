#!/usr/bin/env python3
"""Export an analysis table to Parquet with optional filter joins or row selection.

Two modes of operation:

1. **Join mode** (``--filter_table``): LEFT JOIN all boolean columns from a
   wide filter table onto the analysis table and write the result to Parquet.
   The database table is NOT modified.

2. **Select mode** (``--select_filter``): Use one or more filter names as a
   WHERE clause to subset the analysis table.  Only rows where the specified
   filter(s) are TRUE are exported.  Supports AND/OR logic via ``--select_logic``.

Both modes can be combined: first subset rows via select_filter, then append
all filter columns from filter_table.
"""
from __future__ import annotations

import argparse
import os
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

ELIGIBLE_TYPES = ("merged_scores", "merged_evals", "merged_analysis")
KEY_NAMES = {"chrom", "pos", "ref", "alt", "ensg", "key"}


def _resolve_filter_location(
    con: duckdb.DuckDBPyConnection,
    filter_name: str,
) -> tuple[str, str]:
    """Find which filter table contains a given filter column.

    Returns (table_name, join_column) where join_column is 'key' for
    variant-level or 'ensg' for gene-level filters.
    """
    filter_tables = con.execute(
        "SELECT DISTINCT table_name, analysis_level FROM metadata "
        "WHERE table_type = 'filter'"
    ).fetchall()

    for tbl, level in filter_tables:
        cols = {
            d[0] for d in con.execute(
                f'SELECT * FROM "{tbl}" LIMIT 0'
            ).description
        }
        if filter_name in cols:
            join_col = '"key"' if level == "variant" else "ensg"
            return tbl, join_col

    available = []
    for tbl, _ in filter_tables:
        cols = [
            d[0] for d in con.execute(
                f'SELECT * FROM "{tbl}" LIMIT 0'
            ).description
            if d[0] not in KEY_NAMES
        ]
        available.extend(cols)
    raise ValueError(
        f"Filter {filter_name!r} not found in any filter table. "
        f"Available: {sorted(set(available))}"
    )


def join_filters(
    db_path: str,
    table_name: str,
    output_path: str,
    *,
    filter_table: str | None = None,
    select_filters: list[str] | None = None,
    select_logic: str = "and",
) -> None:
    """Export an analysis table to Parquet with optional filter operations.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    table_name : str
        Name of the source table (merged_scores, merged_evals, or
        merged_analysis).
    output_path : str
        Destination Parquet file path.
    filter_table : str or None
        If provided, LEFT JOIN all filter columns from this table onto the
        output.
    select_filters : list[str] or None
        If provided, only export rows where these filter(s) are TRUE.
    select_logic : str
        ``"and"`` or ``"or"`` — how to combine multiple select_filters.
    """
    output_path = os.path.abspath(output_path)

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Database not found: {db_path}")
    if os.path.exists(output_path):
        raise FileExistsError(f"Output file already exists: {output_path}")
    if filter_table is None and select_filters is None:
        raise ValueError(
            "At least one of --filter_table or --select_filter must be provided."
        )
    if select_logic not in ("and", "or"):
        raise ValueError(
            f"--select_logic must be 'and' or 'or', got {select_logic!r}"
        )

    con = duckdb.connect(db_path, read_only=True)
    try:
        row = con.execute(
            "SELECT table_type, analysis_level FROM metadata "
            "WHERE table_name = ? LIMIT 1",
            [table_name],
        ).fetchone()
        if row is None:
            raise ValueError(f"Table {table_name!r} not found in metadata.")

        table_type, target_level = row
        if table_type not in ELIGIBLE_TYPES:
            raise ValueError(
                f"Table {table_name!r} is type {table_type!r}. "
                f"join_filters only works on merged tables "
                f"({', '.join(ELIGIBLE_TYPES)})."
            )

        target_cols = [
            desc[0]
            for desc in con.execute(
                f'SELECT * FROM "{table_name}" LIMIT 0'
            ).description
        ]

        # Build select clause for the base table
        target_refs = ", ".join(f't."{c}"' for c in target_cols)

        # Handle --select_filter: build WHERE clause and required JOINs
        select_joins: list[str] = []
        where_conditions: list[str] = []
        select_join_aliases: dict[str, str] = {}

        if select_filters:
            for i, sf in enumerate(select_filters):
                ftbl, join_col = _resolve_filter_location(con, sf)
                alias = f"sf{i}"
                if ftbl not in select_join_aliases:
                    select_join_aliases[ftbl] = alias
                    if join_col == "ensg":
                        select_joins.append(
                            f'LEFT JOIN "{ftbl}" {alias} '
                            f"ON t.ensg = {alias}.ensg"
                        )
                    else:
                        select_joins.append(
                            f'LEFT JOIN "{ftbl}" {alias} '
                            f'ON t."key" = {alias}."key"'
                        )
                else:
                    alias = select_join_aliases[ftbl]

                where_conditions.append(f'{alias}."{sf}" = TRUE')

        where_clause = ""
        if where_conditions:
            joiner = " AND " if select_logic == "and" else " OR "
            where_clause = f"WHERE {joiner.join(where_conditions)}"

        # Handle --filter_table: build additional JOIN for column append
        filter_join = ""
        filter_refs = ""

        if filter_table:
            tables = {
                r[0]
                for r in con.execute(
                    "SELECT table_name FROM information_schema.tables "
                    "WHERE table_schema = 'main'"
                ).fetchall()
            }
            if filter_table not in tables:
                raise ValueError(
                    f"Filter table {filter_table!r} does not exist."
                )

            filter_level_row = con.execute(
                "SELECT DISTINCT analysis_level FROM metadata "
                "WHERE table_name = ?",
                [filter_table],
            ).fetchone()
            filter_level = (
                filter_level_row[0] if filter_level_row else target_level
            )

            filter_cols = [
                desc[0]
                for desc in con.execute(
                    f'SELECT * FROM "{filter_table}" LIMIT 0'
                ).description
                if desc[0] not in KEY_NAMES
            ]

            if not filter_cols:
                raise ValueError(
                    f"Filter table {filter_table!r} has no non-key columns."
                )

            overlap = set(filter_cols) & set(target_cols)
            if overlap:
                raise ValueError(
                    f"Target table already contains filter columns: "
                    f"{sorted(overlap)}"
                )

            # Avoid alias collision with select_filter joins
            ft_alias = "ft"
            if filter_level == "gene" and target_level == "variant":
                filter_join = (
                    f'LEFT JOIN "{filter_table}" {ft_alias} '
                    f"ON t.ensg = {ft_alias}.ensg"
                )
            else:
                filter_join = (
                    f'LEFT JOIN "{filter_table}" {ft_alias} '
                    f'ON t."key" = {ft_alias}."key"'
                )

            filter_refs = ", " + ", ".join(
                f'{ft_alias}."{c}"' for c in filter_cols
            )

        # Build and execute the full query
        all_joins = " ".join(select_joins)
        if filter_join:
            all_joins = f"{all_joins} {filter_join}".strip()

        query = (
            f"SELECT {target_refs}{filter_refs} "
            f'FROM "{table_name}" t '
            f"{all_joins} "
            f"{where_clause}"
        )

        con.execute(
            f"COPY ({query}) TO '{output_path}' (FORMAT PARQUET)"
        )

        row_count = con.execute(
            f"SELECT COUNT(*) FROM read_parquet('{output_path}')"
        ).fetchone()[0]

    finally:
        con.close()

    # Log the event (read-write connection)
    con_rw = duckdb.connect(db_path)
    try:
        details_parts = [f"source_table={table_name}"]
        if filter_table:
            details_parts.append(f"filter_table={filter_table}")
        if select_filters:
            details_parts.append(
                f"select_filter={select_filters}, logic={select_logic}"
            )
        details_parts.append(f"rows={row_count}, output={output_path}")

        log_event(
            con_rw, "join_filters", "export", table_name,
            ", ".join(details_parts),
        )
    finally:
        con_rw.close()

    print(
        f"Exported {table_name!r} -> {output_path}: {row_count:,} rows",
        file=sys.stderr,
    )
    if select_filters:
        print(
            f"  select_filter: {select_filters} (logic={select_logic})",
            file=sys.stderr,
        )
    if filter_table:
        filter_col_count = len(filter_refs.split(",")) - 1 if filter_refs else 0
        print(
            f"  filter_table: {filter_table} "
            f"({filter_col_count} columns appended)",
            file=sys.stderr,
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Export an analysis table to Parquet with optional filter "
            "joins or row selection."
        ),
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--table_name", required=True,
        help="Name of the source merged table.",
    )
    parser.add_argument(
        "--output_path", required=True,
        help="Destination Parquet file path.",
    )
    parser.add_argument(
        "--filter_table", default=None,
        help="Wide filter table to LEFT JOIN onto the output "
             "(e.g. 'variant_filters').",
    )
    parser.add_argument(
        "--select_filter", nargs="+", default=None,
        help="Filter name(s) to use as row selection (only TRUE rows exported).",
    )
    parser.add_argument(
        "--select_logic", default="and", choices=["and", "or"],
        help="Logic for combining multiple --select_filter values "
             "(default: and).",
    )
    args = parser.parse_args()
    join_filters(
        args.db, args.table_name, args.output_path,
        filter_table=args.filter_table,
        select_filters=args.select_filter,
        select_logic=args.select_logic,
    )


if __name__ == "__main__":
    main()
