#!/usr/bin/env python3
"""Merge multiple eval tables into a single wide table via FULL OUTER JOIN."""
from __future__ import annotations

import argparse
import sys
from typing import NamedTuple

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

KEY_COLS = ("chrom", "pos", "ref", "alt", "ensg", '"key"')


class _EvalInfo(NamedTuple):
    table_name: str
    source_column: str
    analysis_level: str
    has_bool: bool
    has_counts: bool


def _validate_inputs(
    con: duckdb.DuckDBPyConnection,
    table_names: list[str],
) -> list[_EvalInfo]:
    """Check every table exists, is a deduped eval, and return metadata."""
    infos: list[_EvalInfo] = []
    for tbl in table_names:
        row = con.execute(
            "SELECT table_type, deduped, source_column, analysis_level, "
            "eval_column, case_column "
            "FROM metadata WHERE table_name = ?",
            [tbl],
        ).fetchone()
        if row is None:
            raise ValueError(f"Table {tbl!r} not found in metadata.")
        ttype, deduped, source_column, analysis_level, eval_col, case_col = row
        if ttype != "eval":
            raise ValueError(f"Table {tbl!r} is type {ttype!r}, not 'eval'.")
        if not deduped:
            raise ValueError(
                f"Table {tbl!r} has not been deduped. "
                f"Run remove_duplicates first."
            )
        infos.append(_EvalInfo(
            tbl, source_column, analysis_level,
            has_bool=eval_col is not None,
            has_counts=case_col is not None,
        ))
    return infos


def _data_col_refs(alias: str, info: _EvalInfo) -> list[str]:
    """Return SELECT fragments for the eval data columns of one input table.

    Tables with a boolean column contribute ``is_pos`` (renamed to the
    source label).  Tables with count columns contribute ``n_case`` and
    ``n_ctrl`` (prefixed with the source label).  Tables with both
    contribute all three.
    """
    label = info.source_column
    refs: list[str] = []
    if info.has_bool:
        refs.append(f'{alias}.is_pos AS "{label}"')
    if info.has_counts:
        refs.append(f'{alias}.n_case AS "{label}_n_case"')
        refs.append(f'{alias}.n_ctrl AS "{label}_n_ctrl"')
    return refs


def merge_evals(
    db_path: str,
    table_names: list[str],
    output_table: str,
) -> None:
    """Merge multiple eval tables into a single wide output table via union.

    Variant-level input tables contribute their ``is_pos`` column (renamed
    to the table's ``source_column``).  Gene-level input tables contribute
    ``n_case`` and ``n_ctrl`` columns (prefixed with the source label).
    The union is a ``FULL OUTER JOIN`` on ``key``, so the result contains
    all rows from any input table; values are NULL where a table lacks a key.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    table_names : list[str]
        Names of the input eval tables to merge.
    output_table : str
        Name for the output merged table.
    """
    if len(table_names) < 2:
        raise ValueError("At least two input tables are required.")

    con = duckdb.connect(db_path)
    try:
        existing = {
            r[0]
            for r in con.execute(
                "SELECT table_name FROM information_schema.tables "
                "WHERE table_schema = 'main'"
            ).fetchall()
        }
        if output_table in existing:
            raise ValueError(
                f"Output table {output_table!r} already exists."
            )

        infos = _validate_inputs(con, table_names)

        aliases = [f"t{i}" for i in range(len(infos))]
        first = aliases[0]

        key_coalesce = ", ".join(
            f'COALESCE({", ".join(f"{a}.{k}" for a in aliases)}) AS {k}'
            for k in KEY_COLS
        )

        eval_col_parts: list[str] = []
        for i, info in enumerate(infos):
            eval_col_parts.extend(_data_col_refs(aliases[i], info))
        eval_col_refs = ", ".join(eval_col_parts)

        select_clause = f"{key_coalesce}, {eval_col_refs}"

        subselects = [f'"{info.table_name}"' for info in infos]
        from_clause = f"{subselects[0]} AS {aliases[0]}"
        for i in range(1, len(subselects)):
            from_clause += (
                f'\n    FULL OUTER JOIN {subselects[i]} AS {aliases[i]} '
                f'ON {first}."key" = {aliases[i]}."key"'
            )

        quoted_output = f'"{output_table}"'
        con.execute(
            f"CREATE TABLE {quoted_output} AS\n"
            f"SELECT {select_clause}\nFROM {from_clause}"
        )

        row_count = con.execute(
            f"SELECT COUNT(*) FROM {quoted_output}"
        ).fetchone()[0]
        col_count = len(
            con.execute(
                f"SELECT * FROM {quoted_output} LIMIT 0"
            ).description
        )

        source_names = ", ".join(i.source_column for i in infos)
        analysis_levels = {i.analysis_level for i in infos}
        merged_analysis_level = analysis_levels.pop() if len(analysis_levels) == 1 else "variant"
        con.execute(
            "INSERT INTO metadata "
            "(source_column, source_path, table_name, table_type, analysis_level, deduped) "
            "VALUES (?, 'union', ?, 'merged_evals', ?, TRUE)",
            [source_names, output_table, merged_analysis_level],
        )

        input_tables = ", ".join(i.table_name for i in infos)
        log_event(
            con, "merge_evals", "create_table", output_table,
            f"input_tables=[{input_tables}], rows={row_count}",
        )

        print(
            f"Merged {len(infos)} eval tables (union) → {output_table!r}: "
            f"{row_count:,} rows, {col_count} columns",
            file=sys.stderr,
        )
    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge multiple eval tables into a single wide table via union.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--tables", required=True, nargs="+",
        help="Names of the input eval tables to merge.",
    )
    parser.add_argument(
        "--output_table", required=True,
        help="Name for the output merged eval table.",
    )
    args = parser.parse_args()
    merge_evals(args.db, args.tables, args.output_table)


if __name__ == "__main__":
    main()
