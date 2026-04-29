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
    has_obs_exp: bool


def _validate_inputs(
    con: duckdb.DuckDBPyConnection,
    table_names: list[str],
) -> list[_EvalInfo]:
    """Check every table exists, is a deduped eval, and return metadata."""
    infos: list[_EvalInfo] = []
    for tbl in table_names:
        row = con.execute(
            "SELECT table_type, deduped, source_column, analysis_level, "
            "eval_column, case_column, observed_column "
            "FROM metadata WHERE table_name = ?",
            [tbl],
        ).fetchone()
        if row is None:
            raise ValueError(f"Table {tbl!r} not found in metadata.")
        ttype, deduped, source_column, analysis_level, eval_col, case_col, obs_col = row
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
            has_obs_exp=obs_col is not None,
        ))
    return infos


def _data_col_refs(alias: str, info: _EvalInfo) -> list[str]:
    """Return SELECT fragments for the eval data columns of one input table.

    Tables with a boolean column contribute ``is_pos`` (renamed to the
    source label).  Tables with count columns contribute ``n_case`` and
    ``n_ctrl`` (prefixed with the source label).  Tables with
    observed/expected columns contribute ``observed`` and ``expected``
    (prefixed with the source label).
    """
    label = info.source_column
    refs: list[str] = []
    if info.has_bool:
        refs.append(f'{alias}.is_pos AS "{label}"')
    if info.has_counts:
        refs.append(f'{alias}.n_case AS "{label}_n_case"')
        refs.append(f'{alias}.n_ctrl AS "{label}_n_ctrl"')
    if info.has_obs_exp:
        refs.append(f'{alias}.observed AS "{label}_observed"')
        refs.append(f'{alias}.expected AS "{label}_expected"')
    return refs


def _renamed_col_names(info: _EvalInfo) -> list[str]:
    """Return the output column names that an eval table contributes."""
    label = info.source_column
    names: list[str] = []
    if info.has_bool:
        names.append(label)
    if info.has_counts:
        names.append(f"{label}_n_case")
        names.append(f"{label}_n_ctrl")
    if info.has_obs_exp:
        names.append(f"{label}_observed")
        names.append(f"{label}_expected")
    return names


def merge_evals(
    db_path: str,
    table_names: list[str],
    output_table: str,
    *,
    memory_limit: str | None = None,
) -> None:
    """Merge multiple eval tables into a single wide output table via union.

    Uses cascading 2-way FULL OUTER JOINs to bound memory usage: each step
    joins the running result with one new table, writes an intermediate, and
    drops the previous one.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    table_names : list[str]
        Names of the input eval tables to merge.
    output_table : str
        Name for the output merged table.
    memory_limit : str or None
        Optional DuckDB memory limit (e.g. ``'8GB'``).
    """
    if len(table_names) < 2:
        raise ValueError("At least two input tables are required.")

    con = duckdb.connect(db_path)
    try:
        if memory_limit:
            con.execute(f"SET memory_limit = '{memory_limit}'")
        con.execute("SET preserve_insertion_order = false")

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
        n = len(infos)

        first_info = infos[0]
        first_data_refs = _data_col_refs("src", first_info)
        first_key_refs = ", ".join(f"src.{k}" for k in KEY_COLS)
        first_select = f"{first_key_refs}, {', '.join(first_data_refs)}"

        cascade_name = f"_cascade_{output_table}_0"
        con.execute(
            f'CREATE TABLE "{cascade_name}" AS\n'
            f'SELECT {first_select}\n'
            f'FROM "{first_info.table_name}" src'
        )
        accumulated_cols = _renamed_col_names(first_info)
        print(
            f"  merge step 1/{n}: materialised {first_info.source_column!r}",
            file=sys.stderr,
        )

        for i in range(1, n):
            prev_name = cascade_name
            cascade_name = f"_cascade_{output_table}_{i}"
            info = infos[i]

            key_exprs = [
                f"COALESCE(l.{k}, r.{k}) AS {k}" for k in KEY_COLS
            ]
            left_refs = [f'l."{c}"' for c in accumulated_cols]
            right_refs = _data_col_refs("r", info)
            select = ", ".join(key_exprs + left_refs + right_refs)

            con.execute(
                f'CREATE TABLE "{cascade_name}" AS\n'
                f"SELECT {select}\n"
                f'FROM "{prev_name}" l\n'
                f'FULL OUTER JOIN "{info.table_name}" r '
                f'ON l."key" = r."key"'
            )
            con.execute(f'DROP TABLE "{prev_name}"')
            accumulated_cols.extend(_renamed_col_names(info))

            row_count_so_far = con.execute(
                f'SELECT COUNT(*) FROM "{cascade_name}"'
            ).fetchone()[0]
            print(
                f"  merge step {i + 1}/{n}: joined {info.source_column!r} "
                f"({row_count_so_far:,} rows)",
                file=sys.stderr,
            )

        quoted_output = f'"{output_table}"'
        con.execute(
            f'ALTER TABLE "{cascade_name}" RENAME TO {quoted_output}'
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
    parser.add_argument(
        "--memory_limit", default=None,
        help="DuckDB memory limit (e.g. '8GB'). Default: DuckDB auto.",
    )
    args = parser.parse_args()
    merge_evals(args.db, args.tables, args.output_table,
                memory_limit=args.memory_limit)


if __name__ == "__main__":
    main()
