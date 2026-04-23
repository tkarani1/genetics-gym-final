#!/usr/bin/env python3
"""Merge multiple score tables into a single wide table."""
from __future__ import annotations

import argparse
import re
import sys
from typing import NamedTuple

import duckdb

from initialize_db import log_event


DEFAULT_DB_PATH = "scores.duckdb"

SET_OPERATIONS = ("intersection", "union", "pairwise")
PERCENTILE_MODES = ("pre", "post", "none")

KEY_COLS = ("chrom", "pos", "ref", "alt", "ensg", '"key"')


class _ScoreInfo(NamedTuple):
    table_name: str
    score_name: str
    analysis_level: str


def _percentile_expr(col: str, alias: str) -> str:
    """SQL expression for null-safe CUME_DIST over *col*."""
    return (
        f'CASE WHEN "{col}" IS NOT NULL '
        f"THEN CUME_DIST() OVER ("
        f'PARTITION BY ("{col}" IS NOT NULL) '
        f'ORDER BY "{col}"'
        f") END "
        f'AS "{alias}"'
    )


def _validate_inputs(
    con: duckdb.DuckDBPyConnection,
    table_names: list[str],
    require_variant: bool = False,
) -> list[_ScoreInfo]:
    """Check every table exists, is a deduped score, and return metadata."""
    infos: list[_ScoreInfo] = []
    for tbl in table_names:
        row = con.execute(
            "SELECT table_type, deduped, source_column, analysis_level "
            "FROM metadata WHERE table_name = ?",
            [tbl],
        ).fetchone()
        if row is None:
            raise ValueError(f"Table {tbl!r} not found in metadata.")
        ttype, deduped, score_name, analysis_level = row
        if ttype != "score":
            raise ValueError(f"Table {tbl!r} is type {ttype!r}, not 'score'.")
        if not deduped:
            raise ValueError(
                f"Table {tbl!r} has not been deduped. "
                f"Run remove_duplicates first."
            )
        if require_variant and analysis_level != "variant":
            raise ValueError(
                f"Table {tbl!r} has analysis_level {analysis_level!r}. "
                f"Pairwise operations require all tables to be "
                f"analysis_level 'variant'."
            )
        infos.append(_ScoreInfo(tbl, score_name, analysis_level))
    return infos


def _build_subselect(info: _ScoreInfo, add_percentile: bool) -> str:
    """Build a sub-select that renames ``score`` to the score's own name."""
    quoted = f'"{info.table_name}"'
    cols = ", ".join(KEY_COLS)
    parts = [f"SELECT {cols}, score AS \"{info.score_name}\""]
    if add_percentile:
        pct = _percentile_expr(
            info.score_name, f"{info.score_name}_percentile"
        )
        parts[0] += f", {pct}"
    parts.append(f"FROM {quoted}")
    return "\n".join(parts)



def _non_key_columns(subselect: str) -> list[str]:
    """Extract non-key column names from a sub-select by inspecting its text."""
    key_names = {"chrom", "pos", "ref", "alt", "ensg", "key"}
    cols = []
    for m in re.finditer(r'AS\s+"([^"]+)"', subselect):
        name = m.group(1)
        if name not in key_names:
            cols.append(name)
    return cols


def _merge_intersection_union(
    con: duckdb.DuckDBPyConnection,
    infos: list[_ScoreInfo],
    output_table: str,
    how: str,
    percentile: str,
) -> None:
    """Merge via intersection or union with pre/post/none percentile."""
    add_pre = percentile == "pre"

    subselects = [_build_subselect(info, add_pre) for info in infos]
    aliases = [f"t{i}" for i in range(len(infos))]

    join_type = "INNER" if how == "intersection" else "FULL OUTER"

    # For intersection, the first table's keys are sufficient.
    # For union, we need to COALESCE keys across all tables.
    first = aliases[0]
    key_exprs: list[str] = []
    if join_type == "FULL OUTER":
        for k in KEY_COLS:
            coalesced = ", ".join(f"{a}.{k}" for a in aliases)
            key_exprs.append(f"COALESCE({coalesced}) AS {k}")
    else:
        for k in KEY_COLS:
            key_exprs.append(f"{first}.{k}")

    score_col_refs: list[str] = []
    for i, info in enumerate(infos):
        a = aliases[i]
        non_key = _non_key_columns(subselects[i])
        for c in non_key:
            score_col_refs.append(f'{a}."{c}" AS "{c}"')

    select_clause = ", ".join(key_exprs + score_col_refs)

    from_clause = f"({subselects[0]}) AS {aliases[0]}"
    for i in range(1, len(subselects)):
        from_clause += (
            f'\n    {join_type} JOIN ({subselects[i]}) AS {aliases[i]} '
            f'ON {first}."key" = {aliases[i]}."key"'
        )

    quoted_output = f'"{output_table}"'

    if percentile == "post":
        # Join first, then add percentile columns
        temp_name = f"_tmp_{output_table}"
        quoted_tmp = f'"{temp_name}"'
        con.execute(
            f"CREATE TABLE {quoted_tmp} AS\n"
            f"SELECT {select_clause}\nFROM {from_clause}"
        )
        # Build percentile expressions over the joined result
        pct_exprs = []
        for info in infos:
            pct_exprs.append(
                _percentile_expr(
                    info.score_name, f"{info.score_name}_percentile"
                )
            )
        pct_select = ", ".join(["*"] + pct_exprs)
        con.execute(
            f"CREATE TABLE {quoted_output} AS\n"
            f"SELECT {pct_select} FROM {quoted_tmp}"
        )
        con.execute(f"DROP TABLE {quoted_tmp}")
    else:
        # pre or none: percentile columns (if any) are already in subselects
        con.execute(
            f"CREATE TABLE {quoted_output} AS\n"
            f"SELECT {select_clause}\nFROM {from_clause}"
        )


def _merge_pairwise(
    con: duckdb.DuckDBPyConnection,
    infos: list[_ScoreInfo],
    anchor_table: str,
    output_table: str,
) -> None:
    """Merge via pairwise intersection-based percentiles."""
    anchor_info = next(i for i in infos if i.table_name == anchor_table)
    non_anchor_infos = [i for i in infos if i.table_name != anchor_table]
    anchor_name = anchor_info.score_name
    quoted_anchor = f'"{anchor_table}"'

    pair_tables: list[str] = []
    pair_columns: list[list[str]] = []

    for na_info in non_anchor_infos:
        na_name = na_info.score_name
        quoted_na = f'"{na_info.table_name}"'
        raw_anchor_col = f"_raw_anchor_{na_name}"
        pw_anchor_col = f"{anchor_name}_pairwise_{na_name}"
        pw_na_col = f"{na_name}_pairwise_{anchor_name}"

        pair_tmp = f"_pw_{na_info.table_name}"
        quoted_pair = f'"{pair_tmp}"'

        con.execute(f"""
            CREATE OR REPLACE TABLE {quoted_pair} AS
            WITH intersection AS (
                SELECT
                    a."key",
                    a.score AS anchor_score,
                    b.score AS nonanchor_score
                FROM {quoted_anchor} a
                INNER JOIN {quoted_na} b ON a."key" = b."key"
            )
            SELECT
                "key",
                anchor_score AS "{raw_anchor_col}",
                nonanchor_score AS "{na_name}",
                CASE WHEN anchor_score IS NOT NULL
                     THEN CUME_DIST() OVER (
                         PARTITION BY (anchor_score IS NOT NULL)
                         ORDER BY anchor_score
                     )
                END AS "{pw_anchor_col}",
                CASE WHEN nonanchor_score IS NOT NULL
                     THEN CUME_DIST() OVER (
                         PARTITION BY (nonanchor_score IS NOT NULL)
                         ORDER BY nonanchor_score
                     )
                END AS "{pw_na_col}"
            FROM intersection;
        """)
        pair_tables.append(pair_tmp)
        pair_columns.append(
            [raw_anchor_col, na_name, pw_anchor_col, pw_na_col]
        )

    # Join all pair tables together on key, pulling keys from the first pair
    # and taking the non-key columns from each
    quoted_output = f'"{output_table}"'

    if len(pair_tables) == 1:
        pt = f'"{pair_tables[0]}"'
        # Need to add full key columns from one of the source tables
        con.execute(f"""
            CREATE TABLE {quoted_output} AS
            SELECT
                s.chrom, s.pos, s.ref, s.alt, s.ensg, p."key",
                p.*  EXCLUDE ("key")
            FROM {pt} p
            LEFT JOIN {quoted_anchor} s ON p."key" = s."key";
        """)
    else:
        # Build sequential FULL OUTER JOIN of pair tables
        aliases = [f"p{i}" for i in range(len(pair_tables))]
        first_a = aliases[0]

        key_coalesce = ", ".join(
            f'COALESCE({", ".join(f"{a}."+"\"key\"" for a in aliases)}) AS "key"'
            for _ in [0]
        )

        score_refs: list[str] = []
        for i, cols in enumerate(pair_columns):
            a = aliases[i]
            for c in cols:
                score_refs.append(f'{a}."{c}"')

        select = f"{key_coalesce}, " + ", ".join(score_refs)

        from_part = f'"{pair_tables[0]}" AS {aliases[0]}'
        for i in range(1, len(pair_tables)):
            from_part += (
                f'\nFULL OUTER JOIN "{pair_tables[i]}" AS {aliases[i]} '
                f'ON {first_a}."key" = {aliases[i]}."key"'
            )

        con.execute(f"""
            CREATE TABLE {quoted_output} AS
            SELECT
                s.chrom, s.pos, s.ref, s.alt, s.ensg, merged."key",
                merged.* EXCLUDE ("key")
            FROM (SELECT {select} FROM {from_part}) AS merged
            LEFT JOIN {quoted_anchor} s ON merged."key" = s."key";
        """)

    # Clean up temp pair tables
    for pt in pair_tables:
        con.execute(f'DROP TABLE "{pt}"')


def merge_scores(
    db_path: str,
    table_names: list[str],
    output_table: str,
    set_operation: str,
    percentile: str,
    anchor_table: str | None = None,
) -> None:
    """Merge multiple score tables into a single wide output table.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    table_names : list[str]
        Names of the input score tables to merge.
    output_table : str
        Name for the output merged table.
    set_operation : str
        ``"intersection"``, ``"union"``, or ``"pairwise"``.
    percentile : str
        ``"pre"``, ``"post"``, or ``"none"``.
    anchor_table : str or None
        Required for pairwise; must be one of *table_names*.
    """
    if set_operation not in SET_OPERATIONS:
        raise ValueError(
            f"Unknown set_operation {set_operation!r}. "
            f"Choose from: {SET_OPERATIONS}"
        )
    if percentile not in PERCENTILE_MODES:
        raise ValueError(
            f"Unknown percentile mode {percentile!r}. "
            f"Choose from: {PERCENTILE_MODES}"
        )
    if set_operation == "pairwise" and percentile != "none":
        raise ValueError(
            "Pairwise mode computes its own percentiles. "
            "Set --percentile none when using --set_operation pairwise."
        )
    if set_operation == "pairwise" and anchor_table is None:
        raise ValueError(
            "--anchor_table is required for pairwise set operation."
        )
    if anchor_table is not None and set_operation != "pairwise":
        raise ValueError(
            "--anchor_table is only used with --set_operation pairwise."
        )
    if anchor_table is not None and anchor_table not in table_names:
        raise ValueError(
            f"anchor_table {anchor_table!r} must be one of --tables."
        )
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

        infos = _validate_inputs(
            con, table_names,
            require_variant=(set_operation == "pairwise"),
        )

        if set_operation in ("intersection", "union"):
            _merge_intersection_union(
                con, infos, output_table, set_operation, percentile,
            )
        else:
            _merge_pairwise(con, infos, anchor_table, output_table)

        quoted_output = f'"{output_table}"'
        row_count = con.execute(
            f"SELECT COUNT(*) FROM {quoted_output}"
        ).fetchone()[0]
        col_count = len(
            con.execute(
                f"SELECT * FROM {quoted_output} LIMIT 0"
            ).description
        )

        source_names = ", ".join(i.score_name for i in infos)
        analysis_levels = {i.analysis_level for i in infos}
        merged_analysis_level = analysis_levels.pop() if len(analysis_levels) == 1 else "variant"
        con.execute(
            "INSERT INTO metadata "
            "(source_column, source_path, table_name, table_type, analysis_level, deduped) "
            "VALUES (?, ?, ?, 'merged_scores', ?, TRUE)",
            [source_names, set_operation, output_table, merged_analysis_level],
        )

        input_tables = ", ".join(i.table_name for i in infos)
        log_event(
            con, "merge_scores", "create_table", output_table,
            f"set_operation={set_operation}, percentile={percentile}, "
            f"input_tables=[{input_tables}], rows={row_count}",
        )

        print(
            f"Merged {len(infos)} tables ({set_operation}, "
            f"percentile={percentile}) → {output_table!r}: "
            f"{row_count:,} rows, {col_count} columns",
            file=sys.stderr,
        )
    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge multiple score tables into a single wide table.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--tables", required=True, nargs="+",
        help="Names of the input score tables to merge.",
    )
    parser.add_argument(
        "--output_table", required=True,
        help="Name for the output merged table.",
    )
    parser.add_argument(
        "--set_operation", required=True, choices=list(SET_OPERATIONS),
        help="Set operation: intersection, union, or pairwise.",
    )
    parser.add_argument(
        "--percentile", required=True, choices=list(PERCENTILE_MODES),
        help="Percentile mode: pre, post, or none.",
    )
    parser.add_argument(
        "--anchor_table", default=None,
        help="Anchor table for pairwise mode (must be one of --tables).",
    )
    args = parser.parse_args()
    merge_scores(
        args.db, args.tables, args.output_table,
        args.set_operation, args.percentile, args.anchor_table,
    )


if __name__ == "__main__":
    main()
