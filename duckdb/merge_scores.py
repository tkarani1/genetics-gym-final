#!/usr/bin/env python3
"""Merge score columns from a wide score table and export to Parquet.

With the wide-table model every score column lives in the same DuckDB
table (e.g. ``variant_scores``).  A "merge" is now a single-table scan
with an optional WHERE clause rather than an N-way join.

Pairwise mode produces one Parquet file per anchor-target pair, each
containing keys, both raw scores, both pairwise percentiles, and
(optionally) all eval columns pre-joined.
"""
from __future__ import annotations

import argparse
import os
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

SET_OPERATIONS = ("intersection", "union", "pairwise")
PERCENTILE_MODES = ("pre", "post", "none")

KEY_COLS = ("chrom", "pos", "ref", "alt", "ensg", '"key"')
KEY_NAMES = {"chrom", "pos", "ref", "alt", "ensg", "key"}


def _percentile_expr(col: str, alias: str) -> str:
    """SQL for null-safe CUME_DIST over *col*."""
    return (
        f'CASE WHEN "{col}" IS NOT NULL '
        f"THEN CUME_DIST() OVER ("
        f'PARTITION BY ("{col}" IS NOT NULL) '
        f'ORDER BY "{col}"'
        f") END "
        f'AS "{alias}"'
    )


def _validate_columns(
    con: duckdb.DuckDBPyConnection,
    wide_table: str,
    score_columns: list[str],
) -> str:
    """Verify the wide table exists and all requested columns are present.
    Returns the analysis_level for the table."""
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
    missing = [c for c in score_columns if c not in actual_cols]
    if missing:
        raise ValueError(
            f"Score column(s) not found in {wide_table!r}: {missing}"
        )

    row = con.execute(
        "SELECT DISTINCT analysis_level FROM metadata WHERE table_name = ?",
        [wide_table],
    ).fetchone()
    if row is None:
        raise ValueError(f"No metadata found for {wide_table!r}.")
    return row[0]


def _eval_col_refs(
    con: duckdb.DuckDBPyConnection,
    evals_table: str,
) -> list[str]:
    """Return non-key column names from the merged evals table."""
    cols = [
        desc[0]
        for desc in con.execute(
            f'SELECT * FROM "{evals_table}" LIMIT 0'
        ).description
    ]
    return [c for c in cols if c not in KEY_NAMES]


def _merge_intersection_union(
    con: duckdb.DuckDBPyConnection,
    wide_table: str,
    score_columns: list[str],
    output_path: str,
    how: str,
    percentile: str,
) -> int:
    """Intersection or union merge -> Parquet."""
    quoted = f'"{wide_table}"'
    key_select = ", ".join(KEY_COLS)
    score_select = ", ".join(f'"{c}"' for c in score_columns)

    pct_parts: list[str] = []
    if percentile in ("pre", "post"):
        for c in score_columns:
            pct_parts.append(_percentile_expr(c, f"{c}_percentile"))

    select_parts = [key_select, score_select]
    if pct_parts:
        select_parts.append(", ".join(pct_parts))
    select_clause = ", ".join(select_parts)

    if how == "intersection":
        where_clause = " AND ".join(
            f'"{c}" IS NOT NULL' for c in score_columns
        )
        query = f"SELECT {select_clause} FROM {quoted} WHERE {where_clause}"
    else:
        query = f"SELECT {select_clause} FROM {quoted}"

    con.execute(
        f"COPY ({query}) TO '{output_path}' (FORMAT PARQUET)"
    )
    row_count = con.execute(
        f"SELECT COUNT(*) FROM read_parquet('{output_path}')"
    ).fetchone()[0]
    return row_count


def _merge_pairwise(
    con: duckdb.DuckDBPyConnection,
    wide_table: str,
    score_columns: list[str],
    anchor_column: str,
    output_dir: str,
    evals_table: str | None = None,
    linker_path: str | None = None,
    gene_average: bool = False,
) -> list[str]:
    """Pairwise intersection-based percentiles -> one Parquet per pair.

    Each output file contains keys, the anchor score, the target score,
    2 pairwise percentile columns, and (if *evals_table* is given) all
    eval columns pre-joined.

    When *gene_average* is True, the linker is joined to obtain ``ensg``
    and each pairwise percentile column is supplemented with a
    gene-averaged version (``AVG(...) OVER (PARTITION BY ensg)``).
    The per-variant percentiles are kept alongside the gene-averaged ones.

    Returns the list of output file paths produced.
    """
    quoted = f'"{wide_table}"'
    non_anchor = [c for c in score_columns if c != anchor_column]
    n = len(non_anchor)

    eval_select = ""
    eval_join = ""
    if evals_table is not None:
        eval_data_cols = _eval_col_refs(con, evals_table)
        if eval_data_cols:
            eval_select = ", " + ", ".join(
                f'ev."{c}"' for c in eval_data_cols
            )
            eval_join = (
                f'LEFT JOIN "{evals_table}" ev '
                f'ON base."key" = ev."key"'
            )

    suffix = "_gene_avg" if gene_average else ""
    output_files: list[str] = []

    for idx, na in enumerate(non_anchor):
        both_not_null = (
            f'(base."{anchor_column}" IS NOT NULL '
            f'AND base."{na}" IS NOT NULL)'
        )
        pw_anchor = f"{anchor_column}_pairwise_{na}"
        pw_na = f"{na}_pairwise_{anchor_column}"

        out_path = os.path.join(
            output_dir, f"{anchor_column}_x_{na}{suffix}.parquet"
        )
        if os.path.exists(out_path):
            print(
                f"  pair {idx + 1}/{n}: {anchor_column} x {na} "
                f"— SKIPPED (file exists)",
                file=sys.stderr,
            )
            output_files.append(out_path)
            continue

        key_select = ", ".join(f"base.{k}" for k in KEY_COLS)

        if not gene_average:
            query = f"""
                SELECT {key_select},
                    base."{anchor_column}",
                    base."{na}",
                    CASE WHEN {both_not_null} THEN CUME_DIST() OVER (
                        PARTITION BY {both_not_null}
                        ORDER BY base."{anchor_column}"
                    ) END AS "{pw_anchor}",
                    CASE WHEN {both_not_null} THEN CUME_DIST() OVER (
                        PARTITION BY {both_not_null}
                        ORDER BY base."{na}"
                    ) END AS "{pw_na}"
                    {eval_select}
                FROM {quoted} base
                {eval_join}
            """
        else:
            # Two-level CTE: first compute per-variant pairwise
            # percentiles, then join linker and add gene-averaged columns.
            inner_eval_select = eval_select.replace("ev.", "ev.")
            inner_eval_join = eval_join

            # Passthrough eval columns from the CTE
            eval_passthrough = ""
            if evals_table is not None and eval_data_cols:
                eval_passthrough = ", " + ", ".join(
                    f'p."{c}"' for c in eval_data_cols
                )

            query = f"""
                WITH pairwise AS (
                    SELECT {key_select},
                        base."{anchor_column}",
                        base."{na}",
                        CASE WHEN {both_not_null} THEN CUME_DIST() OVER (
                            PARTITION BY {both_not_null}
                            ORDER BY base."{anchor_column}"
                        ) END AS "{pw_anchor}",
                        CASE WHEN {both_not_null} THEN CUME_DIST() OVER (
                            PARTITION BY {both_not_null}
                            ORDER BY base."{na}"
                        ) END AS "{pw_na}"
                        {inner_eval_select}
                    FROM {quoted} base
                    {inner_eval_join}
                ),
                linker_dedup AS (
                    SELECT
                        hash(chrom || '|' || CAST(pos AS VARCHAR)
                             || '|' || ref || '|' || alt) AS "key",
                        FIRST(ensg) AS ensg
                    FROM read_parquet('{linker_path}')
                    GROUP BY hash(chrom || '|' || CAST(pos AS VARCHAR)
                                  || '|' || ref || '|' || alt)
                )
                SELECT p.chrom, p.pos, p.ref, p.alt,
                    COALESCE(lk.ensg, p.ensg) AS ensg,
                    p."key",
                    p."{anchor_column}",
                    p."{na}",
                    p."{pw_anchor}",
                    p."{pw_na}",
                    AVG(p."{pw_anchor}") OVER (
                        PARTITION BY lk.ensg
                    ) AS "{pw_anchor}_gene_avg",
                    AVG(p."{pw_na}") OVER (
                        PARTITION BY lk.ensg
                    ) AS "{pw_na}_gene_avg"
                    {eval_passthrough}
                FROM pairwise p
                LEFT JOIN linker_dedup lk ON p."key" = lk."key"
            """

        con.execute(
            f"COPY ({query}) TO '{out_path}' (FORMAT PARQUET)"
        )

        row_count = con.execute(
            f"SELECT COUNT(*) FROM read_parquet('{out_path}')"
        ).fetchone()[0]
        output_files.append(out_path)

        label = f"{anchor_column} x {na}"
        if gene_average:
            label += " (gene_avg)"
        print(
            f"  pair {idx + 1}/{n}: {label} "
            f"-> {row_count:,} rows",
            file=sys.stderr,
        )

    return output_files


def merge_scores(
    db_path: str,
    wide_table: str,
    score_columns: list[str],
    output_path: str | None,
    set_operation: str,
    percentile: str,
    anchor_column: str | None = None,
    output_dir: str | None = None,
    evals_table: str | None = None,
    linker_path: str | None = None,
    gene_average: bool = False,
    *,
    memory_limit: str | None = None,
) -> None:
    """Merge score columns from a wide table and write to Parquet.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    wide_table : str
        Name of the wide score table (e.g. ``variant_scores``).
    score_columns : list[str]
        Score columns to include in the merge.
    output_path : str or None
        Destination Parquet file path (for intersection/union).
    set_operation : str
        ``"intersection"``, ``"union"``, or ``"pairwise"``.
    percentile : str
        ``"pre"``, ``"post"``, or ``"none"``.
    anchor_column : str or None
        Required for pairwise; must be one of *score_columns*.
    output_dir : str or None
        Output directory for pairwise (one file per pair).
    evals_table : str or None
        Optional merged-evals table to left-join into pairwise outputs.
    linker_path : str or None
        Path to a linker Parquet (required when *gene_average* is True).
    gene_average : bool
        If True, add gene-averaged pairwise percentile columns by
        joining the linker and computing ``AVG(...) OVER (PARTITION BY ensg)``.
    memory_limit : str or None
        Optional DuckDB memory limit (e.g. ``'8GB'``).
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
    if set_operation == "pairwise" and anchor_column is None:
        raise ValueError(
            "--anchor_column is required for pairwise set operation."
        )
    if anchor_column is not None and set_operation != "pairwise":
        raise ValueError(
            "--anchor_column is only used with --set_operation pairwise."
        )
    if anchor_column is not None and anchor_column not in score_columns:
        raise ValueError(
            f"anchor_column {anchor_column!r} must be one of --columns."
        )
    if len(score_columns) < 2:
        raise ValueError("At least two score columns are required.")
    if gene_average and set_operation != "pairwise":
        raise ValueError(
            "--gene_average is only supported with --set_operation pairwise."
        )
    if gene_average and linker_path is None:
        raise ValueError(
            "--linker_path is required when --gene_average is set."
        )
    if linker_path is not None:
        linker_path = os.path.abspath(linker_path)
        if not os.path.isfile(linker_path):
            raise FileNotFoundError(
                f"Linker Parquet not found: {linker_path}"
            )

    if set_operation == "pairwise":
        if output_dir is None:
            raise ValueError(
                "--output_dir is required for pairwise set operation."
            )
        output_dir = os.path.abspath(output_dir)
        os.makedirs(output_dir, exist_ok=True)
    else:
        if output_path is None:
            raise ValueError(
                "--output_path is required for intersection/union."
            )
        output_path = os.path.abspath(output_path)
        if os.path.exists(output_path):
            raise FileExistsError(
                f"Output file already exists: {output_path}"
            )

    con = duckdb.connect(db_path, read_only=True)
    try:
        if memory_limit:
            con.execute(f"SET memory_limit = '{memory_limit}'")
        con.execute("SET preserve_insertion_order = false")

        _validate_columns(con, wide_table, score_columns)

        if set_operation in ("intersection", "union"):
            row_count = _merge_intersection_union(
                con, wide_table, score_columns, output_path,
                set_operation, percentile,
            )
            print(
                f"Merged {len(score_columns)} columns ({set_operation}, "
                f"percentile={percentile}) -> {output_path}: "
                f"{row_count:,} rows",
                file=sys.stderr,
            )
        else:
            output_files = _merge_pairwise(
                con, wide_table, score_columns,
                anchor_column, output_dir,
                evals_table=evals_table,
                linker_path=linker_path,
                gene_average=gene_average,
            )
            print(
                f"Pairwise complete: {len(output_files)} files in "
                f"{output_dir}",
                file=sys.stderr,
            )
    finally:
        con.close()

    con_rw = duckdb.connect(db_path)
    try:
        source_names = ", ".join(score_columns)
        if set_operation == "pairwise":
            log_event(
                con_rw, "merge_scores", "export_pairwise", wide_table,
                f"anchor={anchor_column}, pairs={len(output_files)}, "
                f"output_dir={output_dir}, "
                f"evals_table={evals_table}, "
                f"gene_average={gene_average}",
            )
        else:
            log_event(
                con_rw, "merge_scores", "export_parquet", wide_table,
                f"set_operation={set_operation}, percentile={percentile}, "
                f"columns=[{source_names}], output={output_path}, "
                f"rows={row_count}",
            )
    finally:
        con_rw.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge score columns from a wide table and export to Parquet.",
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
        "--columns", required=True, nargs="+",
        help="Score columns to include in the merge.",
    )
    parser.add_argument(
        "--output_path", default=None,
        help="Destination Parquet file (for intersection/union).",
    )
    parser.add_argument(
        "--output_dir", default=None,
        help="Output directory for pairwise (one file per pair).",
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
        "--anchor_column", default=None,
        help="Anchor column for pairwise mode (must be one of --columns).",
    )
    parser.add_argument(
        "--evals_table", default=None,
        help="Merged-evals table to left-join into pairwise outputs.",
    )
    parser.add_argument(
        "--linker_path", default=None,
        help="Path to a linker Parquet (required with --gene_average).",
    )
    parser.add_argument(
        "--gene_average", action="store_true",
        help="Add gene-averaged pairwise percentile columns via linker.",
    )
    parser.add_argument(
        "--memory_limit", default=None,
        help="DuckDB memory limit (e.g. '8GB'). Default: DuckDB auto.",
    )
    args = parser.parse_args()
    merge_scores(
        args.db, args.wide_table, args.columns,
        args.output_path,
        args.set_operation, args.percentile, args.anchor_column,
        output_dir=args.output_dir,
        evals_table=args.evals_table,
        linker_path=args.linker_path,
        gene_average=args.gene_average,
        memory_limit=args.memory_limit,
    )


if __name__ == "__main__":
    main()
