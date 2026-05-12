#!/usr/bin/env python3
"""Join a merged_scores table and a merged_evals table into a merged_analysis table."""
from __future__ import annotations

import argparse
import os
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

KEY_COLS = ("chrom", "pos", "ref", "alt", "ensg", '"key"')
KEY_NAMES = {"chrom", "pos", "ref", "alt", "ensg", "key"}


def _data_columns(
    con: duckdb.DuckDBPyConnection,
    table_name: str,
) -> list[str]:
    """Return non-key column names from a table."""
    cols = [
        desc[0]
        for desc in con.execute(
            f'SELECT * FROM "{table_name}" LIMIT 0'
        ).description
    ]
    return [c for c in cols if c not in KEY_NAMES]


def _validate_table(
    con: duckdb.DuckDBPyConnection,
    table_name: str,
    expected_type: str | tuple[str, ...],
) -> str:
    """Validate a table exists with the expected type. Return its analysis_level."""
    row = con.execute(
        "SELECT DISTINCT table_type, analysis_level FROM metadata "
        "WHERE table_name = ?",
        [table_name],
    ).fetchone()
    if row is None:
        raise ValueError(f"Table {table_name!r} not found in metadata.")
    ttype, analysis_level = row
    allowed = (expected_type,) if isinstance(expected_type, str) else expected_type
    if ttype not in allowed:
        raise ValueError(
            f"Table {table_name!r} is type {ttype!r}, expected one of {allowed!r}."
        )
    return analysis_level


def _same_key_join(
    con: duckdb.DuckDBPyConnection,
    scores_table: str,
    evals_table: str,
    output_table: str,
    *,
    score_fields: list[str] | None = None,
    eval_fields: list[str] | None = None,
    drop_null_scores: bool = False,
) -> None:
    """FULL OUTER JOIN two tables that share the same analysis_level."""
    scores_data = score_fields or _data_columns(con, scores_table)
    evals_data = eval_fields or _data_columns(con, evals_table)

    key_coalesce = ", ".join(
        f'COALESCE(s.{k}, e.{k}) AS {k}' for k in KEY_COLS
    )
    score_refs = ", ".join(f's."{c}"' for c in scores_data)
    eval_refs = ", ".join(f'e."{c}"' for c in evals_data)

    select_parts = [key_coalesce]
    if score_refs:
        select_parts.append(score_refs)
    if eval_refs:
        select_parts.append(eval_refs)
    select_clause = ", ".join(select_parts)

    where_clause = ""
    if drop_null_scores:
        conditions = " AND ".join(f's."{c}" IS NOT NULL' for c in scores_data)
        where_clause = f"WHERE {conditions}"

    con.execute(f"""
        CREATE TABLE "{output_table}" AS
        SELECT {select_clause}
        FROM "{scores_table}" s
        FULL OUTER JOIN "{evals_table}" e ON s."key" = e."key"
        {where_clause}
    """)


def _cross_key_join(
    con: duckdb.DuckDBPyConnection,
    scores_table: str,
    evals_table: str,
    output_table: str,
    linker_path: str,
    scores_level: str,
    *,
    score_fields: list[str] | None = None,
    eval_fields: list[str] | None = None,
    drop_null_scores: bool = False,
) -> None:
    """Join variant-level and gene-level tables via a linker Parquet."""
    parquet_columns = {
        r[0]
        for r in con.execute(
            f"SELECT name FROM parquet_schema('{linker_path}')"
        ).fetchall()
    }
    required = {"chrom", "pos", "ref", "alt", "ensg"}
    missing = required - parquet_columns
    if missing:
        raise ValueError(
            f"Linker Parquet is missing required columns: {sorted(missing)}"
        )

    if scores_level == "variant":
        variant_table, gene_table = scores_table, evals_table
        variant_is_scores = True
    else:
        variant_table, gene_table = evals_table, scores_table
        variant_is_scores = False

    variant_data = (
        (score_fields if variant_is_scores else eval_fields)
        or _data_columns(con, variant_table)
    )
    gene_data = (
        (eval_fields if variant_is_scores else score_fields)
        or _data_columns(con, gene_table)
    )

    con.execute(f"""
        CREATE TEMP TABLE _linker AS
        SELECT
            chrom, pos::BIGINT AS pos, ref, alt, ensg,
            hash(chrom || '|' || CAST(pos AS VARCHAR)
                 || '|' || ref || '|' || alt) AS "key"
        FROM read_parquet('{linker_path}')
    """)

    vt_data_refs = ", ".join(f'vt."{c}"' for c in variant_data)
    gt_data_refs = ", ".join(f'gt."{c}"' for c in gene_data)

    key_select = (
        'COALESCE(lk.chrom, vt.chrom) AS chrom, '
        'COALESCE(lk.pos, vt.pos) AS pos, '
        'COALESCE(lk.ref, vt.ref) AS ref, '
        'COALESCE(lk.alt, vt.alt) AS alt, '
        'COALESCE(lk.ensg, gt.ensg) AS ensg, '
        'COALESCE(lk."key", vt."key") AS "key"'
    )

    select_parts = [key_select]
    if variant_is_scores:
        if vt_data_refs:
            select_parts.append(vt_data_refs)
        if gt_data_refs:
            select_parts.append(gt_data_refs)
    else:
        if gt_data_refs:
            select_parts.append(gt_data_refs)
        if vt_data_refs:
            select_parts.append(vt_data_refs)

    select_clause = ", ".join(select_parts)

    where_clause = ""
    if drop_null_scores:
        score_cols = score_fields or (
            variant_data if variant_is_scores else gene_data
        )
        conditions = " AND ".join(
            f'vt."{c}" IS NOT NULL' if variant_is_scores
            else f'gt."{c}" IS NOT NULL'
            for c in score_cols
        )
        where_clause = f"WHERE {conditions}"

    con.execute(f"""
        CREATE TABLE "{output_table}" AS
        SELECT {select_clause}
        FROM _linker lk
        LEFT JOIN "{variant_table}" vt ON lk."key" = vt."key"
        LEFT JOIN "{gene_table}" gt ON lk.ensg = gt.ensg
        {where_clause}
    """)

    con.execute("DROP TABLE _linker")


def create_analysis_table(
    db_path: str,
    scores_table: str,
    evals_table: str,
    output_table: str,
    linker_path: str | None = None,
    *,
    score_fields: list[str] | None = None,
    eval_fields: list[str] | None = None,
    drop_null_scores: bool = False,
) -> None:
    """Join a merged_scores and a merged_evals table into a merged_analysis table.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    scores_table : str
        Name of a ``merged_scores`` or ``score`` table.
    evals_table : str
        Name of a ``merged_evals`` or ``eval`` table.
    output_table : str
        Name for the output merged analysis table.
    linker_path : str or None
        Path to a linker Parquet file (required when the two input tables
        have different ``analysis_level`` values).
    score_fields : list[str] or None
        Optional subset of score columns to include.  If None, all non-key
        columns from the scores table are included.
    eval_fields : list[str] or None
        Optional subset of eval columns to include.  If None, all non-key
        columns from the evals table are included.
    drop_null_scores : bool
        If True, rows where any selected score column is NULL are excluded
        before the join (intersection semantics for scores).
    """
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

        scores_level = _validate_table(
            con, scores_table, ("merged_scores", "score")
        )
        evals_level = _validate_table(
            con, evals_table, ("merged_evals", "eval")
        )

        if score_fields:
            available = set(_data_columns(con, scores_table))
            bad = [f for f in score_fields if f not in available]
            if bad:
                raise ValueError(
                    f"--score_fields not found in {scores_table!r}: {bad}"
                )
        if eval_fields:
            available = set(_data_columns(con, evals_table))
            bad = [f for f in eval_fields if f not in available]
            if bad:
                raise ValueError(
                    f"--eval_fields not found in {evals_table!r}: {bad}"
                )

        same_level = scores_level == evals_level

        if not same_level and linker_path is None:
            raise ValueError(
                f"Scores table {scores_table!r} is {scores_level!r} and "
                f"evals table {evals_table!r} is {evals_level!r}. "
                f"A --linker_path is required to join tables with "
                f"different analysis levels."
            )
        if same_level and linker_path is not None:
            raise ValueError(
                f"Both tables share analysis_level {scores_level!r}. "
                f"A --linker_path is not needed and should not be provided."
            )

        if linker_path is not None:
            linker_path = os.path.abspath(linker_path)
            if not os.path.isfile(linker_path):
                raise FileNotFoundError(
                    f"Linker Parquet not found: {linker_path}"
                )

        if same_level:
            _same_key_join(
                con, scores_table, evals_table, output_table,
                score_fields=score_fields,
                eval_fields=eval_fields,
                drop_null_scores=drop_null_scores,
            )
            output_level = scores_level
        else:
            _cross_key_join(
                con, scores_table, evals_table, output_table,
                linker_path, scores_level,
                score_fields=score_fields,
                eval_fields=eval_fields,
                drop_null_scores=drop_null_scores,
            )
            output_level = "variant"

        quoted_output = f'"{output_table}"'
        row_count = con.execute(
            f"SELECT COUNT(*) FROM {quoted_output}"
        ).fetchone()[0]
        col_count = len(
            con.execute(
                f"SELECT * FROM {quoted_output} LIMIT 0"
            ).description
        )

        scores_meta_rows = con.execute(
            "SELECT source_column FROM metadata WHERE table_name = ?",
            [scores_table],
        ).fetchall()
        scores_meta = ", ".join(r[0] for r in scores_meta_rows)
        evals_meta_rows = con.execute(
            "SELECT source_column FROM metadata WHERE table_name = ?",
            [evals_table],
        ).fetchall()
        evals_meta = ", ".join(r[0] for r in evals_meta_rows)
        source_names = f"{scores_meta}; {evals_meta}"

        con.execute(
            "INSERT INTO metadata "
            "(source_column, source_path, table_name, table_type, "
            "analysis_level, deduped) "
            "VALUES (?, ?, ?, 'merged_analysis', ?, TRUE)",
            [source_names, f"{scores_table}+{evals_table}",
             output_table, output_level],
        )

        linker_note = f", linker_path={linker_path}" if linker_path else ""
        drop_note = ", drop_null_scores=True" if drop_null_scores else ""
        log_event(
            con, "create_analysis_table", "create_table", output_table,
            f"scores_table={scores_table}, evals_table={evals_table}, "
            f"analysis_level={output_level}, rows={row_count}"
            f"{linker_note}{drop_note}",
        )

        print(
            f"Created analysis table {output_table!r}: "
            f"{row_count:,} rows, {col_count} columns "
            f"(analysis_level={output_level!r})",
            file=sys.stderr,
        )
    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Join a merged_scores table and a merged_evals table "
            "into a merged_analysis table."
        ),
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--scores_table", required=True,
        help="Name of the merged_scores or score table.",
    )
    parser.add_argument(
        "--evals_table", required=True,
        help="Name of the merged_evals or eval table.",
    )
    parser.add_argument(
        "--output_table", required=True,
        help="Name for the output merged_analysis table.",
    )
    parser.add_argument(
        "--linker_path", default=None,
        help="Path to a linker Parquet (required for cross-key joins).",
    )
    parser.add_argument(
        "--score_fields", nargs="+", default=None,
        help="Optional subset of score columns to include.",
    )
    parser.add_argument(
        "--eval_fields", nargs="+", default=None,
        help="Optional subset of eval columns to include.",
    )
    parser.add_argument(
        "--drop_null_scores", action="store_true",
        help="Exclude rows where any selected score column is NULL "
             "(intersection semantics).",
    )
    args = parser.parse_args()
    create_analysis_table(
        args.db, args.scores_table, args.evals_table,
        args.output_table, args.linker_path,
        score_fields=args.score_fields,
        eval_fields=args.eval_fields,
        drop_null_scores=args.drop_null_scores,
    )


if __name__ == "__main__":
    main()
