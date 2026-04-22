#!/usr/bin/env python3
"""Join a merged_scores table and a merged_evals table into a merged_analysis table."""
from __future__ import annotations

import argparse
import os
import sys

import duckdb


DEFAULT_DB_PATH = "gg_data.duckdb"

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
    expected_type: str,
) -> str:
    """Validate a table exists with the expected type. Return its analysis_level."""
    row = con.execute(
        "SELECT table_type, analysis_level FROM metadata WHERE table_name = ?",
        [table_name],
    ).fetchone()
    if row is None:
        raise ValueError(f"Table {table_name!r} not found in metadata.")
    ttype, analysis_level = row
    if ttype != expected_type:
        raise ValueError(
            f"Table {table_name!r} is type {ttype!r}, expected {expected_type!r}."
        )
    return analysis_level


def _same_key_join(
    con: duckdb.DuckDBPyConnection,
    scores_table: str,
    evals_table: str,
    output_table: str,
) -> None:
    """FULL OUTER JOIN two tables that share the same analysis_level."""
    scores_data = _data_columns(con, scores_table)
    evals_data = _data_columns(con, evals_table)

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

    con.execute(f"""
        CREATE TABLE "{output_table}" AS
        SELECT {select_clause}
        FROM "{scores_table}" s
        FULL OUTER JOIN "{evals_table}" e ON s."key" = e."key"
    """)


def _cross_key_join(
    con: duckdb.DuckDBPyConnection,
    scores_table: str,
    evals_table: str,
    output_table: str,
    linker_path: str,
    scores_level: str,
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

    # Determine which input is variant-level and which is gene-level
    if scores_level == "variant":
        variant_table, gene_table = scores_table, evals_table
        variant_is_scores = True
    else:
        variant_table, gene_table = evals_table, scores_table
        variant_is_scores = False

    variant_data = _data_columns(con, variant_table)
    gene_data = _data_columns(con, gene_table)

    # Load linker into a temp table with a computed variant key hash
    con.execute(f"""
        CREATE TEMP TABLE _linker AS
        SELECT
            chrom, pos::BIGINT AS pos, ref, alt, ensg,
            hash(chrom || '|' || CAST(pos AS VARCHAR)
                 || '|' || ref || '|' || alt) AS "key"
        FROM read_parquet('{linker_path}')
    """)

    # Build the three-way join:
    # variant_table JOIN linker ON variant key -> linker provides ensg
    # linker+variant JOIN gene_table ON ensg
    vt_data_refs = ", ".join(f'vt."{c}"' for c in variant_data)
    gt_data_refs = ", ".join(f'gt."{c}"' for c in gene_data)

    # Key columns come from the linker (has both variant and gene keys)
    # coalesce with the variant table for coverage
    key_select = (
        'COALESCE(lk.chrom, vt.chrom) AS chrom, '
        'COALESCE(lk.pos, vt.pos) AS pos, '
        'COALESCE(lk.ref, vt.ref) AS ref, '
        'COALESCE(lk.alt, vt.alt) AS alt, '
        'COALESCE(lk.ensg, gt.ensg) AS ensg, '
        'COALESCE(lk."key", vt."key") AS "key"'
    )

    select_parts = [key_select]
    # Add score columns then eval columns, regardless of which is variant/gene
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

    con.execute(f"""
        CREATE TABLE "{output_table}" AS
        SELECT {select_clause}
        FROM _linker lk
        LEFT JOIN "{variant_table}" vt ON lk."key" = vt."key"
        LEFT JOIN "{gene_table}" gt ON lk.ensg = gt.ensg
    """)

    con.execute("DROP TABLE _linker")


def create_analysis_table(
    db_path: str,
    scores_table: str,
    evals_table: str,
    output_table: str,
    linker_path: str | None = None,
) -> None:
    """Join a merged_scores and a merged_evals table into a merged_analysis table.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    scores_table : str
        Name of a ``merged_scores`` table.
    evals_table : str
        Name of a ``merged_evals`` table.
    output_table : str
        Name for the output merged analysis table.
    linker_path : str or None
        Path to a linker Parquet file (required when the two input tables
        have different ``analysis_level`` values).
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

        scores_level = _validate_table(con, scores_table, "merged_scores")
        evals_level = _validate_table(con, evals_table, "merged_evals")

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
            _same_key_join(con, scores_table, evals_table, output_table)
            output_level = scores_level
        else:
            _cross_key_join(
                con, scores_table, evals_table, output_table,
                linker_path, scores_level,
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

        scores_meta = con.execute(
            "SELECT source_column FROM metadata WHERE table_name = ?",
            [scores_table],
        ).fetchone()[0]
        evals_meta = con.execute(
            "SELECT source_column FROM metadata WHERE table_name = ?",
            [evals_table],
        ).fetchone()[0]
        source_names = f"{scores_meta}; {evals_meta}"

        con.execute(
            "INSERT INTO metadata "
            "(source_column, source_path, table_name, table_type, "
            "analysis_level, deduped) "
            "VALUES (?, ?, ?, 'merged_analysis', ?, TRUE)",
            [source_names, f"{scores_table}+{evals_table}",
             output_table, output_level],
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
        help="Name of the merged_scores table.",
    )
    parser.add_argument(
        "--evals_table", required=True,
        help="Name of the merged_evals table.",
    )
    parser.add_argument(
        "--output_table", required=True,
        help="Name for the output merged_analysis table.",
    )
    parser.add_argument(
        "--linker_path", default=None,
        help="Path to a linker Parquet (required for cross-key joins).",
    )
    args = parser.parse_args()
    create_analysis_table(
        args.db, args.scores_table, args.evals_table,
        args.output_table, args.linker_path,
    )


if __name__ == "__main__":
    main()
