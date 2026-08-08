#!/usr/bin/env python3
"""Enrich a merged table with key columns from a linker Parquet file."""
from __future__ import annotations

import argparse
import os
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

MERGED_TYPES = ("merged_scores", "merged_evals", "merged_analysis")
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


def _enrich_variant(
    con: duckdb.DuckDBPyConnection,
    table_name: str,
) -> None:
    """Add ensg to a variant-level table via the _linker temp table.

    A variant may map to multiple genes in the linker.  To keep the row
    count stable we pick one arbitrary ensg per variant key.
    """
    con.execute("""
        CREATE TEMP TABLE _linker_dedup AS
        SELECT "key", FIRST(ensg) AS ensg
        FROM _linker
        GROUP BY "key"
    """)

    data_cols = _data_columns(con, table_name)
    data_refs = ", ".join(f't."{c}"' for c in data_cols)

    select_parts = [
        't.chrom',
        't.pos',
        't.ref',
        't.alt',
        'COALESCE(ld.ensg, t.ensg) AS ensg',
        't."key"',
    ]
    if data_refs:
        select_parts.append(data_refs)
    select_clause = ", ".join(select_parts)

    con.execute(f"""
        CREATE OR REPLACE TABLE "{table_name}" AS
        SELECT {select_clause}
        FROM "{table_name}" t
        LEFT JOIN _linker_dedup ld ON t."key" = ld."key"
    """)

    con.execute("DROP TABLE _linker_dedup")


def _enrich_gene(
    con: duckdb.DuckDBPyConnection,
    table_name: str,
) -> None:
    """Add chrom/pos/ref/alt/key to a gene-level table via the _linker temp table.

    This is a one-to-many join that expands the row count.
    """
    data_cols = _data_columns(con, table_name)
    data_refs = ", ".join(f't."{c}"' for c in data_cols)

    select_parts = [
        'lk.chrom',
        'lk.pos',
        'lk.ref',
        'lk.alt',
        'COALESCE(lk.ensg, t.ensg) AS ensg',
        'lk."key"',
    ]
    if data_refs:
        select_parts.append(data_refs)
    select_clause = ", ".join(select_parts)

    con.execute(f"""
        CREATE OR REPLACE TABLE "{table_name}" AS
        SELECT {select_clause}
        FROM "{table_name}" t
        LEFT JOIN _linker lk ON t.ensg = lk.ensg
    """)


def join_linker(
    db_path: str,
    table_name: str,
    linker_path: str,
) -> None:
    """Enrich a merged table with the 'other side' key columns from a linker.

    For variant-level tables, fills in the ``ensg`` column.
    For gene-level tables, expands rows to variant granularity by populating
    ``chrom``, ``pos``, ``ref``, ``alt``, and ``key``.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    table_name : str
        Name of a merged table (merged_scores, merged_evals, or merged_analysis).
    linker_path : str
        Path to a linker Parquet containing ``chrom``, ``pos``, ``ref``,
        ``alt``, and ``ensg`` columns.
    """
    linker_path = os.path.abspath(linker_path)

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Database not found: {db_path}")
    if not os.path.isfile(linker_path):
        raise FileNotFoundError(f"Linker Parquet not found: {linker_path}")

    con = duckdb.connect(db_path)
    try:
        row = con.execute(
            "SELECT table_type, analysis_level FROM metadata "
            "WHERE table_name = ?",
            [table_name],
        ).fetchone()
        if row is None:
            raise ValueError(f"Table {table_name!r} not found in metadata.")

        table_type, analysis_level = row
        if table_type not in MERGED_TYPES:
            raise ValueError(
                f"Table {table_name!r} is type {table_type!r}. "
                f"join_linker only works on merged tables "
                f"({', '.join(MERGED_TYPES)})."
            )

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
                f"Linker Parquet is missing required columns: "
                f"{sorted(missing)}"
            )

        before_count = con.execute(
            f'SELECT COUNT(*) FROM "{table_name}"'
        ).fetchone()[0]

        con.execute(f"""
            CREATE TEMP TABLE _linker AS
            SELECT
                chrom, pos::BIGINT AS pos, ref, alt, ensg,
                hash(chrom || '|' || CAST(pos AS VARCHAR)
                     || '|' || ref || '|' || alt) AS "key"
            FROM read_parquet('{linker_path}')
        """)

        if analysis_level == "variant":
            _enrich_variant(con, table_name)
            new_level = "variant"
        else:
            _enrich_gene(con, table_name)
            new_level = "variant"
            con.execute(
                "UPDATE metadata SET analysis_level = 'variant' "
                "WHERE table_name = ?",
                [table_name],
            )

        con.execute("DROP TABLE _linker")

        after_count = con.execute(
            f'SELECT COUNT(*) FROM "{table_name}"'
        ).fetchone()[0]

        if analysis_level == "variant":
            ensg_filled = con.execute(
                f'SELECT COUNT(*) FROM "{table_name}" '
                f'WHERE ensg IS NOT NULL'
            ).fetchone()[0]
            log_event(
                con, "join_linker", "enrich_variant", table_name,
                f"linker_path={linker_path}, ensg_filled={ensg_filled}",
            )
            print(
                f"Enriched {table_name!r} (variant): "
                f"{before_count:,} rows, "
                f"{ensg_filled:,} rows now have ensg "
                f"(analysis_level unchanged)",
                file=sys.stderr,
            )
        else:
            log_event(
                con, "join_linker", "enrich_gene_to_variant", table_name,
                f"linker_path={linker_path}, rows_before={before_count}, "
                f"rows_after={after_count}",
            )
            print(
                f"Enriched {table_name!r} (gene -> variant): "
                f"{before_count:,} -> {after_count:,} rows "
                f"(analysis_level changed to 'variant')",
                file=sys.stderr,
            )
    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Enrich a merged table with key columns from a linker "
            "Parquet file."
        ),
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--table_name", required=True,
        help="Name of the merged table to enrich.",
    )
    parser.add_argument(
        "--linker_path", required=True,
        help="Path to a linker Parquet with chrom/pos/ref/alt/ensg.",
    )
    args = parser.parse_args()
    join_linker(args.db, args.table_name, args.linker_path)


if __name__ == "__main__":
    main()
