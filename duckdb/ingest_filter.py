#!/usr/bin/env python3
"""Ingest a variant-level or gene-level filter into a wide DuckDB filter table.

Two ingestion modes are supported:

1. **Membership mode** (default): every row present in the source file is
   treated as TRUE for the filter.  The file is collapsed to distinct keys.

2. **Column mode** (``--source_column``): a specific boolean column in the
   source file provides the filter value.  Only rows where that column is
   TRUE are marked; all others are NULL.  This is useful for pre-pivoted
   wide filter files that already contain multiple boolean columns.

Variant filter files are expected to have at minimum (chrom, pos, ref, alt).
Gene filter files require an (ensg) column.
"""
from __future__ import annotations

import argparse
import os
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

WIDE_TABLE_DEFAULTS = {"variant": "variant_filters", "gene": "gene_filters"}
KEY_NAMES = {"chrom", "pos", "ref", "alt", "ensg", "key"}


def ingest_filter(
    db_path: str,
    filter_name: str,
    filter_path: str,
    table_name: str,
    analysis_level: str,
    *,
    source_column: str | None = None,
) -> None:
    """Add a boolean filter column to a wide filter table.

    If the wide table does not exist, it is created with canonical key columns
    plus the first filter column.  If it already exists, the new filter column
    is added (ALTER TABLE) and populated via UPDATE/INSERT.

    The source file is collapsed to distinct keys before insertion — for
    variant-level data this means deduplication from variant-transcript
    granularity down to (chrom, pos, ref, alt).

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file (must already be initialised).
    filter_name : str
        Name for the boolean column in the wide table (e.g. ``"buried"``).
    filter_path : str
        Path to the source TSV/CSV/Parquet file.
    table_name : str
        Name for the wide filter table (e.g. ``"variant_filters"``).
    analysis_level : str
        ``"variant"`` or ``"gene"``.
    source_column : str or None
        If provided, read this boolean column from the source file instead of
        treating all rows as TRUE (column mode).  Only rows where the column
        is TRUE are included.
    """
    filter_path = os.path.abspath(filter_path)

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Database not found: {db_path}")
    if not os.path.isfile(filter_path):
        raise FileNotFoundError(f"Source file not found: {filter_path}")
    if analysis_level not in ("variant", "gene"):
        raise ValueError(
            f"analysis_level must be 'variant' or 'gene', got {analysis_level!r}"
        )

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

        read_fn = _read_expr(filter_path)

        source_columns = {
            d[0]
            for d in con.execute(
                f"SELECT * FROM {read_fn} LIMIT 0"
            ).description
        }

        if source_column and source_column not in source_columns:
            raise ValueError(
                f"--source_column {source_column!r} not found in source. "
                f"Available: {sorted(source_columns)}"
            )

        if analysis_level == "variant":
            missing = [
                k for k in ("chrom", "pos", "ref", "alt")
                if k not in source_columns
            ]
            if missing:
                raise ValueError(
                    f"Variant key column(s) missing from source: {missing}"
                )
            key_select = (
                "chrom, "
                "pos::BIGINT AS pos, "
                "ref, "
                "alt, "
                "NULL::VARCHAR AS ensg"
            )
            hash_expr = (
                "hash(chrom || '|' || CAST(pos AS VARCHAR) "
                "|| '|' || ref || '|' || alt)"
            )
        else:
            if "ensg" not in source_columns:
                raise ValueError("Key column 'ensg' not found in source.")
            key_select = (
                "NULL::VARCHAR AS chrom, "
                "NULL::BIGINT AS pos, "
                "NULL::VARCHAR AS ref, "
                "NULL::VARCHAR AS alt, "
                "ensg"
            )
            hash_expr = "hash(ensg)"

        quoted_table = f'"{table_name}"'
        quoted_filter = f'"{filter_name}"'

        where_filter = (
            f'WHERE "{source_column}" = TRUE' if source_column else ""
        )
        deduped_src = f"""
            (SELECT DISTINCT {key_select}, {hash_expr} AS "key"
             FROM {read_fn} {where_filter})
        """

        table_exists = table_name in tables

        if not table_exists:
            con.execute(
                f"CREATE TABLE {quoted_table} AS\n"
                f"SELECT chrom, pos, ref, alt, ensg, \"key\", "
                f"TRUE AS {quoted_filter}\n"
                f"FROM {deduped_src}"
            )
            row_count = con.execute(
                f"SELECT COUNT(*) FROM {quoted_table}"
            ).fetchone()[0]
            print(
                f"Created {table_name!r} with filter {filter_name!r}: "
                f"{row_count:,} rows (all TRUE)",
                file=sys.stderr,
            )
        else:
            existing_cols = {
                desc[0]
                for desc in con.execute(
                    f"SELECT * FROM {quoted_table} LIMIT 0"
                ).description
            }
            if filter_name in existing_cols:
                raise ValueError(
                    f"Filter column {filter_name!r} already exists in "
                    f"{table_name!r}."
                )

            con.execute(
                f"ALTER TABLE {quoted_table} ADD COLUMN {quoted_filter} BOOLEAN"
            )

            con.execute(f"""
                UPDATE {quoted_table} AS t
                SET {quoted_filter} = TRUE
                FROM {deduped_src} AS src
                WHERE t."key" = src."key"
            """)

            con.execute(f"""
                INSERT INTO {quoted_table}
                SELECT src.chrom, src.pos, src.ref, src.alt, src.ensg,
                       src."key", {_null_fills(con, table_name, filter_name)},
                       TRUE AS {quoted_filter}
                FROM {deduped_src} AS src
                WHERE src."key" NOT IN (SELECT "key" FROM {quoted_table})
            """)

            row_count = con.execute(
                f"SELECT COUNT(*) FROM {quoted_table}"
            ).fetchone()[0]
            true_count = con.execute(
                f"SELECT COUNT({quoted_filter}) FROM {quoted_table}"
            ).fetchone()[0]
            print(
                f"Added filter {filter_name!r} to {table_name!r}: "
                f"{row_count:,} total rows ({true_count:,} TRUE, "
                f"{row_count - true_count:,} NULL)",
                file=sys.stderr,
            )

        con.execute(
            "INSERT INTO metadata "
            "(source_column, source_path, table_name, table_type, "
            "analysis_level, deduped) "
            "VALUES (?, ?, ?, 'filter', ?, TRUE)",
            [filter_name, filter_path, table_name, analysis_level],
        )
        log_event(
            con, "ingest_filter", "add_column", table_name,
            f"filter_name={filter_name}, source_path={filter_path}, "
            f"analysis_level={analysis_level}, rows={row_count}",
        )

    finally:
        con.close()


def _read_expr(path: str) -> str:
    """Return the DuckDB read expression appropriate for the file type."""
    lower = path.lower()
    if lower.endswith(".parquet"):
        return f"read_parquet('{path}')"
    return f"read_csv_auto('{path}', delim='\\t')"


def _null_fills(
    con: duckdb.DuckDBPyConnection,
    table_name: str,
    new_filter: str,
) -> str:
    """Build NULL fill expressions for existing filter columns (excluding keys
    and the new filter) for the INSERT of new rows."""
    cols = [
        desc[0]
        for desc in con.execute(
            f'SELECT * FROM "{table_name}" LIMIT 0'
        ).description
    ]
    fills = []
    for c in cols:
        if c in KEY_NAMES or c == new_filter:
            continue
        fills.append(f'NULL AS "{c}"')
    return ", ".join(fills)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ingest a filter column into a wide DuckDB filter table.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--filter_name", required=True,
        help="Name for the boolean filter column (e.g. 'buried').",
    )
    parser.add_argument(
        "--filter_path", required=True,
        help="Path to the source TSV/Parquet file.",
    )
    parser.add_argument(
        "--table_name", default=None,
        help="Name for the wide filter table "
             "(default: variant_filters or gene_filters).",
    )
    parser.add_argument(
        "--analysis_level", required=True, choices=["variant", "gene"],
        help="Key type: 'variant' (chrom/pos/ref/alt) or 'gene' (ensg).",
    )
    parser.add_argument(
        "--source_column", default=None,
        help="Read this boolean column from the source instead of treating "
             "all rows as TRUE (for pre-pivoted wide filter files).",
    )
    args = parser.parse_args()
    table_name = args.table_name or WIDE_TABLE_DEFAULTS[args.analysis_level]
    ingest_filter(args.db, args.filter_name, args.filter_path,
                  table_name, args.analysis_level,
                  source_column=args.source_column)


if __name__ == "__main__":
    main()
