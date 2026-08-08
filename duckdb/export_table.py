#!/usr/bin/env python3
"""Export a DuckDB table to a Parquet file."""
from __future__ import annotations

import argparse
import os
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event


def export_table(
    db_path: str,
    table_name: str,
    output_path: str,
) -> None:
    """Write the contents of *table_name* to a Parquet file at *output_path*.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    table_name : str
        Name of the table to export.
    output_path : str
        Destination file path for the Parquet output.
    """
    con = duckdb.connect(db_path, read_only=True)
    try:
        exists = con.execute(
            "SELECT COUNT(*) FROM information_schema.tables "
            "WHERE table_name = ? AND table_schema = 'main'",
            [table_name],
        ).fetchone()[0]
        if not exists:
            raise ValueError(f"Table {table_name!r} does not exist.")

        if os.path.exists(output_path):
            raise FileExistsError(
                f"Output file already exists: {output_path}"
            )

        con.execute(
            f"COPY \"{table_name}\" TO '{output_path}' (FORMAT PARQUET)"
        )

        row_count = con.execute(
            f'SELECT COUNT(*) FROM "{table_name}"'
        ).fetchone()[0]

        print(
            f"Exported {table_name!r} → {output_path} "
            f"({row_count:,} rows)",
            file=sys.stderr,
        )
    finally:
        con.close()

    con_rw = duckdb.connect(db_path)
    try:
        log_event(
            con_rw, "export_table", "export_parquet", table_name,
            f"output_path={output_path}, rows={row_count}",
        )
    finally:
        con_rw.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export a DuckDB table to a Parquet file.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--table_name", required=True,
        help="Name of the table to export.",
    )
    parser.add_argument(
        "--output_path", required=True,
        help="Destination file path for the Parquet output.",
    )
    args = parser.parse_args()
    export_table(args.db, args.table_name, args.output_path)


if __name__ == "__main__":
    main()
