#!/usr/bin/env python3
"""Remove a table (score or eval) and its metadata row from a DuckDB database."""
from __future__ import annotations

import argparse
import sys

import duckdb

from initialize_db import log_event


DEFAULT_DB_PATH = "scores.duckdb"


def eject_table(db_path: str, table_name: str) -> None:
    """Drop the table *table_name* and delete its metadata row.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    table_name : str
        Name of the table to remove.
    """
    con = duckdb.connect(db_path)
    try:
        exists = con.execute(
            "SELECT COUNT(*) FROM information_schema.tables "
            "WHERE table_name = ? AND table_schema = 'main'",
            [table_name],
        ).fetchone()[0]
        if not exists:
            raise ValueError(f"Table {table_name!r} does not exist.")

        con.execute(f'DROP TABLE "{table_name}"')
        con.execute(
            "DELETE FROM metadata WHERE table_name = ?",
            [table_name],
        )
        log_event(con, "eject_table", "drop_table", table_name)
        print(f"Ejected table {table_name!r}", file=sys.stderr)
    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove a table (score or eval) and its metadata from a DuckDB database.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--table_name", required=True,
        help="Name of the table to remove.",
    )
    args = parser.parse_args()
    eject_table(args.db, args.table_name)


if __name__ == "__main__":
    main()
