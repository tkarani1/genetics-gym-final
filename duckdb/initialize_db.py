#!/usr/bin/env python3
"""Create a persistent DuckDB database with an empty metadata table."""
from __future__ import annotations

import argparse
import os
import sys

import duckdb


DEFAULT_DB_PATH = "gg_data.duckdb"

METADATA_DDL = """\
CREATE TABLE metadata (
    source_column  VARCHAR NOT NULL,
    source_path    VARCHAR NOT NULL,
    table_name     VARCHAR NOT NULL PRIMARY KEY,
    table_type      VARCHAR NOT NULL CHECK (table_type IN ('score', 'eval', 'merged_scores', 'merged_evals')),
    analysis_level  VARCHAR NOT NULL CHECK (analysis_level IN ('variant', 'gene')),
    deduped        BOOLEAN NOT NULL DEFAULT FALSE
);
"""


def initialize_db(db_path: str) -> None:
    """Create a new .duckdb file containing only the *metadata* table.

    Raises ``FileExistsError`` if *db_path* already exists to prevent
    accidentally overwriting a populated database.
    """
    if os.path.exists(db_path):
        raise FileExistsError(f"Database already exists: {db_path}")

    con = duckdb.connect(db_path)
    try:
        con.execute(METADATA_DDL)
        print(f"Initialized database: {db_path}", file=sys.stderr)
    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a persistent DuckDB database with a metadata table.",
    )
    parser.add_argument(
        "--db",
        default=DEFAULT_DB_PATH,
        help=f"Path for the new .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    args = parser.parse_args()
    initialize_db(args.db)


if __name__ == "__main__":
    main()
