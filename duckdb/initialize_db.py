#!/usr/bin/env python3
"""Create a persistent DuckDB database with metadata and audit_log tables."""
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
    table_name     VARCHAR NOT NULL,
    table_type      VARCHAR NOT NULL CHECK (table_type IN ('score', 'eval', 'merged_scores', 'merged_evals', 'merged_analysis')),
    analysis_level  VARCHAR NOT NULL CHECK (analysis_level IN ('variant', 'gene')),
    deduped        BOOLEAN NOT NULL DEFAULT FALSE,
    eval_column    VARCHAR,
    case_column    VARCHAR,
    ctrl_column    VARCHAR,
    UNIQUE (table_name, source_column)
);
"""

AUDIT_LOG_DDL = """\
CREATE TABLE audit_log (
    ts          TIMESTAMP NOT NULL DEFAULT current_timestamp,
    module      VARCHAR NOT NULL,
    action      VARCHAR NOT NULL,
    table_name  VARCHAR,
    details     VARCHAR
);
"""


def log_event(
    con: duckdb.DuckDBPyConnection,
    module: str,
    action: str,
    table_name: str | None = None,
    details: str | None = None,
) -> None:
    """Append a row to the audit_log table.

    Safe to call even if the database pre-dates the audit_log table;
    the INSERT is silently skipped in that case.
    """
    try:
        con.execute(
            "INSERT INTO audit_log (module, action, table_name, details) "
            "VALUES (?, ?, ?, ?)",
            [module, action, table_name, details],
        )
    except duckdb.CatalogException:
        pass


def initialize_db(db_path: str) -> None:
    """Create a new .duckdb file containing *metadata* and *audit_log* tables.

    Raises ``FileExistsError`` if *db_path* already exists to prevent
    accidentally overwriting a populated database.
    """
    if os.path.exists(db_path):
        raise FileExistsError(f"Database already exists: {db_path}")

    con = duckdb.connect(db_path)
    try:
        con.execute(METADATA_DDL)
        con.execute(AUDIT_LOG_DDL)
        log_event(con, "initialize_db", "create_database",
                  details=f"db_path={db_path}")
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
