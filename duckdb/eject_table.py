#!/usr/bin/env python3
"""Remove a table (or a score column from a wide table) and its metadata."""
from __future__ import annotations

import argparse
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

KEY_NAMES = {"chrom", "pos", "ref", "alt", "ensg", "key"}


def eject_table(db_path: str, table_name: str, *,
                score_column: str | None = None) -> None:
    """Drop a table or a score column from a wide score table.

    For eval / merged tables (or when *score_column* is not given), the
    entire table is dropped along with all of its metadata rows.

    For score tables, if *score_column* is provided, only that column is
    removed from the wide table and the corresponding metadata row is
    deleted.  If no score columns remain after the drop, the whole table
    is removed.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    table_name : str
        Name of the table to modify / remove.
    score_column : str or None
        If specified, drop only this score column from a wide table.
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

        if score_column is not None:
            existing_cols = {
                desc[0]
                for desc in con.execute(
                    f'SELECT * FROM "{table_name}" LIMIT 0'
                ).description
            }
            if score_column not in existing_cols:
                raise ValueError(
                    f"Column {score_column!r} not found in {table_name!r}."
                )

            remaining_score_cols = existing_cols - KEY_NAMES - {score_column}
            if remaining_score_cols:
                con.execute(
                    f'ALTER TABLE "{table_name}" DROP COLUMN "{score_column}"'
                )
                con.execute(
                    "DELETE FROM metadata "
                    "WHERE table_name = ? AND source_column = ?",
                    [table_name, score_column],
                )
                log_event(con, "eject_table", "drop_column", table_name,
                          f"column={score_column}")
                print(
                    f"Dropped column {score_column!r} from {table_name!r}",
                    file=sys.stderr,
                )
            else:
                con.execute(f'DROP TABLE "{table_name}"')
                con.execute(
                    "DELETE FROM metadata WHERE table_name = ?",
                    [table_name],
                )
                log_event(con, "eject_table", "drop_table", table_name,
                          f"last_column={score_column}")
                print(
                    f"Dropped last column {score_column!r} → "
                    f"removed table {table_name!r}",
                    file=sys.stderr,
                )
        else:
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
        description="Remove a table or a score column from a DuckDB database.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--table_name", required=True,
        help="Name of the table to modify / remove.",
    )
    parser.add_argument(
        "--score_column", default=None,
        help="If given, drop only this column from a wide score table.",
    )
    args = parser.parse_args()
    eject_table(args.db, args.table_name, score_column=args.score_column)


if __name__ == "__main__":
    main()
