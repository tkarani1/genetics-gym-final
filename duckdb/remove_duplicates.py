#!/usr/bin/env python3
"""Remove duplicate-key rows from a DuckDB score or eval table."""
from __future__ import annotations

import argparse
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

STRATEGIES = ("keep_random", "prefer_scored")


def remove_duplicates(
    db_path: str,
    table_name: str,
    strategy: str = "keep_random",
) -> None:
    """De-duplicate rows in *table_name* based on the ``key`` column.

    After deduplication the ``deduped`` flag in the metadata table is
    set to ``TRUE``.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    table_name : str
        Name of the table to deduplicate.
    strategy : str
        ``"keep_random"`` — for each set of duplicate keys, one row is
        kept at random and the rest are deleted.
        ``"prefer_scored"`` — (score tables only) keeps the row with a
        non-null score when possible, breaking ties at random.
    """
    if strategy not in STRATEGIES:
        raise ValueError(
            f"Unknown strategy {strategy!r}. Choose from: {STRATEGIES}"
        )

    con = duckdb.connect(db_path)
    try:
        exists = con.execute(
            "SELECT COUNT(*) FROM information_schema.tables "
            "WHERE table_name = ? AND table_schema = 'main'",
            [table_name],
        ).fetchone()[0]
        if not exists:
            raise ValueError(f"Table {table_name!r} does not exist.")

        quoted = f'"{table_name}"'

        before = con.execute(
            f"SELECT COUNT(*) FROM {quoted}"
        ).fetchone()[0]
        unique_keys = con.execute(
            f'SELECT COUNT(DISTINCT "key") FROM {quoted}'
        ).fetchone()[0]
        dup_count = before - unique_keys

        if dup_count == 0:
            print(
                f"No duplicates in {table_name!r} ({before} rows, "
                f"{unique_keys} unique keys). Marking as deduped.",
                file=sys.stderr,
            )
            con.execute(
                "UPDATE metadata SET deduped = TRUE WHERE table_name = ?",
                [table_name],
            )
            return

        table_type = con.execute(
            "SELECT table_type FROM metadata WHERE table_name = ?",
            [table_name],
        ).fetchone()
        is_score = table_type and table_type[0] == "score"

        if strategy == "prefer_scored" and not is_score:
            raise ValueError(
                f"Strategy 'prefer_scored' is only valid for score tables, "
                f"but {table_name!r} is an eval table."
            )

        if strategy == "prefer_scored":
            window_order = "(score IS NOT NULL) DESC, random()"
        else:
            window_order = "random()"

        order_clause = "ORDER BY score NULLS LAST" if is_score else ""
        con.execute(f"""
            CREATE OR REPLACE TABLE {quoted} AS
            SELECT * EXCLUDE (rn) FROM (
                SELECT *, ROW_NUMBER() OVER (
                    PARTITION BY "key" ORDER BY {window_order}
                ) AS rn
                FROM {quoted}
            )
            WHERE rn = 1
            {order_clause};
        """)

        con.execute(
            "UPDATE metadata SET deduped = TRUE WHERE table_name = ?",
            [table_name],
        )

        after = con.execute(
            f"SELECT COUNT(*) FROM {quoted}"
        ).fetchone()[0]
        log_event(
            con, "remove_duplicates", "deduplicate", table_name,
            f"strategy={strategy}, rows_before={before}, rows_after={after}",
        )
        print(
            f"Deduped {table_name!r}: {before} → {after} rows "
            f"({before - after} duplicates removed)",
            file=sys.stderr,
        )
    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove duplicate-key rows from a DuckDB table.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--table_name", required=True,
        help="Name of the table to deduplicate.",
    )
    parser.add_argument(
        "--strategy", default="keep_random", choices=list(STRATEGIES),
        help="Deduplication strategy (default: keep_random).",
    )
    args = parser.parse_args()
    remove_duplicates(args.db, args.table_name, args.strategy)


if __name__ == "__main__":
    main()
