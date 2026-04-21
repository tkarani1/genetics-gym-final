#!/usr/bin/env python3
"""Compute pairwise percentiles between an anchor and a non-anchor score table."""
from __future__ import annotations

import argparse
import sys

import duckdb


DEFAULT_DB_PATH = "scores.duckdb"


def pairwise_percentile(
    db_path: str,
    anchor_table: str,
    table_name: str,
) -> None:
    """Compute intersection-based pairwise percentiles and write them into
    the non-anchor table's ``temp_1`` and ``temp_2`` columns.

    ``temp_1`` receives the percentile of the *anchor* table's scores
    among the intersection, and ``temp_2`` receives the percentile of the
    *non-anchor* table's scores among the intersection.  Rows not in the
    intersection retain ``NULL`` in both columns.  Null scores within the
    intersection also receive ``NULL`` percentiles and are excluded from
    the ranking denominator.

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file.
    anchor_table : str
        Name of the anchor score table.
    table_name : str
        Name of the non-anchor score table whose temp columns are written.
    """
    con = duckdb.connect(db_path)
    try:
        for tbl in (anchor_table, table_name):
            meta = con.execute(
                "SELECT table_type, deduped, analysis_level "
                "FROM metadata WHERE table_name = ?",
                [tbl],
            ).fetchone()
            if meta is None:
                raise ValueError(f"Table {tbl!r} not found in metadata.")
            ttype, deduped, analysis_level = meta
            if ttype != "score":
                raise ValueError(
                    f"Table {tbl!r} is type {ttype!r}, not 'score'."
                )
            if not deduped:
                raise ValueError(
                    f"Table {tbl!r} has not been deduped. "
                    f"Run remove_duplicates first."
                )
            if analysis_level != "variant":
                raise ValueError(
                    f"Table {tbl!r} has analysis_level {analysis_level!r}. "
                    f"Pairwise operations require all tables to be "
                    f"analysis_level 'variant'."
                )

        quoted_anchor = f'"{anchor_table}"'
        quoted_target = f'"{table_name}"'

        # Reset temp columns before computing
        con.execute(
            f"UPDATE {quoted_target} SET temp_1 = NULL, temp_2 = NULL"
        )

        con.execute(f"""
            WITH intersection AS (
                SELECT
                    a."key",
                    a.score AS anchor_score,
                    b.score AS nonanchor_score
                FROM {quoted_anchor} a
                INNER JOIN {quoted_target} b ON a."key" = b."key"
            ),
            ranked AS (
                SELECT
                    "key",
                    CASE WHEN anchor_score IS NOT NULL
                         THEN PERCENT_RANK() OVER (
                             PARTITION BY (anchor_score IS NOT NULL)
                             ORDER BY anchor_score
                         )
                    END AS temp_1,
                    CASE WHEN nonanchor_score IS NOT NULL
                         THEN PERCENT_RANK() OVER (
                             PARTITION BY (nonanchor_score IS NOT NULL)
                             ORDER BY nonanchor_score
                         )
                    END AS temp_2
                FROM intersection
            )
            UPDATE {quoted_target} t
            SET temp_1 = r.temp_1, temp_2 = r.temp_2
            FROM ranked r
            WHERE t."key" = r."key";
        """)

        stats = con.execute(f"""
            SELECT
                COUNT(*) FILTER (WHERE temp_1 IS NOT NULL OR temp_2 IS NOT NULL),
                COUNT(temp_1),
                COUNT(temp_2)
            FROM {quoted_target}
        """).fetchone()
        intersection_size, anchor_scored, target_scored = stats

        print(
            f"Pairwise percentiles {anchor_table!r} x {table_name!r}: "
            f"{intersection_size} intersection rows "
            f"({anchor_scored} anchor scored, {target_scored} non-anchor scored)",
            file=sys.stderr,
        )
    finally:
        con.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute pairwise percentiles between an anchor and a "
            "non-anchor score table."
        ),
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--anchor_table", required=True,
        help="Name of the anchor score table.",
    )
    parser.add_argument(
        "--table_name", required=True,
        help="Name of the non-anchor score table to write percentiles into.",
    )
    args = parser.parse_args()
    pairwise_percentile(args.db, args.anchor_table, args.table_name)


if __name__ == "__main__":
    main()
