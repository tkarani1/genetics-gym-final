#!/usr/bin/env python3
"""Ingest a single score column into a wide DuckDB score table."""
from __future__ import annotations

import argparse
import os
import sys

import duckdb

from initialize_db import DEFAULT_DB_PATH, log_event

WIDE_TABLE_DEFAULTS = {"variant": "variant_scores", "gene": "gene_scores"}
KEY_NAMES = {"chrom", "pos", "ref", "alt", "ensg", "key"}


def ingest_score(
    db_path: str,
    score_name: str,
    score_path: str,
    table_name: str,
    analysis_level: str,
) -> None:
    """Add a score column to a wide score table from a Parquet source.

    If the wide table does not exist yet, it is created with canonical key
    columns plus the first score column.  If it already exists, the new
    score column is added and populated.  Rows in the source that are new
    to the table are inserted; existing rows are updated.

    Incoming data is deduplicated on the key before insertion (keeps the
    row where the score is non-null, breaking ties randomly).

    Parameters
    ----------
    db_path : str
        Path to the persistent ``.duckdb`` file (must already be initialised).
    score_name : str
        Column name in the parquet file to use as the score value.
    score_path : str
        Path to the source parquet file.
    table_name : str
        Name for the wide score table (e.g. ``variant_scores``).
    analysis_level : str
        ``"variant"`` or ``"gene"``.
    """
    score_path = os.path.abspath(score_path)

    if not os.path.isfile(db_path):
        raise FileNotFoundError(f"Database not found: {db_path}")
    if not os.path.isfile(score_path):
        raise FileNotFoundError(f"Parquet file not found: {score_path}")
    if analysis_level not in ("variant", "gene"):
        raise ValueError(f"analysis_level must be 'variant' or 'gene', got {analysis_level!r}")

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

        parquet_columns = {
            r[0]
            for r in con.execute(
                f"SELECT name FROM parquet_schema('{score_path}')"
            ).fetchall()
        }
        if score_name not in parquet_columns:
            raise ValueError(
                f"Column {score_name!r} not found in parquet. "
                f"Available: {sorted(parquet_columns)}"
            )

        if analysis_level == "variant":
            missing = [k for k in ("chrom", "pos", "ref", "alt") if k not in parquet_columns]
            if missing:
                raise ValueError(
                    f"Variant key column(s) missing from parquet: {missing}"
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
            if "ensg" not in parquet_columns:
                raise ValueError("Key column 'ensg' not found in parquet.")
            key_select = (
                "NULL::VARCHAR AS chrom, "
                "NULL::BIGINT AS pos, "
                "NULL::VARCHAR AS ref, "
                "NULL::VARCHAR AS alt, "
                "ensg"
            )
            hash_expr = "hash(ensg)"

        quoted_table = f'"{table_name}"'
        quoted_score = f'"{score_name}"'

        deduped_src = f"""
            (SELECT * EXCLUDE (rn) FROM (
                SELECT
                    {key_select},
                    {hash_expr} AS "key",
                    {quoted_score}::DOUBLE AS score,
                    ROW_NUMBER() OVER (
                        PARTITION BY {hash_expr}
                        ORDER BY ({quoted_score} IS NOT NULL) DESC, random()
                    ) AS rn
                FROM read_parquet('{score_path}')
            ) WHERE rn = 1)
        """

        table_exists = table_name in tables

        if not table_exists:
            con.execute(
                f"CREATE TABLE {quoted_table} AS\n"
                f"SELECT chrom, pos, ref, alt, ensg, \"key\", "
                f"score AS {quoted_score}\n"
                f"FROM {deduped_src}"
            )
            row_count = con.execute(
                f"SELECT COUNT(*) FROM {quoted_table}"
            ).fetchone()[0]
            scored_count = con.execute(
                f"SELECT COUNT({quoted_score}) FROM {quoted_table}"
            ).fetchone()[0]
            print(
                f"Created {table_name!r} with {score_name!r}: "
                f"{row_count:,} rows ({scored_count:,} scored, "
                f"{row_count - scored_count:,} null)",
                file=sys.stderr,
            )
        else:
            existing_cols = {
                desc[0]
                for desc in con.execute(
                    f"SELECT * FROM {quoted_table} LIMIT 0"
                ).description
            }
            if score_name in existing_cols:
                raise ValueError(
                    f"Score column {score_name!r} already exists in {table_name!r}."
                )

            con.execute(
                f"ALTER TABLE {quoted_table} ADD COLUMN {quoted_score} DOUBLE"
            )

            updated = con.execute(f"""
                UPDATE {quoted_table} AS t
                SET {quoted_score} = src.score
                FROM {deduped_src} AS src
                WHERE t."key" = src."key"
            """).fetchone()

            new_rows = con.execute(f"""
                INSERT INTO {quoted_table}
                SELECT src.chrom, src.pos, src.ref, src.alt, src.ensg,
                       src."key", {_null_fills(con, table_name, score_name)},
                       src.score AS {quoted_score}
                FROM {deduped_src} AS src
                WHERE src."key" NOT IN (SELECT "key" FROM {quoted_table})
            """).fetchone()

            row_count = con.execute(
                f"SELECT COUNT(*) FROM {quoted_table}"
            ).fetchone()[0]
            scored_count = con.execute(
                f"SELECT COUNT({quoted_score}) FROM {quoted_table}"
            ).fetchone()[0]
            new_row_count = row_count - (
                con.execute(
                    f'SELECT COUNT(*) FROM {quoted_table} WHERE {quoted_score} IS NULL'
                ).fetchone()[0]
                + scored_count
                - row_count
            ) if False else 0  # placeholder; we report from scored_count

            print(
                f"Added {score_name!r} to {table_name!r}: "
                f"{row_count:,} total rows ({scored_count:,} scored, "
                f"{row_count - scored_count:,} null)",
                file=sys.stderr,
            )

        con.execute(
            "INSERT INTO metadata "
            "(source_column, source_path, table_name, table_type, "
            "analysis_level, deduped) "
            "VALUES (?, ?, ?, 'score', ?, TRUE)",
            [score_name, score_path, table_name, analysis_level],
        )
        log_event(
            con, "ingest_score", "add_column", table_name,
            f"source_column={score_name}, source_path={score_path}, "
            f"analysis_level={analysis_level}, rows={row_count}",
        )

    finally:
        con.close()


def _null_fills(
    con: duckdb.DuckDBPyConnection,
    table_name: str,
    new_score: str,
) -> str:
    """Build NULL fill expressions for all existing score columns except keys
    and the new score being added, for the INSERT of new rows."""
    cols = [
        desc[0]
        for desc in con.execute(
            f'SELECT * FROM "{table_name}" LIMIT 0'
        ).description
    ]
    fills = []
    for c in cols:
        if c in KEY_NAMES or c == new_score:
            continue
        fills.append(f'NULL AS "{c}"')
    return ", ".join(fills)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ingest a score column into a wide DuckDB score table.",
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB_PATH,
        help=f"Path to the .duckdb file (default: {DEFAULT_DB_PATH}).",
    )
    parser.add_argument(
        "--score_name", required=True,
        help="Name of the score column in the parquet file.",
    )
    parser.add_argument(
        "--score_path", required=True,
        help="Path to the source parquet file.",
    )
    parser.add_argument(
        "--table_name", default=None,
        help="Name for the wide score table (default: variant_scores or gene_scores).",
    )
    parser.add_argument(
        "--analysis_level", required=True, choices=["variant", "gene"],
        help="Key type: 'variant' (chrom/pos/ref/alt) or 'gene' (ensg).",
    )
    args = parser.parse_args()
    table_name = args.table_name or WIDE_TABLE_DEFAULTS[args.analysis_level]
    ingest_score(args.db, args.score_name, args.score_path,
                 table_name, args.analysis_level)


if __name__ == "__main__":
    main()
