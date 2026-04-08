#!/usr/bin/env python3
"""
Join a new score table onto a precomputed merged prediction parquet,
avoiding a full pipeline re-run.

Guarantees: merge(A, B, C) == join(merge(A, B), C) for inner/outer
joins on unique-keyed tables.
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time

import polars as pl

from table_io import ensure_parquet, normalize_chrom_key
from merge import JOIN_KEYS
from percentile import add_percentile_columns


def _score_columns(lf: pl.LazyFrame) -> list[str]:
    schema = lf.collect_schema()
    return [c for c, dtype in schema.items()
            if c not in JOIN_KEYS and dtype == pl.Float64]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Add a new score column to a precomputed merged prediction "
            "table via a single join, without re-running the full pipeline."
        ),
    )
    parser.add_argument(
        "--base",
        required=True,
        help="Path to the precomputed merged prediction parquet.",
    )
    parser.add_argument(
        "--new_table",
        required=True,
        help=(
            "Path to the new score table to join in "
            "(parquet, .tsv, .tsv.bgz, .tsv.gz)."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination path for the merged output parquet.",
    )
    parser.add_argument(
        "--join_type",
        choices=["inner", "outer"],
        required=True,
        help="Join strategy: 'inner' or 'outer'.",
    )
    parser.add_argument(
        "--score_columns",
        default=None,
        help=(
            "Comma-separated list of score columns to select from the new "
            "table. If omitted, all Float64 non-key columns are included."
        ),
    )
    parser.add_argument(
        "--percentile_order",
        choices=["pre", "post", "none"],
        default="none",
        help=(
            "Compute percentile ranks for the new score column(s). "
            "'pre' = compute on the new table's own population before "
            "joining (matches --percentile_order pre in the full pipeline). "
            "'post' = compute on the joined population after merging "
            "(matches --percentile_order post). "
            "'none' = skip percentile computation (default)."
        ),
    )
    parser.add_argument(
        "--use_cache",
        action="store_true",
        default=False,
        help="Re-use cached TSV-to-parquet conversions when available.",
    )
    parser.add_argument(
        "--store_cache",
        action="store_true",
        default=False,
        help="Persist TSV-to-parquet conversions for future reuse.",
    )

    args = parser.parse_args()
    start = time.perf_counter()

    cache_dir = os.path.join(tempfile.gettempdir(), "vsm_table_cache")

    print(f"  Loading base table: {args.base}", file=sys.stderr)
    base_lf = normalize_chrom_key(pl.scan_parquet(args.base))

    print(f"  Loading new table: {args.new_table}", file=sys.stderr)
    pq_path = ensure_parquet(
        args.new_table, cache_dir,
        use_cache=args.use_cache, store_cache=args.store_cache,
    )
    new_lf = normalize_chrom_key(pl.scan_parquet(pq_path))

    if args.score_columns:
        score_cols = [c.strip() for c in args.score_columns.split(",") if c.strip()]
    else:
        score_cols = _score_columns(new_lf)

    if not score_cols:
        raise ValueError(
            f"No score columns found in {args.new_table}. "
            "Specify them explicitly with --score_columns."
        )

    new_lf = new_lf.select(JOIN_KEYS + score_cols)
    print(f"  New score columns: {score_cols}", file=sys.stderr)

    if args.percentile_order == "pre":
        print("  Computing percentiles BEFORE join ...", file=sys.stderr)
        new_lf = add_percentile_columns(new_lf, score_cols)

    print(
        f"  Joining ({args.join_type}) onto base "
        f"({len(base_lf.collect_schema().names())} cols) ...",
        file=sys.stderr,
    )
    result = base_lf.join(
        new_lf, on=JOIN_KEYS, how=args.join_type, coalesce=True,
    )

    if args.percentile_order == "post":
        print("  Computing percentiles AFTER join ...", file=sys.stderr)
        result = add_percentile_columns(result, score_cols)

    print(f"  Writing output to {args.output} ...", file=sys.stderr)
    result.sink_parquet(args.output, compression="zstd")

    elapsed = time.perf_counter() - start
    print(f"Done in {elapsed:.1f}s.", file=sys.stderr)


if __name__ == "__main__":
    main()
