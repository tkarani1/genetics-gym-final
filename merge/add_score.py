#!/usr/bin/env python3
"""
Join a new score table onto a precomputed merged prediction parquet,
avoiding a full pipeline re-run.

Guarantees: merge(A, B, C) == join(merge(A, B), C) for inner/outer
joins on unique-keyed tables.

For pairwise tables, use --pairwise_anchor to also create the
pairwise anchor-masked columns for the new score(s), matching the
behaviour of merge_tables_pairwise in the full pipeline.
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
        "--pairwise_anchor",
        default=None,
        help=(
            "Name of the pairwise anchor score (e.g. 'polyphen_score'). "
            "When set, the base table must contain a '_raw_anchor_{name}' "
            "column (written by create_vsm_table.py for pairwise subtables). "
            "Pairwise anchor-masked columns are created for each new score, "
            "replicating the merge_tables_pairwise behaviour."
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

    if args.pairwise_anchor:
        raw_anchor_col = f"_raw_anchor_{args.pairwise_anchor}"
        base_schema = base_lf.collect_schema()
        if raw_anchor_col not in base_schema.names():
            raise ValueError(
                f"Base table does not contain '{raw_anchor_col}'. "
                f"Pairwise add-one requires the base table to have been "
                f"written with raw anchor retention (create_vsm_table.py "
                f"pairwise subtable mode)."
            )

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

    # --- Pairwise column creation -----------------------------------------
    pairwise_cols: list[str] = []
    if args.pairwise_anchor:
        anchor = args.pairwise_anchor
        raw_anchor_col = f"_raw_anchor_{anchor}"
        print(
            f"  Creating pairwise columns using {raw_anchor_col} ...",
            file=sys.stderr,
        )
        pairwise_exprs: list[pl.Expr] = []
        for c in score_cols:
            both_non_null = pl.col(c).is_not_null() & pl.col(raw_anchor_col).is_not_null()
            pw_name = f"{anchor}_pairwise_{c}"
            pairwise_exprs.append(
                pl.when(both_non_null).then(pl.col(raw_anchor_col)).alias(pw_name)
            )
            pairwise_exprs.append(
                pl.when(both_non_null).then(pl.col(c)).alias(c)
            )
            pairwise_cols.append(pw_name)
        result = result.with_columns(pairwise_exprs)

    # --- Percentile computation (post) ------------------------------------
    if args.percentile_order == "post":
        print("  Computing percentiles AFTER join ...", file=sys.stderr)
        all_new_cols = score_cols + pairwise_cols
        result = add_percentile_columns(result, all_new_cols)
        if pairwise_cols:
            anchor = args.pairwise_anchor
            pct_renames: dict[str, str] = {}
            for c in score_cols:
                pct_renames[f"{c}_percentile"] = f"{c}_percentile_with_anchor"
                pct_renames[f"{anchor}_pairwise_{c}_percentile"] = (
                    f"{anchor}_anchor_percentile_with_{c}"
                )
            result = result.rename(pct_renames)

    print(f"  Writing output to {args.output} ...", file=sys.stderr)
    result.sink_parquet(args.output, compression="zstd")

    elapsed = time.perf_counter() - start
    print(f"Done in {elapsed:.1f}s.", file=sys.stderr)


if __name__ == "__main__":
    main()
