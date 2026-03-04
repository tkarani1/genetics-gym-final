#!/usr/bin/env python3
"""
CLI tool for creating consolidated VSM (Variant Scoring Method) tables.

Reads prediction and evaluation tables from GCS or local paths, merges them
on genomic keys (chrom, pos, ref, alt), optionally computes percentile ranks
for prediction score columns, and writes a single output parquet file.
"""
import argparse
import os
import sys
import tempfile
import time

import polars as pl

from table_io import ensure_parquet, write_parquet
from merge import merge_tables, JOIN_KEYS
from percentile import add_percentile_columns
from apply_filters import apply_filters


def _parse_uri_list(raw: str) -> list[str]:
    return [s.strip() for s in raw.split(",") if s.strip()]


def _score_columns(lf: pl.LazyFrame) -> list[str]:
    """Return all columns in *lf* that are not join keys."""
    return [c for c in lf.collect_schema().names() if c not in JOIN_KEYS]


def run_pipeline(
    prediction_uris: list[str],
    evaluation_uris: list[str],
    output_uri: str,
    join_type: str = "inner",
    percentile_order: str = "post",
    filter_uris: list[str] | None = None,
    keep_raw_scores: bool = False,
) -> None:
    """
    Core pipeline -- decoupled from CLI for reuse (e.g. future resources.json).
    """
    start = time.perf_counter()

    cache_dir = os.path.join(tempfile.gettempdir(), "vsm_table_cache")

    # --- Load prediction tables -------------------------------------------
    pred_frames: list[pl.LazyFrame] = []
    all_score_cols: list[str] = []

    for uri in prediction_uris:
        pq_path = ensure_parquet(uri, cache_dir)
        lf = pl.scan_parquet(
            pq_path,
            **({"storage_options": {"token": "google_default"}}
               if pq_path.startswith("gs://") else {}),
        )
        score_cols = _score_columns(lf)
        all_score_cols.extend(score_cols)
        print(
            f"  Prediction table {uri}: score columns = {score_cols}",
            file=sys.stderr,
        )
        pred_frames.append(lf)

    # --- Load evaluation tables -------------------------------------------
    eval_frames: list[pl.LazyFrame] = []

    for uri in evaluation_uris:
        pq_path = ensure_parquet(uri, cache_dir)
        lf = pl.scan_parquet(
            pq_path,
            **({"storage_options": {"token": "google_default"}}
               if pq_path.startswith("gs://") else {}),
        )
        print(
            f"  Evaluation table {uri}: columns = {_score_columns(lf)}",
            file=sys.stderr,
        )
        eval_frames.append(lf)

    # --- Pre-merge percentile (if requested) ------------------------------
    if percentile_order == "pre" and join_type == "inner":
        print("  Calculating percentiles BEFORE merge ...", file=sys.stderr)
        pred_frames = [
            add_percentile_columns(lf, _score_columns(lf))
            for lf in pred_frames
        ]

    # --- Merge all tables -------------------------------------------------
    all_frames = pred_frames + eval_frames
    print(
        f"  Merging {len(all_frames)} tables ({join_type} join) ...",
        file=sys.stderr,
    )
    merged = merge_tables(all_frames, join_type=join_type)

    # --- Post-merge percentile --------------------------------------------
    if percentile_order == "post" or join_type == "outer":
        print("  Calculating percentiles AFTER merge ...", file=sys.stderr)
        merged = add_percentile_columns(merged, all_score_cols)

    # --- Drop raw score columns (if requested) -------------------------------
    if not keep_raw_scores and all_score_cols:
        print("  Dropping raw score columns ...", file=sys.stderr)
        merged = merged.drop(all_score_cols)

    # --- Apply filter columns (if requested) --------------------------------
    if filter_uris:
        print(
            f"  Applying {len(filter_uris)} filter table(s) ...",
            file=sys.stderr,
        )
        merged = apply_filters(merged, filter_uris, cache_dir)

    # --- Write output -----------------------------------------------------
    print(f"  Writing output to {output_uri} ...", file=sys.stderr)
    write_parquet(merged, output_uri)

    elapsed = time.perf_counter() - start
    print(f"Done in {elapsed:.1f}s.", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a consolidated VSM table from prediction and evaluation data.",
    )

    parser.add_argument(
        "--prediction_tables",
        required=True,
        help="Comma-separated URIs (GCS or local) of prediction parquet/tsv files.",
    )
    parser.add_argument(
        "--evaluation_tables",
        required=True,
        help="Comma-separated URIs (GCS or local) of evaluation parquet/tsv files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination URI for the merged parquet file (GCS or local).",
    )
    parser.add_argument(
        "--join_type",
        choices=["inner", "outer"],
        default="inner",
        help="Join strategy across all tables (default: inner).",
    )
    parser.add_argument(
        "--percentile_order",
        choices=["pre", "post"],
        default="post",
        help=(
            "'pre' = compute percentile ranks on each prediction table before "
            "merging; 'post' = compute after merging (default: post). "
            "For outer joins the order is irrelevant; post is always used."
        ),
    )
    parser.add_argument(
        "--keep_raw_scores",
        action="store_true",
        default=False,
        help=(
            "Keep the original score columns alongside their percentile "
            "counterparts. By default, raw score columns are dropped and "
            "only the percentile columns are retained."
        ),
    )
    parser.add_argument(
        "--filter_tables",
        default=None,
        help=(
            "Comma-separated URIs of filter tables. Each adds a boolean "
            "column indicating whether the variant key is present in that "
            "filter table."
        ),
    )

    args = parser.parse_args()

    run_pipeline(
        prediction_uris=_parse_uri_list(args.prediction_tables),
        evaluation_uris=_parse_uri_list(args.evaluation_tables),
        output_uri=args.output,
        join_type=args.join_type,
        percentile_order=args.percentile_order,
        filter_uris=_parse_uri_list(args.filter_tables) if args.filter_tables else None,
        keep_raw_scores=args.keep_raw_scores,
    )


if __name__ == "__main__":
    main()
