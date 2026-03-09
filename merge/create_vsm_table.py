#!/usr/bin/env python3
"""
CLI tool for creating consolidated VSM (Variant Scoring Method) tables.

Reads prediction and evaluation tables from GCS or local paths, merges them
on genomic keys (chrom, pos, ref, alt), optionally computes percentile ranks
for prediction score columns, and writes a single output parquet file.
"""
from __future__ import annotations

import argparse
import os
import posixpath
import sys
import tempfile
import time
from collections import Counter

import polars as pl

from table_io import ensure_parquet, write_parquet
from merge import merge_tables, JOIN_KEYS
from percentile import add_percentile_columns
from apply_filters import apply_filters


def _parse_uri_list(raw: str) -> list[str]:
    return [s.strip() for s in raw.split(",") if s.strip()]


_KNOWN_EXTENSIONS = (".tsv.bgz", ".tsv.gz", ".tsv", ".parquet")


def _derive_stem(uri: str) -> str:
    """Filename without directory or known extensions."""
    basename = posixpath.basename(uri.rstrip("/"))
    if not basename:
        basename = os.path.basename(uri.rstrip(os.sep))
    lower = basename.lower()
    for ext in _KNOWN_EXTENSIONS:
        if lower.endswith(ext):
            return basename[: len(basename) - len(ext)]
    return os.path.splitext(basename)[0]


def _score_columns(lf: pl.LazyFrame) -> list[str]:
    """Return Float64 non-key columns in *lf* (the actual numeric scores)."""
    schema = lf.collect_schema()
    return [c for c, dtype in schema.items()
            if c not in JOIN_KEYS and dtype == pl.Float64]


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
        lf = lf.select(JOIN_KEYS + score_cols)
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
        schema = lf.collect_schema()
        label_cols = [c for c, dtype in schema.items()
                      if c not in JOIN_KEYS and dtype == pl.Boolean]
        lf = lf.select(JOIN_KEYS + label_cols)
        print(
            f"  Evaluation table {uri}: label columns = {label_cols}",
            file=sys.stderr,
        )
        eval_frames.append(lf)

    # --- Auto-suffix colliding eval label columns -------------------------
    all_label_names = [
        c
        for lf in eval_frames
        for c in lf.collect_schema().names()
        if c not in JOIN_KEYS
    ]
    collisions = {
        name for name, count in Counter(all_label_names).items() if count > 1
    }
    if collisions:
        print(
            f"  Renaming colliding label columns: {collisions}",
            file=sys.stderr,
        )
        for i, uri in enumerate(evaluation_uris):
            stem = _derive_stem(uri)
            frame_cols = set(eval_frames[i].collect_schema().names())
            renames = {
                col: f"{col}__{stem}"
                for col in collisions
                if col in frame_cols
            }
            if renames:
                eval_frames[i] = eval_frames[i].rename(renames)

    # --- Phase 1: Outer-join eval tables (always) --------------------------
    print(
        f"  Merging {len(eval_frames)} eval table(s) (outer join) ...",
        file=sys.stderr,
    )
    merged_eval = merge_tables(eval_frames, join_type="outer")

    # --- Phase 2: Join pred tables with percentiles -----------------------
    if percentile_order == "pre":
        print("  Calculating percentiles BEFORE pred merge ...", file=sys.stderr)
        pred_frames = [
            add_percentile_columns(lf, _score_columns(lf))
            for lf in pred_frames
        ]

    print(
        f"  Merging {len(pred_frames)} pred table(s) ({join_type} join) ...",
        file=sys.stderr,
    )
    merged_pred = merge_tables(pred_frames, join_type=join_type)

    if percentile_order == "post":
        print("  Calculating percentiles AFTER pred merge ...", file=sys.stderr)
        merged_pred = add_percentile_columns(merged_pred, all_score_cols)

    # --- Phase 3: Left-join eval onto pred --------------------------------
    print("  Left-joining eval and pred ...", file=sys.stderr)
    merged = merged_eval.join(
        merged_pred, on=JOIN_KEYS, how="left", coalesce=True,
    )

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
    # write_parquet(merged, output_uri)
    merged.sink_parquet(
        output_uri,
        compression="zstd"
    )

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
        help=(
            "Join strategy for prediction tables (default: inner). "
            "Evaluation tables are always outer-joined."
        ),
    )
    parser.add_argument(
        "--percentile_order",
        choices=["pre", "post"],
        default="post",
        help=(
            "'pre' = compute percentile ranks on each prediction table before "
            "merging; 'post' = compute after merging prediction tables "
            "(default: post)."
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
