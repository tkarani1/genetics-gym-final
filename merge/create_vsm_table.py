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
from functools import reduce

import polars as pl

from table_io import ensure_parquet, normalize_chrom_key, write_parquet
from merge import merge_tables, merge_tables_pairwise, aggregate_by_gene, JOIN_KEYS
from percentile import add_percentile_columns
from smooth import add_smoothed_columns
from apply_filters import apply_filters
from negate import compute_negations, negate_scores


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
    reference_score: str | None = "AM",
    output_table_fields: list[str] | None = None,
    anchor: str | None = None,
    linker_uri: str | None = None,
    smooth: bool = False,
    smooth_reference_dir: str | None = None,
    smooth_sigma: float = 10.0,
) -> None:
    """
    Core pipeline -- decoupled from CLI for reuse (e.g. future resources.json).
    """
    start = time.perf_counter()

    cache_dir = os.path.join(tempfile.gettempdir(), "vsm_table_cache")
    linker_mode = linker_uri is not None

    if linker_mode and percentile_order != "post":
        print(
            "  WARNING: --percentile_order is ignored in linker mode; "
            "percentiles are computed on gene-level aggregates.",
            file=sys.stderr,
        )
    if linker_mode:
        percentile_order = "none"

    # --- Validate output_table_fields against anchor / reference_score ----
    fields_set: set[str] | None = None
    if output_table_fields is not None:
        fields_set = set(output_table_fields)
        if join_type == "pairwise" and anchor not in fields_set:
            raise ValueError(
                f"--output_table_fields must include the anchor column "
                f"'{anchor}' when --join_type is 'pairwise'."
            )
        if reference_score is not None and reference_score not in fields_set:
            raise ValueError(
                f"--output_table_fields must include the reference score "
                f"column '{reference_score}' when score negation is enabled. "
                f"Set --reference_score none to disable negation."
            )

    # --- Load prediction tables -------------------------------------------
    pred_frames: list[pl.LazyFrame] = []
    all_score_cols: list[str] = []

    for uri in prediction_uris:
        pq_path = ensure_parquet(uri, cache_dir)
        lf = normalize_chrom_key(pl.scan_parquet(pq_path))
        score_cols = _score_columns(lf)
        if fields_set is not None:
            score_cols = [c for c in score_cols if c in fields_set]
        if not score_cols:
            print(
                f"  Prediction table {uri}: no matching score columns, skipping.",
                file=sys.stderr,
            )
            continue
        all_score_cols.extend(score_cols)
        lf = lf.select(JOIN_KEYS + score_cols)
        print(
            f"  Prediction table {uri}: score columns = {score_cols}",
            file=sys.stderr,
        )
        pred_frames.append(lf)

    if not pred_frames:
        specified = ", ".join(output_table_fields) if output_table_fields else "(none)"
        raise ValueError(
            f"None of the --output_table_fields ({specified}) matched any "
            f"score columns in the prediction tables."
        )

    # --- Load evaluation tables -------------------------------------------
    eval_keys = ["ensg"] if linker_mode else JOIN_KEYS
    eval_frames: list[pl.LazyFrame] = []

    for uri in evaluation_uris:
        pq_path = ensure_parquet(uri, cache_dir)
        lf = normalize_chrom_key(pl.scan_parquet(pq_path))
        schema = lf.collect_schema()
        if linker_mode:
            label_cols = [c for c in schema.names()
                          if c not in JOIN_KEYS and c != "ensg"]
            lf = lf.select(eval_keys + label_cols)
            print(
                f"  Evaluation table {uri}: label columns = {label_cols}",
                file=sys.stderr,
            )
        else:
            stem = _derive_stem(uri)
            is_pos_cols = [c for c in schema.names() if "is_pos" in c]
            if not is_pos_cols:
                raise ValueError(
                    f"Evaluation table {uri} has no column containing 'is_pos'."
                )
            is_pos_col = is_pos_cols[0]
            target_name = f"is_pos_{stem}"
            lf = lf.select(JOIN_KEYS + [is_pos_col]).rename({is_pos_col: target_name})
            print(
                f"  Evaluation table {uri}: {is_pos_col} -> {target_name}",
                file=sys.stderr,
            )
        eval_frames.append(lf)

    # --- Auto-suffix colliding eval label columns -------------------------
    if linker_mode:
        all_label_names = [
            c
            for lf in eval_frames
            for c in lf.collect_schema().names()
            if c not in eval_keys
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
                    col: f"{col}_{stem}"
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
    if linker_mode:
        def _eval_join(left: pl.LazyFrame, right: pl.LazyFrame) -> pl.LazyFrame:
            return left.join(right, on=eval_keys, how="outer", coalesce=True)

        if len(eval_frames) == 1:
            merged_eval = eval_frames[0]
        else:
            merged_eval = reduce(_eval_join, eval_frames)
    else:
        merged_eval = merge_tables(eval_frames, join_type="outer")

    # --- Phase 2: Join pred tables with percentiles -----------------------
    negate_enabled = reference_score is not None

    if percentile_order == "pre" and (negate_enabled or join_type == "pairwise"):
        reason = (
            "score negation (requires merged frame)"
            if negate_enabled
            else "pairwise join (restructures columns)"
        )
        print(
            f"  WARNING: --percentile_order 'pre' is incompatible with "
            f"{reason}; falling back to 'post'.",
            file=sys.stderr,
        )
        percentile_order = "post"

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
    non_anchor_cols: list[str] = []
    if join_type == "pairwise":
        if anchor not in all_score_cols:
            available = ", ".join(all_score_cols)
            raise ValueError(
                f"Anchor column '{anchor}' not found in prediction tables. "
                f"Available score columns: {available}"
            )
        non_anchor_cols = [c for c in all_score_cols if c != anchor]
        merged_pred, all_score_cols = merge_tables_pairwise(
            pred_frames, all_score_cols, anchor=anchor,
        )
    else:
        merged_pred = merge_tables(pred_frames, join_type=join_type)

    if join_type == "inner" and all_score_cols:
        print("  Dropping rows with any null score column ...", file=sys.stderr)
        merged_pred = merged_pred.drop_nulls(subset=all_score_cols)

    # --- Negate score columns (if reference provided) ---------------------
    if negate_enabled:
        if reference_score not in all_score_cols:
            available = ", ".join(all_score_cols)
            raise ValueError(
                f"Reference score column '{reference_score}' not found in "
                f"prediction tables. Available score columns: {available}"
            )
        print(
            f"  Computing correlations with reference '{reference_score}' ...",
            file=sys.stderr,
        )
        cols_to_negate = compute_negations(
            merged_pred, all_score_cols, reference_score,
        )
        if cols_to_negate:
            print(
                f"  Negating {len(cols_to_negate)} column(s): {cols_to_negate}",
                file=sys.stderr,
            )
            merged_pred = negate_scores(merged_pred, cols_to_negate)
        else:
            print("  All scores already aligned; no negation needed.",
                  file=sys.stderr)

    # --- Linker mode: join onto linker, aggregate, percentile ---------------
    if linker_mode:
        print(f"  Loading linker table {linker_uri} ...", file=sys.stderr)
        linker_pq = ensure_parquet(linker_uri, cache_dir)
        linker_lf = normalize_chrom_key(pl.scan_parquet(linker_pq))
        linker_lf = linker_lf.select(JOIN_KEYS + ["ensg"])

        print("  Left-joining linker onto merged predictions ...", file=sys.stderr)
        joined = linker_lf.join(
            merged_pred, on=JOIN_KEYS, how="left", coalesce=True,
        )

        print("  Aggregating scores by ENSG (mean + max) ...", file=sys.stderr)
        merged_pred, agg_score_cols = aggregate_by_gene(joined, all_score_cols)

        print(
            f"  Computing percentiles on {len(agg_score_cols)} aggregated columns ...",
            file=sys.stderr,
        )
        merged_pred = add_percentile_columns(merged_pred, agg_score_cols)

        if join_type == "pairwise" and non_anchor_cols:
            pct_renames: dict[str, str] = {}
            for suffix in ("mean", "max"):
                pct_renames[f"{anchor}_{suffix}_percentile"] = (
                    f"{anchor}_anchor_{suffix}_percentile"
                )
                for c in non_anchor_cols:
                    pct_renames[f"{c}_{suffix}_percentile"] = (
                        f"{c}_{suffix}_percentile_with_anchor"
                    )
                    pct_renames[f"{anchor}_pairwise_{c}_{suffix}_percentile"] = (
                        f"{anchor}_anchor_{suffix}_percentile_with_{c}"
                    )
            merged_pred = merged_pred.rename(pct_renames)

        drop_cols = agg_score_cols
    else:
        if percentile_order == "post":
            print("  Calculating percentiles AFTER pred merge ...", file=sys.stderr)
            merged_pred = add_percentile_columns(merged_pred, all_score_cols)

        if join_type == "pairwise" and non_anchor_cols and percentile_order != "none":
            pct_renames = {
                f"{anchor}_percentile": f"{anchor}_anchor_percentile",
            }
            for c in non_anchor_cols:
                pct_renames[f"{c}_percentile"] = f"{c}_percentile_with_anchor"
                pct_renames[f"{anchor}_pairwise_{c}_percentile"] = (
                    f"{anchor}_anchor_percentile_with_{c}"
                )
            merged_pred = merged_pred.rename(pct_renames)

        drop_cols = all_score_cols

    # --- Spatial smoothing (if requested) ---------------------------------
    if smooth:
        if percentile_order == "none":
            print(
                "  WARNING: --smooth requires percentile-ranked scores; "
                "--percentile_order is 'none' so smoothing will operate on "
                "raw scores instead.",
                file=sys.stderr,
            )
        if smooth_reference_dir is None:
            raise ValueError(
                "--smooth_reference_dir is required when --smooth is set."
            )
        cols_to_smooth = (
            [f"{c}_percentile" for c in all_score_cols]
            if percentile_order != "none"
            else all_score_cols
        )
        print(
            f"  Spatially smoothing {len(cols_to_smooth)} column(s) "
            f"(sigma={smooth_sigma} Å) ...",
            file=sys.stderr,
        )
        merged_pred = add_smoothed_columns(
            merged_pred,
            cols_to_smooth,
            reference_dir=smooth_reference_dir,
            sigma=smooth_sigma,
        )

    # --- Phase 3: Left-join eval onto pred --------------------------------
    print("  Left-joining eval and pred ...", file=sys.stderr)
    merged = merged_eval.join(
        merged_pred, on=eval_keys, how="left", coalesce=True,
    )

    # --- Drop raw score columns (if requested) -------------------------------
    has_percentiles = linker_mode or percentile_order != "none"
    if not keep_raw_scores and drop_cols and has_percentiles:
        print("  Dropping raw score columns ...", file=sys.stderr)
        merged = merged.drop(drop_cols)

    # --- Apply filter columns (if requested) --------------------------------
    if filter_uris:
        if linker_mode:
            print(
                "  WARNING: --filter_tables ignored in linker mode "
                "(gene-level output has no variant keys).",
                file=sys.stderr,
            )
        else:
            print(
                f"  Applying {len(filter_uris)} filter table(s) ...",
                file=sys.stderr,
            )
            merged = apply_filters(merged, filter_uris, cache_dir)

    # --- Write output -----------------------------------------------------
    print(f"  Writing output to {output_uri} ...", file=sys.stderr)
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
        choices=["inner", "outer", "pairwise"],
        default="inner",
        help=(
            "Join strategy for prediction tables (default: inner). "
            "'inner' keeps only rows where every score column is non-null "
            "across all prediction tables. 'outer' keeps all rows. "
            "'pairwise' inner-joins the column specified by --anchor with "
            "each other score column to form pairs, then outer-joins all "
            "pairs while preserving the full anchor column. "
            "Evaluation tables are always outer-joined."
        ),
    )
    parser.add_argument(
        "--anchor",
        default=None,
        help=(
            "Score column to use as the anchor for pairwise joining. "
            "Required when --join_type is 'pairwise'. The anchor column "
            "is inner-joined with each other score column to form pairs."
        ),
    )
    parser.add_argument(
        "--percentile_order",
        choices=["pre", "post", "none"],
        default="post",
        help=(
            "'pre' = compute percentile ranks on each prediction table before "
            "merging; 'post' = compute after merging prediction tables; "
            "'none' = skip percentile computation and retain raw scores "
            "(default: post). 'pre' is incompatible with score negation and "
            "will fall back to 'post' when --reference_score is set."
        ),
    )
    parser.add_argument(
        "--reference_score",
        default="AM",
        help=(
            "Reference score column for directional alignment. Scores with "
            "negative Pearson correlation to this column are negated so all "
            "scores point in the same direction. Set to 'none' to disable "
            "negation (default: AM)."
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
    parser.add_argument(
        "--output_table_fields",
        default=None,
        help=(
            "Comma-separated list of prediction score column names to include "
            "from the --prediction_tables. Score columns are filtered to this "
            "set BEFORE merging, so only the specified columns participate in "
            "joins. Derived columns (percentiles, pairwise pairs) are "
            "automatically generated for the included set. Evaluation label "
            "columns are not affected and are always included in full. The "
            "genomic key columns (chrom, pos, ref, alt) are always included. "
            "When using --join_type pairwise, the --anchor column must be "
            "listed. When score negation is enabled, the --reference_score "
            "column must be listed. If omitted, all score columns are included."
        ),
    )

    parser.add_argument(
        "--linker_table",
        default=None,
        help=(
            "URI of a linker parquet table mapping variant keys "
            "(chrom, pos, ref, alt) to ENSG gene IDs. When provided, "
            "prediction scores are left-joined onto the linker, aggregated "
            "to gene level (mean + max per ENSG), and percentiles are "
            "computed on the aggregated values. Eval tables are expected "
            "to be at gene level and are joined on ENSG."
        ),
    )
    parser.add_argument(
        "--smooth",
        action="store_true",
        default=False,
        help=(
            "Apply spatial smoothing to percentile-ranked scores using a "
            "Gaussian kernel over 3D protein structure distances. Requires "
            "--smooth_reference_dir."
        ),
    )
    parser.add_argument(
        "--smooth_reference_dir",
        default=None,
        help=(
            "Path to the sir-reference-data directory containing "
            "all_missense_variants_gr38.h5, pdb_pae_file_pos_guide.tsv, "
            "pdb_files/, and pae_files/. Required when --smooth is set. "
            "Reference data can be downloaded from: "
            "https://www.dropbox.com/scl/fi/t4it7sa9lkdx9maj0vois/sir-data.tar.gz"
            "?rlkey=flvsvmzyopj1cbn6gya0c3am0&st=uyk0l7iw&dl=0"
        ),
    )
    parser.add_argument(
        "--smooth_sigma",
        type=float,
        default=10.0,
        help=(
            "Gaussian kernel scale in Ångströms for spatial smoothing "
            "(default: 10.0). Only used when --smooth is set."
        ),
    )

    args = parser.parse_args()

    if args.join_type == "pairwise" and args.anchor is None:
        parser.error("--anchor is required when --join_type is 'pairwise'")
    if args.join_type != "pairwise" and args.anchor is not None:
        print(
            "  WARNING: --anchor is ignored when --join_type is not 'pairwise'.",
            file=sys.stderr,
        )

    ref = args.reference_score
    if ref and ref.lower() == "none":
        ref = None

    run_pipeline(
        prediction_uris=_parse_uri_list(args.prediction_tables),
        evaluation_uris=_parse_uri_list(args.evaluation_tables),
        output_uri=args.output,
        join_type=args.join_type,
        percentile_order=args.percentile_order,
        filter_uris=_parse_uri_list(args.filter_tables) if args.filter_tables else None,
        keep_raw_scores=args.keep_raw_scores,
        reference_score=ref,
        output_table_fields=(
            [f.strip() for f in args.output_table_fields.split(",") if f.strip()]
            if args.output_table_fields else None
        ),
        anchor=args.anchor,
        linker_uri=args.linker_table,
        smooth=args.smooth,
        smooth_reference_dir=args.smooth_reference_dir,
        smooth_sigma=args.smooth_sigma,
    )


if __name__ == "__main__":
    main()
