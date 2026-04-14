#!/usr/bin/env python3
"""
CLI tool for creating consolidated VSM (Variant Scoring Method) tables.

Reads prediction and evaluation tables from GCS or local paths, merges them
on genomic keys (chrom, pos, ref, alt), optionally computes percentile ranks
for prediction score columns, and writes a single output parquet file.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import time

from .table_io import ensure_parquet, normalize_chrom_key
from .merge import JOIN_KEYS
from .apply_filters import apply_filters
from .row_counts import RowCountsCollector, count_parquet_rows, write_report
from .paths import CACHE_DIR, parse_uri_list, derive_stem
from .pipeline import (
    load_inputs,
    merge_predictions,
    apply_gene_aggregation,
    apply_post_processing,
    join_and_write,
)


def run_pipeline(
    prediction_uris: list[str],
    evaluation_uris: list[str],
    output_uri: str,
    join_type: str = "inner",
    percentile_order: str = "post",
    filter_uris: list[str] | None = None,
    keep_raw_scores: bool = False,
    reference_score: str | None = "AM_score",
    output_table_fields: list[str] | None = None,
    anchor: str | None = None,
    linker_uri: str | None = None,
    smooth_order: str = "none",
    smooth_reference_dir: str | None = None,
    smooth_sigma: float = 10.0,
    aggregate_genes: bool = False,
    collapse_genes: bool = True,
    use_cache: bool = False,
    store_cache: bool = False,
    row_counts_output: str | None = None,
    percentile_thresholds: list[float] | None = None,
    subtable: str | None = None,
    precomputed_prediction: str | None = None,
    precomputed_evaluation: str | None = None,
    retain_raw_anchor: bool = False,
    dry_run: bool = False,
) -> None:
    """
    Core pipeline -- decoupled from CLI for reuse (e.g. future resources.json).

    Orchestrates phased functions from merge.pipeline to load inputs,
    merge predictions, optionally aggregate to gene level, apply
    post-processing, and write the final output.

    When *dry_run* is True, loads inputs and resolves schemas/column
    lists, prints what the output would contain, and returns without
    materializing or writing anything.
    """
    start = time.perf_counter()

    # --- Validate smooth_order -------------------------------------------
    if smooth_order != "none":
        if smooth_reference_dir is None:
            raise ValueError(
                "--smooth_reference_dir is required when --smooth_order is set."
            )
        if smooth_order == "pre" and percentile_order != "pre":
            raise ValueError(
                "--smooth_order 'pre' requires --percentile_order 'pre' "
                "(percentiles must be computed before smoothing)."
            )
        if smooth_order == "post" and percentile_order == "none":
            raise ValueError(
                "--smooth_order 'post' requires --percentile_order 'pre' or "
                "'post' (smoothing operates on percentile-ranked scores)."
            )

    # --- Row-count collector setup ----------------------------------------
    collector: RowCountsCollector | None = None
    if row_counts_output is not None:
        collector = RowCountsCollector()
        collector.record_config("output_uri", output_uri)
        collector.record_config("join_type", join_type)
        collector.record_config("percentile_order", percentile_order)
        if reference_score is not None:
            collector.record_config("reference_score", reference_score)
        if linker_uri is not None:
            collector.record_config("linker_table", derive_stem(linker_uri))
        collector.record_config("aggregate_genes", str(aggregate_genes))
        if aggregate_genes:
            collector.record_config("collapse_genes", str(collapse_genes))
        if smooth_order != "none":
            collector.record_config("smooth_order", smooth_order)
        if filter_uris:
            collector.record_config("filter_count", str(len(filter_uris)))

    if aggregate_genes and percentile_order != "post":
        print(
            "  WARNING: --percentile_order is ignored when --aggregate_genes "
            "is set; percentiles are computed on gene-level aggregates.",
            file=sys.stderr,
        )
    if aggregate_genes:
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

    # --- Fast path: both prediction and eval are precomputed --------------
    if precomputed_prediction is not None and precomputed_evaluation is not None:
        cache_dir = CACHE_DIR
        eval_keys = ["ensg"] if aggregate_genes else JOIN_KEYS
        pq_eval = ensure_parquet(
            precomputed_evaluation, cache_dir,
            use_cache=use_cache, store_cache=store_cache,
        )
        pq_pred = ensure_parquet(
            precomputed_prediction, cache_dir,
            use_cache=use_cache, store_cache=store_cache,
        )
        import polars as pl

        merged_eval = normalize_chrom_key(pl.scan_parquet(pq_eval))
        merged_pred = normalize_chrom_key(pl.scan_parquet(pq_pred))
        if linker_uri:
            print(f"  Left-joining linker {linker_uri} onto precomputed pred ...", file=sys.stderr)
            linker_lf = normalize_chrom_key(pl.scan_parquet(linker_uri))
            merged_pred = merged_pred.join(linker_lf, on=JOIN_KEYS, how="left", coalesce=True)
        print("  Left-joining precomputed eval and pred ...", file=sys.stderr)
        merged = merged_eval.join(
            merged_pred, on=eval_keys, how="left", coalesce=True,
        )
        if filter_uris:
            print(
                f"  Applying {len(filter_uris)} filter table(s) ...",
                file=sys.stderr,
            )
            merged = apply_filters(
                merged, filter_uris, cache_dir,
                use_cache=use_cache, store_cache=store_cache,
            )
        print(f"  Writing output to {output_uri} ...", file=sys.stderr)
        merged.sink_parquet(output_uri, compression="zstd")
        elapsed = time.perf_counter() - start
        print(f"Done in {elapsed:.1f}s.", file=sys.stderr)
        return

    # === Phased pipeline ==================================================

    inputs = load_inputs(
        prediction_uris,
        evaluation_uris,
        aggregate_genes=aggregate_genes,
        fields_set=fields_set,
        precomputed_evaluation=precomputed_evaluation,
        subtable=subtable,
        use_cache=use_cache,
        store_cache=store_cache,
        collector=collector,
    )

    # --- Dry-run: print schema summary and exit ---------------------------
    if dry_run:
        print("=== DRY RUN ===", file=sys.stderr)
        print(f"  Prediction score columns: {inputs.all_score_cols}", file=sys.stderr)
        if inputs.merged_eval is not None:
            eval_schema = inputs.merged_eval.collect_schema()
            eval_cols = [c for c in eval_schema.names() if c not in inputs.eval_keys]
            print(f"  Eval label columns: {eval_cols}", file=sys.stderr)
        else:
            print("  Eval: (none)", file=sys.stderr)
        print(f"  Join type: {join_type}", file=sys.stderr)
        print(f"  Percentile order: {percentile_order}", file=sys.stderr)
        print(f"  Aggregate genes: {aggregate_genes}", file=sys.stderr)
        print(f"  Output: {output_uri}", file=sys.stderr)
        return

    # --- Eval-only subtable early exit ------------------------------------
    if subtable == "eval":
        if inputs.merged_eval is None:
            raise ValueError("No evaluation tables to merge.")
        print(f"  Writing eval sub-table to {output_uri} ...", file=sys.stderr)
        inputs.merged_eval.sink_parquet(output_uri, compression="zstd")
        elapsed = time.perf_counter() - start
        print(f"Done in {elapsed:.1f}s.", file=sys.stderr)
        return

    # --- Merge predictions ------------------------------------------------
    pred = merge_predictions(
        inputs,
        join_type=join_type,
        percentile_order=percentile_order,
        reference_score=reference_score,
        anchor=anchor,
        smooth_order=smooth_order,
        smooth_reference_dir=smooth_reference_dir,
        smooth_sigma=smooth_sigma,
        collector=collector,
    )

    # --- Gene-level aggregation or variant-level post-processing ----------
    if aggregate_genes:
        pred = apply_gene_aggregation(
            pred,
            linker_uri=linker_uri,
            anchor=anchor,
            join_type=join_type,
            collapse_genes=collapse_genes,
            use_cache=use_cache,
            store_cache=store_cache,
            merged_eval=inputs.merged_eval,
            eval_keys=inputs.eval_keys,
            collector=collector,
        )
    else:
        pred = apply_post_processing(
            pred,
            join_type=join_type,
            percentile_order=percentile_order,
            anchor=anchor,
            aggregate_genes=aggregate_genes,
            smooth_order=smooth_order,
            smooth_reference_dir=smooth_reference_dir,
            smooth_sigma=smooth_sigma,
            retain_raw_anchor=retain_raw_anchor,
            linker_uri=linker_uri,
            use_cache=use_cache,
            store_cache=store_cache,
            merged_eval=inputs.merged_eval,
            eval_keys=inputs.eval_keys,
            collector=collector,
        )

    # --- Final join, filters, and output ----------------------------------
    result = join_and_write(
        pred,
        inputs.merged_eval,
        output_uri,
        eval_keys=inputs.eval_keys,
        subtable=subtable,
        keep_raw_scores=keep_raw_scores,
        aggregate_genes=aggregate_genes,
        percentile_order=percentile_order,
        join_type=join_type,
        filter_uris=filter_uris,
        percentile_thresholds=percentile_thresholds,
        use_cache=use_cache,
        store_cache=store_cache,
        collector=collector,
    )

    # --- Row count report -------------------------------------------------
    if collector is not None and row_counts_output is not None:
        try:
            collector.record("output", "rows", count_parquet_rows(output_uri))
        except Exception:
            pass
        write_report(collector, row_counts_output)

    elapsed = time.perf_counter() - start
    print(f"Done in {elapsed:.1f}s.", file=sys.stderr)


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="  %(message)s",
        stream=sys.stderr,
    )

    parser = argparse.ArgumentParser(
        description="Create a consolidated VSM table from prediction and evaluation data.",
    )

    parser.add_argument(
        "--prediction_tables",
        default=None,
        help=(
            "Comma-separated URIs (GCS or local) of prediction parquet/tsv "
            "files. An entry may also be a directory path, in which case all "
            "files with recognised extensions (.parquet, .tsv, .tsv.gz, "
            ".tsv.bgz) inside it are used. For GCS directories, append a "
            "trailing slash (e.g. gs://bucket/pred/). "
            "Required unless --subtable eval is set."
        ),
    )
    parser.add_argument(
        "--evaluation_tables",
        default=None,
        help=(
            "Comma-separated URIs (GCS or local) of evaluation parquet/tsv "
            "files. An entry may also be a directory path, in which case all "
            "files with recognised extensions (.parquet, .tsv, .tsv.gz, "
            ".tsv.bgz) inside it are used. For GCS directories, append a "
            "trailing slash (e.g. gs://bucket/eval/). "
            "Required unless --subtable pred is set."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination URI for the merged parquet file (GCS or local).",
    )
    parser.add_argument(
        "--subtable",
        choices=["pred", "eval"],
        default=None,
        help=(
            "Write only a prediction or evaluation sub-table instead of the "
            "full merged output. 'pred' merges and writes only prediction "
            "tables (with optional percentiles, negation, filters, linker). "
            "'eval' merges and writes only evaluation tables (outer join)."
        ),
    )
    parser.add_argument(
        "--precomputed_prediction",
        default=None,
        help=(
            "Path to a precomputed prediction parquet file (e.g. output of "
            "--subtable pred). Bypasses prediction table loading, merging, "
            "negation, percentiles, linker, and gene aggregation. Must be "
            "used together with --precomputed_evaluation."
        ),
    )
    parser.add_argument(
        "--precomputed_evaluation",
        default=None,
        help=(
            "Path to a precomputed evaluation parquet file (e.g. output of "
            "--subtable eval). Bypasses evaluation table loading and "
            "merging. Can be used alone (with fresh --prediction_tables) "
            "or together with --precomputed_prediction."
        ),
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
        default="AM_score",
        help=(
            "Reference score column for directional alignment. Scores with "
            "negative Pearson correlation to this column are negated so all "
            "scores point in the same direction. Set to 'none' to disable "
            "negation (default: AM_score)."
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
            "filter table. An entry may also be a directory path, in which "
            "case all files with recognised extensions inside it are used."
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
        "--smooth_order",
        choices=["pre", "post", "none"],
        default="none",
        help=(
            "When to apply spatial smoothing to percentile-ranked scores using "
            "a Gaussian kernel over 3D protein structure distances. "
            "'pre' = smooth each VSM's percentile-ranked scores before merging "
            "(requires --percentile_order pre); "
            "'post' = smooth after merging prediction tables (requires "
            "--percentile_order pre or post); "
            "'none' = skip smoothing (default). "
            "Requires --smooth_reference_dir when not 'none'."
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
    parser.add_argument(
        "--aggregate_genes",
        action="store_true",
        default=False,
        help=(
            "Aggregate prediction scores to gene level (mean + max per "
            "ensg). Requires that the prediction frame contains an 'ensg' "
            "column (via --linker_table or already present in prediction "
            "tables). Evaluation tables must be keyed on ensg when this "
            "flag is set."
        ),
    )
    parser.add_argument(
        "--collapse_genes",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "When --aggregate_genes is set, controls whether gene-level "
            "aggregation collapses to one row per gene (default) or "
            "preserves original variant rows with repeated mean/max values. "
            "Use --no-collapse_genes to preserve rows. Only valid with "
            "--aggregate_genes."
        ),
    )
    parser.add_argument(
        "--use_cache",
        action="store_true",
        default=False,
        help=(
            "Re-use previously cached TSV-to-Parquet conversions from "
            "$TMPDIR/vsm_table_cache/ when available (default: off)."
        ),
    )
    parser.add_argument(
        "--store_cache",
        action="store_true",
        default=False,
        help=(
            "Persist TSV-to-Parquet conversions in $TMPDIR/vsm_table_cache/ "
            "so subsequent runs can reuse them with --use_cache (default: off)."
        ),
    )
    parser.add_argument(
        "--row_counts",
        default=None,
        help=(
            "Write a comprehensive markdown row-count report tracking data "
            "through every pipeline stage. If the path ends with '.md' it "
            "is used as-is; otherwise it is treated as a directory and the "
            "filename is derived from the output table name."
        )
    )
    parser.add_argument(
        "--percentile_thresholds",
        default=None,
        help=(
            "Comma-separated percentile cutoffs (e.g. 0.85,0.9,0.95,0.995). "
            "Writes a TSV of raw score values at each threshold next to the "
            "output parquet (e.g. output.percentile_thresholds.tsv)."
        ),
    )
    parser.add_argument(
        "--retain_raw_anchor",
        action="store_true",
        default=False,
        help=(
            "When --join_type is 'pairwise', retain the original (raw, "
            "non-percentiled) anchor score column as '_raw_anchor_{name}' "
            "in the output. This column is required by merge-add-score for "
            "incremental pairwise score additions. Has no effect when "
            "--join_type is not 'pairwise'."
        ),
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        default=False,
        help=(
            "Resolve input schemas and column lists, print what the output "
            "would contain, then exit without materializing or writing data. "
            "Useful for validating arguments and previewing the pipeline."
        ),
    )

    args = parser.parse_args()

    if args.precomputed_prediction and not args.precomputed_evaluation:
        parser.error("--precomputed_prediction requires --precomputed_evaluation")
    if args.prediction_tables and args.precomputed_prediction:
        parser.error("--prediction_tables and --precomputed_prediction are mutually exclusive")
    if args.evaluation_tables and args.precomputed_evaluation:
        parser.error("--evaluation_tables and --precomputed_evaluation are mutually exclusive")

    if (args.subtable != "eval"
            and not args.prediction_tables
            and not args.precomputed_prediction):
        parser.error(
            "--prediction_tables is required (unless --subtable eval "
            "or --precomputed_prediction)"
        )
    if (args.subtable != "pred"
            and not args.evaluation_tables
            and not args.precomputed_evaluation):
        parser.error(
            "--evaluation_tables is required (unless --subtable pred "
            "or --precomputed_evaluation)"
        )

    if args.join_type == "pairwise" and args.anchor is None:
        parser.error("--anchor is required when --join_type is 'pairwise'")
    if args.join_type != "pairwise" and args.anchor is not None:
        print(
            "  WARNING: --anchor is ignored when --join_type is not 'pairwise'.",
            file=sys.stderr,
        )
    if not args.aggregate_genes and not args.collapse_genes:
        parser.error("--no-collapse_genes is only valid with --aggregate_genes")

    ref = args.reference_score
    if ref and ref.lower() == "none":
        ref = None

    row_counts_output = None
    if args.row_counts is not None:
        rc_path = args.row_counts
        if not rc_path.endswith(".md"):
            rc_path = os.path.join(rc_path, f"{derive_stem(args.output)}_counts.md")
        row_counts_output = rc_path
      
    pct_thresholds = None
    if args.percentile_thresholds is not None:
        pct_thresholds = [
            float(v.strip()) for v in args.percentile_thresholds.split(",") if v.strip()
        ]

    run_pipeline(
        prediction_uris=parse_uri_list(args.prediction_tables) if args.prediction_tables else [],
        evaluation_uris=parse_uri_list(args.evaluation_tables) if args.evaluation_tables else [],
        output_uri=args.output,
        join_type=args.join_type,
        percentile_order=args.percentile_order,
        filter_uris=parse_uri_list(args.filter_tables) if args.filter_tables else None,
        keep_raw_scores=args.keep_raw_scores,
        reference_score=ref,
        output_table_fields=(
            [f.strip() for f in args.output_table_fields.split(",") if f.strip()]
            if args.output_table_fields else None
        ),
        anchor=args.anchor,
        linker_uri=args.linker_table,
        smooth_order=args.smooth_order,
        smooth_reference_dir=args.smooth_reference_dir,
        smooth_sigma=args.smooth_sigma,
        aggregate_genes=args.aggregate_genes,
        collapse_genes=args.collapse_genes,
        use_cache=args.use_cache,
        store_cache=args.store_cache,
        row_counts_output=row_counts_output,
        percentile_thresholds=pct_thresholds,
        subtable=args.subtable,
        precomputed_prediction=args.precomputed_prediction,
        precomputed_evaluation=args.precomputed_evaluation,
        retain_raw_anchor=args.retain_raw_anchor,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
