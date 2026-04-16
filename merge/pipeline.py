"""
Phased pipeline functions for the merge pipeline.

Each function corresponds to a logical phase of the pipeline and returns
a structured dataclass carrying the LazyFrame plus metadata.  These
functions are designed to be called both from the CLI orchestrator
(create_vsm_table.run_pipeline) and from a future interactive session layer.

Progress messages go through the ``merge.pipeline`` logger at INFO level.
CLI callers can attach a StreamHandler(sys.stderr) to see progress;
library callers can silence output or redirect it via standard logging
configuration.
"""
from __future__ import annotations

import logging
import re
from collections import Counter
from functools import reduce

import polars as pl

from .table_io import ensure_parquet, normalize_chrom_key
from .merge import merge_tables, merge_tables_pairwise, aggregate_by_gene, JOIN_KEYS
from .percentile import add_percentile_columns
from .smooth import add_smoothed_columns
from .apply_filters import apply_filters
from .negate import compute_negations, negate_scores
from .row_counts import RowCountsCollector, count_parquet_rows, count_lazy
from .paths import CACHE_DIR, derive_stem, score_columns
from .types import LoadedInputs, MergedPrediction, PipelineResult

logger = logging.getLogger(__name__)


def load_inputs(
    prediction_uris: list[str],
    evaluation_uris: list[str],
    *,
    gene_prediction_uris: list[str] | None = None,
    aggregate_genes: bool = False,
    fields_set: set[str] | None = None,
    precomputed_evaluation: str | None = None,
    subtable: str | None = None,
    use_cache: bool = False,
    store_cache: bool = False,
    collector: RowCountsCollector | None = None,
) -> LoadedInputs:
    """Load and prepare prediction and evaluation tables.

    Scans prediction and evaluation URIs, resolves column lists, handles
    eval column collision renaming, and outer-joins evaluation frames.
    Also loads gene-level prediction tables (ensg-keyed) when provided.
    """
    cache_dir = CACHE_DIR
    eval_keys = ["ensg"] if aggregate_genes else JOIN_KEYS

    # --- Load prediction tables -------------------------------------------
    pred_frames: list[pl.LazyFrame] = []
    all_score_cols: list[str] = []

    for uri in prediction_uris:
        pq_path = ensure_parquet(
            uri, cache_dir, use_cache=use_cache, store_cache=store_cache,
        )
        if collector is not None:
            collector.record_input("pred", uri, count_parquet_rows(pq_path))
        lf = normalize_chrom_key(pl.scan_parquet(pq_path))
        s_cols = score_columns(lf)
        if fields_set is not None:
            s_cols = [c for c in s_cols if c in fields_set]
        if not s_cols:
            logger.info(
                "Prediction table %s: no matching score columns, skipping.",
                uri,
            )
            continue
        all_score_cols.extend(s_cols)
        lf = lf.select(JOIN_KEYS + s_cols)
        logger.info("Prediction table %s: score columns = %s", uri, s_cols)
        pred_frames.append(lf)

    if not pred_frames and subtable != "eval":
        specified = ", ".join(sorted(fields_set)) if fields_set else "(none)"
        raise ValueError(
            f"None of the --output_table_fields ({specified}) matched any "
            f"score columns in the prediction tables."
        )

    # --- Load evaluation tables -------------------------------------------
    eval_frames: list[pl.LazyFrame] = []

    for uri in evaluation_uris:
        pq_path = ensure_parquet(
            uri, cache_dir, use_cache=use_cache, store_cache=store_cache,
        )
        if collector is not None:
            collector.record_input("eval", uri, count_parquet_rows(pq_path))
        lf = normalize_chrom_key(pl.scan_parquet(pq_path))
        schema = lf.collect_schema()
        if aggregate_genes:
            label_cols = [
                c for c in schema.names() if c not in JOIN_KEYS and c != "ensg"
            ]
            if not label_cols:
                raise ValueError(
                    f"Evaluation table {uri} has no non-key label columns."
                )
            lf = lf.select(eval_keys + label_cols)
            logger.info("Evaluation table %s: labels = %s", uri, label_cols)
        else:
            stem = derive_stem(uri)
            is_pos_cols = [c for c in schema.names() if "is_pos" in c]
            if not is_pos_cols:
                raise ValueError(
                    f"Evaluation table {uri} has no column containing 'is_pos'."
                )
            is_pos_col = "is_pos" if "is_pos" in is_pos_cols else is_pos_cols[0]
            target_name = f"is_pos_{stem}"
            lf = lf.select(JOIN_KEYS + [is_pos_col]).rename(
                {is_pos_col: target_name}
            )
            logger.info(
                "Evaluation table %s: %s -> %s", uri, is_pos_col, target_name
            )
        eval_frames.append(lf)

    # --- Auto-suffix colliding eval label columns -------------------------
    if aggregate_genes:
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
            logger.info("Renaming colliding label columns: %s", collisions)
            for i, uri in enumerate(evaluation_uris):
                stem = derive_stem(uri)
                frame_cols = set(eval_frames[i].collect_schema().names())
                renames = {
                    col: f"{col}_{stem}"
                    for col in collisions
                    if col in frame_cols
                }
                if renames:
                    eval_frames[i] = eval_frames[i].rename(renames)

    # --- Outer-join eval tables -------------------------------------------
    merged_eval: pl.LazyFrame | None = None
    if eval_frames:
        logger.info(
            "Merging %d eval table(s) (outer join) ...", len(eval_frames)
        )
        if aggregate_genes:

            def _eval_join(
                left: pl.LazyFrame, right: pl.LazyFrame
            ) -> pl.LazyFrame:
                return left.join(
                    right, on=eval_keys, how="full", coalesce=True
                )

            if len(eval_frames) == 1:
                merged_eval = eval_frames[0]
            else:
                merged_eval = reduce(_eval_join, eval_frames)
        else:
            merged_eval = merge_tables(eval_frames, join_type="outer")

        if collector is not None:
            collector.record("eval_merge", "rows", count_lazy(merged_eval))

    if precomputed_evaluation is not None:
        pq_eval = ensure_parquet(
            precomputed_evaluation,
            cache_dir,
            use_cache=use_cache,
            store_cache=store_cache,
        )
        merged_eval = normalize_chrom_key(pl.scan_parquet(pq_eval))

    # --- Load gene-level prediction tables --------------------------------
    gene_pred_frames: list[pl.LazyFrame] = []
    gene_score_cols: list[str] = []

    for uri in gene_prediction_uris or []:
        pq_path = ensure_parquet(
            uri, cache_dir, use_cache=use_cache, store_cache=store_cache,
        )
        if collector is not None:
            collector.record_input("gene_pred", uri, count_parquet_rows(pq_path))
        lf = pl.scan_parquet(pq_path)
        s_cols = score_columns(lf)
        if fields_set is not None:
            s_cols = [c for c in s_cols if c in fields_set]
        if not s_cols:
            logger.info(
                "Gene prediction table %s: no matching score columns, skipping.",
                uri,
            )
            continue
        gene_score_cols.extend(s_cols)
        lf = lf.select(["ensg"] + s_cols)
        logger.info(
            "Gene prediction table %s: score columns = %s", uri, s_cols
        )
        gene_pred_frames.append(lf)

    return LoadedInputs(
        pred_frames=pred_frames,
        all_score_cols=all_score_cols,
        merged_eval=merged_eval,
        eval_keys=eval_keys,
        fields_set=fields_set,
        gene_pred_frames=gene_pred_frames,
        gene_score_cols=gene_score_cols,
    )


def merge_predictions(
    inputs: LoadedInputs,
    *,
    join_type: str = "inner",
    percentile_order: str = "post",
    reference_score: str | None = "AM_score",
    anchor: str | None = None,
    smooth_order: str = "none",
    smooth_reference_dir: str | None = None,
    smooth_sigma: float = 10.0,
    collector: RowCountsCollector | None = None,
) -> MergedPrediction:
    """Merge prediction frames with optional negation and pre-merge percentiles.

    Returns a MergedPrediction carrying the merged LazyFrame plus metadata
    about score columns, pairwise structure, and negation results.
    """
    pred_frames = list(inputs.pred_frames)
    all_score_cols = list(inputs.all_score_cols)
    negate_enabled = reference_score is not None
    negated_cols: list[str] = []

    if percentile_order == "pre" and (negate_enabled or join_type == "pairwise"):
        reason = (
            "score negation (requires merged frame)"
            if negate_enabled
            else "pairwise join (restructures columns)"
        )
        logger.warning(
            "--percentile_order 'pre' is incompatible with %s; "
            "falling back to 'post'.",
            reason,
        )
        percentile_order = "post"

    # --- Pre-merge negation discovery pass --------------------------------
    if negate_enabled and percentile_order == "pre":
        if reference_score not in all_score_cols:
            available = ", ".join(all_score_cols)
            raise ValueError(
                f"Reference score column '{reference_score}' not found in "
                f"prediction tables. Available score columns: {available}"
            )
        logger.info(
            "Discovery pass: computing correlations with '%s' ...",
            reference_score,
        )
        temp_merged = merge_tables(pred_frames, join_type="inner")
        cols_to_negate = compute_negations(
            temp_merged, all_score_cols, reference_score,
        )
        if cols_to_negate:
            logger.info(
                "Negating %d column(s) on individual frames: %s",
                len(cols_to_negate),
                cols_to_negate,
            )
            pred_frames = [
                negate_scores(
                    lf,
                    [c for c in cols_to_negate if c in lf.collect_schema().names()],
                )
                for lf in pred_frames
            ]
            negated_cols = cols_to_negate
        else:
            logger.info("All scores already aligned; no negation needed.")
        negate_enabled = False

    if percentile_order == "pre":
        logger.info("Calculating percentiles BEFORE pred merge ...")
        new_frames = []
        for lf in pred_frames:
            s_cols = score_columns(lf)
            lf = add_percentile_columns(lf, s_cols)
            if smooth_order == "pre":
                assert smooth_reference_dir is not None
                pct_cols = [f"{c}_percentile" for c in s_cols]
                logger.info(
                    "Spatially smoothing %d column(s) before merge "
                    "(sigma=%.1f Å) ...",
                    len(pct_cols),
                    smooth_sigma,
                )
                lf = add_smoothed_columns(
                    lf,
                    pct_cols,
                    reference_dir=smooth_reference_dir,
                    sigma=smooth_sigma,
                )
            new_frames.append(lf)
        pred_frames = new_frames

    logger.info(
        "Merging %d pred table(s) (%s join) ...", len(pred_frames), join_type
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

    if collector is not None:
        collector.record("pred_merge", "rows", count_lazy(merged_pred))

    if join_type == "inner" and all_score_cols:
        logger.info("Dropping rows with any null score column ...")
        merged_pred = merged_pred.drop_nulls(subset=all_score_cols)
        if collector is not None:
            collector.record(
                "pred_merge", "after_null_drop", count_lazy(merged_pred)
            )

    # --- Post-merge negation ----------------------------------------------
    if negate_enabled:
        if reference_score not in all_score_cols:
            available = ", ".join(all_score_cols)
            raise ValueError(
                f"Reference score column '{reference_score}' not found in "
                f"prediction tables. Available score columns: {available}"
            )
        logger.info(
            "Computing correlations with reference '%s' ...", reference_score
        )
        cols_to_negate = compute_negations(
            merged_pred, all_score_cols, reference_score,
        )
        if cols_to_negate:
            logger.info(
                "Negating %d column(s): %s", len(cols_to_negate), cols_to_negate
            )
            merged_pred = negate_scores(merged_pred, cols_to_negate)
            negated_cols = cols_to_negate
        else:
            logger.info("All scores already aligned; no negation needed.")

    return MergedPrediction(
        frame=merged_pred,
        score_cols=all_score_cols,
        non_anchor_cols=non_anchor_cols,
        drop_cols=list(all_score_cols),
        negated_cols=negated_cols,
    )


def _load_linker(
    linker_uri: str,
    *,
    use_cache: bool = False,
    store_cache: bool = False,
    merged_eval: pl.LazyFrame | None = None,
    eval_keys: list[str] | None = None,
    collector: RowCountsCollector | None = None,
) -> pl.LazyFrame:
    """Load and prepare a linker table, recording diagnostics if a collector is active."""
    cache_dir = CACHE_DIR
    if eval_keys is None:
        eval_keys = JOIN_KEYS

    logger.info("Loading linker table %s ...", linker_uri)
    pq = ensure_parquet(
        linker_uri, cache_dir, use_cache=use_cache, store_cache=store_cache,
    )
    linker_lf = (
        normalize_chrom_key(pl.scan_parquet(pq))
        .select(JOIN_KEYS + ["ensg"])
        .unique()
    )

    if collector is not None:
        collector.record_input("linker", linker_uri, count_parquet_rows(pq))
        linker_keys = linker_lf.select(JOIN_KEYS).unique()
        collector.record("linker", "unique_variants", count_lazy(linker_keys))
        collector.record(
            "linker",
            "unique_ensgs",
            linker_lf.select("ensg")
            .unique()
            .select(pl.len())
            .collect()
            .item(),
        )
        if merged_eval is not None:
            eval_not_in_linker = (
                count_lazy(
                    merged_eval.join(linker_keys, on=eval_keys, how="anti")
                )
                if eval_keys == JOIN_KEYS
                else count_lazy(
                    merged_eval.join(
                        linker_lf.select("ensg").unique(),
                        on=["ensg"],
                        how="anti",
                    )
                )
            )
            collector.record(
                "coverage", "eval_not_in_linker", eval_not_in_linker
            )
            eval_total = collector.stages["eval_merge"]["rows"]
            collector.record(
                "coverage", "eval_in_linker", eval_total - eval_not_in_linker
            )

    return linker_lf


def _join_linker(
    merged_pred: pl.LazyFrame,
    linker_lf: pl.LazyFrame,
    *,
    collector: RowCountsCollector | None = None,
) -> pl.LazyFrame:
    """Left-join linker onto predictions and record diagnostics."""
    logger.info("Left-joining linker onto merged predictions ...")
    merged_pred = merged_pred.join(
        linker_lf, on=JOIN_KEYS, how="left", coalesce=True,
    )

    if collector is not None:
        collector.record("linker", "pred_after_join", count_lazy(merged_pred))
        pred_no_ensg = (
            merged_pred.filter(pl.col("ensg").is_null())
            .select(pl.len())
            .collect()
            .item()
        )
        collector.record("linker", "pred_no_ensg", pred_no_ensg)

    return merged_pred


def apply_gene_aggregation(
    pred: MergedPrediction,
    *,
    linker_uri: str | None = None,
    anchor: str | None = None,
    join_type: str = "inner",
    collapse_genes: bool = True,
    use_cache: bool = False,
    store_cache: bool = False,
    merged_eval: pl.LazyFrame | None = None,
    eval_keys: list[str] | None = None,
    gene_pred_frames: list[pl.LazyFrame] | None = None,
    gene_score_cols: list[str] | None = None,
    collector: RowCountsCollector | None = None,
) -> MergedPrediction:
    """Apply linker join and gene-level aggregation to merged predictions.

    Handles loading the linker table, joining it onto predictions,
    aggregating scores by gene (mean + max), computing gene-level
    percentiles, and renaming pairwise percentile columns.

    If ``gene_pred_frames`` is provided, those ensg-keyed frames are
    outer-joined after aggregation but before percentile computation so
    that all gene-level scores share the same percentile denominator.
    """
    merged_pred = pred.frame
    all_score_cols = list(pred.score_cols)
    non_anchor_cols = list(pred.non_anchor_cols)
    if eval_keys is None:
        eval_keys = JOIN_KEYS

    # --- Linker join ------------------------------------------------------
    if linker_uri:
        linker_lf = _load_linker(
            linker_uri,
            use_cache=use_cache,
            store_cache=store_cache,
            merged_eval=merged_eval,
            eval_keys=eval_keys,
            collector=collector,
        )
        merged_pred = _join_linker(
            merged_pred, linker_lf, collector=collector
        )

    # --- Gene-level aggregation -------------------------------------------
    pred_schema = merged_pred.collect_schema()
    if "ensg" not in pred_schema.names():
        raise ValueError(
            "Gene-level aggregation requires an 'ensg' column. "
            "Provide --linker_table or ensure prediction tables "
            "contain 'ensg'."
        )

    logger.info("Aggregating scores by ensg (mean + max) ...")
    merged_pred, agg_score_cols = aggregate_by_gene(
        merged_pred, all_score_cols, collapse=collapse_genes,
    )
    if collector is not None:
        collector.record("gene_agg", "rows", count_lazy(merged_pred))

    # --- Join native gene-level prediction tables -------------------------
    _gene_frames = gene_pred_frames or []
    _gene_cols = list(gene_score_cols or [])
    if _gene_frames:
        for gf in _gene_frames:
            merged_pred = merged_pred.join(
                gf, on="ensg", how="full", coalesce=True,
            )
        agg_score_cols.extend(_gene_cols)
        logger.info(
            "Joined %d gene-level prediction table(s): %s",
            len(_gene_frames),
            _gene_cols,
        )
        if collector is not None:
            collector.record(
                "gene_pred_join", "rows", count_lazy(merged_pred)
            )

    logger.info(
        "Computing percentiles on %d aggregated columns ...",
        len(agg_score_cols),
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

    return MergedPrediction(
        frame=merged_pred,
        score_cols=all_score_cols,
        non_anchor_cols=non_anchor_cols,
        drop_cols=agg_score_cols,
        negated_cols=pred.negated_cols,
    )


def apply_post_processing(
    pred: MergedPrediction,
    *,
    join_type: str = "inner",
    percentile_order: str = "post",
    anchor: str | None = None,
    aggregate_genes: bool = False,
    smooth_order: str = "none",
    smooth_reference_dir: str | None = None,
    smooth_sigma: float = 10.0,
    retain_raw_anchor: bool = False,
    linker_uri: str | None = None,
    use_cache: bool = False,
    store_cache: bool = False,
    merged_eval: pl.LazyFrame | None = None,
    eval_keys: list[str] | None = None,
    collector: RowCountsCollector | None = None,
) -> MergedPrediction:
    """Apply variant-level post-processing: percentiles, smoothing, anchor retention.

    Also handles the linker join for non-gene-aggregated tables (when a
    linker table is provided but aggregate_genes is False).
    """
    merged_pred = pred.frame
    all_score_cols = list(pred.score_cols)
    non_anchor_cols = list(pred.non_anchor_cols)
    drop_cols = list(pred.drop_cols)
    if eval_keys is None:
        eval_keys = JOIN_KEYS

    # --- Linker join (non-aggregated path) --------------------------------
    if linker_uri:
        linker_lf = _load_linker(
            linker_uri,
            use_cache=use_cache,
            store_cache=store_cache,
            merged_eval=merged_eval,
            eval_keys=eval_keys,
            collector=collector,
        )
        merged_pred = _join_linker(
            merged_pred, linker_lf, collector=collector
        )

    # --- Post-merge percentiles -------------------------------------------
    if percentile_order == "post":
        logger.info("Calculating percentiles AFTER pred merge ...")
        merged_pred = add_percentile_columns(merged_pred, all_score_cols)

    if (
        join_type == "pairwise"
        and non_anchor_cols
        and percentile_order != "none"
    ):
        pct_renames = {
            f"{anchor}_percentile": f"{anchor}_anchor_percentile",
        }
        for c in non_anchor_cols:
            pct_renames[f"{c}_percentile"] = f"{c}_percentile_with_anchor"
            pct_renames[f"{anchor}_pairwise_{c}_percentile"] = (
                f"{anchor}_anchor_percentile_with_{c}"
            )
        merged_pred = merged_pred.rename(pct_renames)

    drop_cols = list(all_score_cols)

    # --- Post-merge smoothing ---------------------------------------------
    if smooth_order == "post":
        assert smooth_reference_dir is not None
        cols_to_smooth = [f"{c}_percentile" for c in all_score_cols]
        logger.info(
            "Spatially smoothing %d column(s) after merge (sigma=%.1f Å) ...",
            len(cols_to_smooth),
            smooth_sigma,
        )
        merged_pred = add_smoothed_columns(
            merged_pred,
            cols_to_smooth,
            reference_dir=smooth_reference_dir,
            sigma=smooth_sigma,
        )

    # --- Raw anchor retention ---------------------------------------------
    if retain_raw_anchor and join_type == "pairwise" and anchor:
        raw_anchor_name = f"_raw_anchor_{anchor}"
        merged_pred = merged_pred.rename({anchor: raw_anchor_name})
        all_score_cols = [
            raw_anchor_name if c == anchor else c for c in all_score_cols
        ]
        drop_cols = [c for c in drop_cols if c != anchor]

    return MergedPrediction(
        frame=merged_pred,
        score_cols=all_score_cols,
        non_anchor_cols=non_anchor_cols,
        drop_cols=drop_cols,
        negated_cols=pred.negated_cols,
    )


def join_and_write(
    pred: MergedPrediction,
    merged_eval: pl.LazyFrame | None,
    output_uri: str,
    *,
    eval_keys: list[str] | None = None,
    subtable: str | None = None,
    keep_raw_scores: bool = False,
    aggregate_genes: bool = False,
    percentile_order: str = "post",
    join_type: str = "inner",
    filter_uris: list[str] | None = None,
    percentile_thresholds: list[float] | None = None,
    use_cache: bool = False,
    store_cache: bool = False,
    collector: RowCountsCollector | None = None,
) -> PipelineResult:
    """Join eval onto pred, apply filters, and write output.

    Returns a PipelineResult with the final frame and metadata.
    """
    cache_dir = CACHE_DIR
    merged_pred = pred.frame
    all_score_cols = pred.score_cols
    drop_cols = pred.drop_cols
    if eval_keys is None:
        eval_keys = JOIN_KEYS

    # --- Join eval and pred -----------------------------------------------
    if subtable == "pred" or merged_eval is None:
        merged = merged_pred
    else:
        if collector is not None:
            eval_not_in_pred = count_lazy(
                merged_eval.join(merged_pred, on=eval_keys, how="anti")
            )
            eval_total = collector.stages["eval_merge"]["rows"]
            collector.record("coverage", "eval_not_in_pred", eval_not_in_pred)
            collector.record(
                "coverage", "eval_in_pred", eval_total - eval_not_in_pred
            )

        logger.info("Left-joining eval and pred ...")
        merged = merged_eval.join(
            merged_pred, on=eval_keys, how="left", coalesce=True,
        )

    # --- Per-score null counts --------------------------------------------
    if collector is not None:
        score_cols_to_check = (
            [f"{c}_percentile" for c in drop_cols]
            if aggregate_genes
            else (
                [f"{c}_percentile" for c in all_score_cols]
                if percentile_order != "none"
                else list(all_score_cols)
            )
        )
        present = set(merged.collect_schema().names())
        score_cols_to_check = [c for c in score_cols_to_check if c in present]
        if score_cols_to_check:
            null_counts = (
                merged.select(
                    [
                        pl.col(c).null_count().alias(c)
                        for c in score_cols_to_check
                    ]
                )
                .collect()
                .row(0, named=True)
            )
            for col, n in null_counts.items():
                collector.record("per_score_nulls", col, n)

    # --- Drop raw score columns -------------------------------------------
    has_percentiles = aggregate_genes or percentile_order != "none"
    if not keep_raw_scores and drop_cols and has_percentiles:
        logger.info("Dropping raw score columns ...")
        merged = merged.drop(drop_cols)

    # --- Apply filters ----------------------------------------------------
    if filter_uris:
        logger.info("Applying %d filter table(s) ...", len(filter_uris))
        merged = apply_filters(
            merged,
            filter_uris,
            cache_dir,
            use_cache=use_cache,
            store_cache=store_cache,
        )

    # --- Percentile thresholds TSV ----------------------------------------
    if percentile_thresholds:
        threshold_score_cols = drop_cols if aggregate_genes else all_score_cols
        thresholds_df = (
            merged_pred.select(
                [
                    pl.col(c)
                    .quantile(q, interpolation="nearest")
                    .alias(f"{c}_p{q}")
                    for c in threshold_score_cols
                    for q in percentile_thresholds
                ]
            ).collect()
        )

        rows: list[dict[str, object]] = []
        for c in threshold_score_cols:
            row: dict[str, object] = {"score": c}
            for q in percentile_thresholds:
                row[f"p{q}"] = thresholds_df[f"{c}_p{q}"][0]
            rows.append(row)
        result_df = pl.DataFrame(rows)

        tsv_path = (
            re.sub(r"\.parquet$", "", output_uri)
            + ".percentile_thresholds.tsv"
        )
        result_df.write_csv(tsv_path, separator="\t")
        logger.info("Percentile thresholds written to %s", tsv_path)

    # --- Write output -----------------------------------------------------
    logger.info("Writing output to %s ...", output_uri)
    merged.sink_parquet(output_uri, compression="zstd")

    # --- Determine eval columns for result metadata -----------------------
    eval_cols: list[str] = []
    if merged_eval is not None and subtable != "pred":
        eval_schema = merged_eval.collect_schema()
        eval_cols = [c for c in eval_schema.names() if c not in eval_keys]

    return PipelineResult(
        frame=merged,
        score_cols=list(all_score_cols),
        eval_cols=eval_cols,
        join_type=join_type,
        percentile_order=percentile_order,
        has_gene_aggregation=aggregate_genes,
    )
