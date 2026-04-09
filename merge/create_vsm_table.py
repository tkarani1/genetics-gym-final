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
import re
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
from row_counts import RowCountsCollector, count_parquet_rows, count_lazy, write_report


_KNOWN_EXTENSIONS = (".tsv.bgz", ".tsv.gz", ".tsv", ".parquet")


def _expand_path(entry: str) -> list[str]:
    """If *entry* is a directory, return all files inside it with a known
    table extension.  Works for both local paths and GCS URIs (the latter
    require a trailing ``/`` to signal directory intent)."""
    if entry.startswith("gs://"):
        if not entry.endswith("/"):
            return [entry]
        import gcsfs
        fs = gcsfs.GCSFileSystem()
        blobs = fs.ls(entry, detail=False)
        found = sorted(
            f"gs://{b}" for b in blobs
            if any(b.lower().endswith(ext) for ext in _KNOWN_EXTENSIONS)
        )
        if not found:
            print(
                f"  WARNING: GCS directory '{entry}' contains no files with "
                f"recognised extensions {_KNOWN_EXTENSIONS}.",
                file=sys.stderr,
            )
        return found

    if os.path.isdir(entry):
        found = sorted(
            os.path.join(entry, name)
            for name in os.listdir(entry)
            if any(name.lower().endswith(ext) for ext in _KNOWN_EXTENSIONS)
        )
        if not found:
            print(
                f"  WARNING: Directory '{entry}' contains no files with "
                f"recognised extensions {_KNOWN_EXTENSIONS}.",
                file=sys.stderr,
            )
        return found

    return [entry]


def _parse_uri_list(raw: str) -> list[str]:
    entries = [s.strip() for s in raw.split(",") if s.strip()]
    expanded: list[str] = []
    for entry in entries:
        expanded.extend(_expand_path(entry))
    return expanded


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
    aggregate_genes: bool = False,
    collapse_genes: bool = True,
    use_cache: bool = False,
    store_cache: bool = False,
    row_counts_output: str | None = None,
    percentile_thresholds: list[float] | None = None,
    subtable: str | None = None,
    precomputed_prediction: str | None = None,
    precomputed_evaluation: str | None = None,
) -> None:
    """
    Core pipeline -- decoupled from CLI for reuse (e.g. future resources.json).
    """
    start = time.perf_counter()

    collector: RowCountsCollector | None = None
    if row_counts_output is not None:
        collector = RowCountsCollector()
        collector.record_config("output_uri", output_uri)
        collector.record_config("join_type", join_type)
        collector.record_config("percentile_order", percentile_order)
        if reference_score is not None:
            collector.record_config("reference_score", reference_score)
        if linker_uri is not None:
            collector.record_config("linker_table", _derive_stem(linker_uri))
        collector.record_config("aggregate_genes", str(aggregate_genes))
        if aggregate_genes:
            collector.record_config("collapse_genes", str(collapse_genes))
        if smooth:
            collector.record_config("smooth", "True")
        if filter_uris:
            collector.record_config("filter_count", str(len(filter_uris)))

    cache_dir = os.path.join(tempfile.gettempdir(), "vsm_table_cache")

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

    # --- Fast path: both prediction and eval are precomputed ---------------
    if precomputed_prediction is not None and precomputed_evaluation is not None:
        eval_keys = ["ensg"] if aggregate_genes else JOIN_KEYS
        pq_eval = ensure_parquet(
            precomputed_evaluation, cache_dir,
            use_cache=use_cache, store_cache=store_cache,
        )
        pq_pred = ensure_parquet(
            precomputed_prediction, cache_dir,
            use_cache=use_cache, store_cache=store_cache,
        )
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

    # --- Load prediction tables -------------------------------------------
    pred_frames: list[pl.LazyFrame] = []
    all_score_cols: list[str] = []

    for uri in prediction_uris:
        pq_path = ensure_parquet(uri, cache_dir, use_cache=use_cache, store_cache=store_cache)
        if collector is not None:
            collector.record_input("pred", uri, count_parquet_rows(pq_path))
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

    if not pred_frames and subtable != "eval":
        specified = ", ".join(output_table_fields) if output_table_fields else "(none)"
        raise ValueError(
            f"None of the --output_table_fields ({specified}) matched any "
            f"score columns in the prediction tables."
        )

    # --- Load evaluation tables -------------------------------------------
    eval_keys = ["ensg"] if aggregate_genes else JOIN_KEYS
    eval_frames: list[pl.LazyFrame] = []

    for uri in evaluation_uris:
        pq_path = ensure_parquet(uri, cache_dir, use_cache=use_cache, store_cache=store_cache)
        if collector is not None:
            collector.record_input("eval", uri, count_parquet_rows(pq_path))
        lf = normalize_chrom_key(pl.scan_parquet(pq_path))
        schema = lf.collect_schema()
        if aggregate_genes:
            label_cols = [c for c in schema.names()
                          if c not in JOIN_KEYS and c != "ensg"]
            if not label_cols:
                raise ValueError(
                    f"Evaluation table {uri} has no non-key label columns."
                )
            lf = lf.select(eval_keys + label_cols)
            print(
                f"  Evaluation table {uri}: labels = {label_cols}",
                file=sys.stderr,
            )
        else:
            stem = _derive_stem(uri)
            is_pos_cols = [c for c in schema.names() if "is_pos" in c]
            if not is_pos_cols:
                raise ValueError(
                    f"Evaluation table {uri} has no column containing 'is_pos'."
                )
            is_pos_col = "is_pos" if "is_pos" in is_pos_cols else is_pos_cols[0]
            target_name = f"is_pos_{stem}"
            lf = lf.select(JOIN_KEYS + [is_pos_col]).rename({is_pos_col: target_name})
            print(
                f"  Evaluation table {uri}: {is_pos_col} -> {target_name}",
                file=sys.stderr,
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

    # --- Phase 1: Outer-join eval tables ------------------------------------
    merged_eval: pl.LazyFrame | None = None
    if eval_frames:
        print(
            f"  Merging {len(eval_frames)} eval table(s) (outer join) ...",
            file=sys.stderr,
        )
        if aggregate_genes:
            def _eval_join(left: pl.LazyFrame, right: pl.LazyFrame) -> pl.LazyFrame:
                return left.join(right, on=eval_keys, how="outer", coalesce=True)

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
            precomputed_evaluation, cache_dir,
            use_cache=use_cache, store_cache=store_cache,
        )
        merged_eval = normalize_chrom_key(pl.scan_parquet(pq_eval))

    if subtable == "eval":
        if merged_eval is None:
            raise ValueError("No evaluation tables to merge.")
        print(f"  Writing eval sub-table to {output_uri} ...", file=sys.stderr)
        merged_eval.sink_parquet(output_uri, compression="zstd")
        elapsed = time.perf_counter() - start
        print(f"Done in {elapsed:.1f}s.", file=sys.stderr)
        return

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

    # --- Pre-merge negation discovery pass --------------------------------
    if negate_enabled and percentile_order == "pre":
        if reference_score not in all_score_cols:
            available = ", ".join(all_score_cols)
            raise ValueError(
                f"Reference score column '{reference_score}' not found in "
                f"prediction tables. Available score columns: {available}"
            )
        print(
            f"  Discovery pass: computing correlations with "
            f"'{reference_score}' ...",
            file=sys.stderr,
        )
        temp_merged = merge_tables(pred_frames, join_type="inner")
        cols_to_negate = compute_negations(
            temp_merged, all_score_cols, reference_score,
        )
        if cols_to_negate:
            print(
                f"  Negating {len(cols_to_negate)} column(s) on individual "
                f"frames: {cols_to_negate}",
                file=sys.stderr,
            )
            pred_frames = [
                negate_scores(
                    lf,
                    [c for c in cols_to_negate
                     if c in lf.collect_schema().names()],
                )
                for lf in pred_frames
            ]
        else:
            print("  All scores already aligned; no negation needed.",
                  file=sys.stderr)
        negate_enabled = False

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

    if collector is not None:
        collector.record("pred_merge", "rows", count_lazy(merged_pred))

    if join_type == "inner" and all_score_cols:
        print("  Dropping rows with any null score column ...", file=sys.stderr)
        merged_pred = merged_pred.drop_nulls(subset=all_score_cols)
        if collector is not None:
            collector.record("pred_merge", "after_null_drop", count_lazy(merged_pred))

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

    # --- Linker table: add ensg column to predictions -----------------------
    linker_lf: pl.LazyFrame | None = None
    if linker_uri:
        print(f"  Loading linker table {linker_uri} ...", file=sys.stderr)
        pq = ensure_parquet(linker_uri, cache_dir, use_cache=use_cache, store_cache=store_cache)
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
                "linker", "unique_ensgs",
                linker_lf.select("ensg").unique().select(pl.len()).collect().item(),
            )
            if merged_eval is not None:
                eval_not_in_linker = count_lazy(
                    merged_eval.join(linker_keys, on=eval_keys, how="anti")
                ) if eval_keys == JOIN_KEYS else count_lazy(
                    merged_eval.join(
                        linker_lf.select("ensg").unique(), on=["ensg"], how="anti",
                    )
                )
                collector.record("coverage", "eval_not_in_linker", eval_not_in_linker)
                eval_total = collector.stages["eval_merge"]["rows"]
                collector.record("coverage", "eval_in_linker", eval_total - eval_not_in_linker)

        print("  Left-joining linker onto merged predictions ...", file=sys.stderr)
        merged_pred = merged_pred.join(
            linker_lf, on=JOIN_KEYS, how="left", coalesce=True,
        )

        if collector is not None:
            collector.record("linker", "pred_after_join", count_lazy(merged_pred))
            pred_no_ensg = merged_pred.filter(pl.col("ensg").is_null()).select(pl.len()).collect().item()
            collector.record("linker", "pred_no_ensg", pred_no_ensg)

    # --- Gene-level aggregation (if requested) ----------------------------
    if aggregate_genes:
        pred_schema = merged_pred.collect_schema()
        if "ensg" not in pred_schema.names():
            raise ValueError(
                "Gene-level aggregation requires an 'ensg' column. "
                "Provide --linker_table or ensure prediction tables "
                "contain 'ensg'."
            )

        print("  Aggregating scores by ensg (mean + max) ...", file=sys.stderr)
        merged_pred, agg_score_cols = aggregate_by_gene(
            merged_pred, all_score_cols, collapse=collapse_genes,
        )
        if collector is not None:
            collector.record("gene_agg", "rows", count_lazy(merged_pred))

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

    # --- Retain raw anchor for future pairwise additions --------------------
    if join_type == "pairwise" and anchor and not aggregate_genes:
        raw_anchor_name = f"_raw_anchor_{anchor}"
        merged_pred = merged_pred.rename({anchor: raw_anchor_name})
        all_score_cols = [raw_anchor_name if c == anchor else c for c in all_score_cols]
        drop_cols = [c for c in drop_cols if c != anchor]

    # --- Phase 3: Left-join eval onto pred --------------------------------
    if subtable == "pred" or merged_eval is None:
        merged = merged_pred
    else:
        if collector is not None:
            eval_not_in_pred = count_lazy(
                merged_eval.join(merged_pred, on=eval_keys, how="anti")
            )
            eval_total = collector.stages["eval_merge"]["rows"]
            collector.record("coverage", "eval_not_in_pred", eval_not_in_pred)
            collector.record("coverage", "eval_in_pred", eval_total - eval_not_in_pred)

        print("  Left-joining eval and pred ...", file=sys.stderr)
        merged = merged_eval.join(
            merged_pred, on=eval_keys, how="left", coalesce=True,
        )

    # --- Per-score null counts (after Phase 3 join) -----------------------
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
                merged.select([
                    pl.col(c).null_count().alias(c) for c in score_cols_to_check
                ]).collect().row(0, named=True)
            )
            for col, n in null_counts.items():
                collector.record("per_score_nulls", col, n)

    # --- Drop raw score columns (if requested) -------------------------------
    has_percentiles = aggregate_genes or percentile_order != "none"
    if not keep_raw_scores and drop_cols and has_percentiles:
        print("  Dropping raw score columns ...", file=sys.stderr)
        merged = merged.drop(drop_cols)

    # --- Apply filter columns (if requested) --------------------------------
    if filter_uris:
        print(
            f"  Applying {len(filter_uris)} filter table(s) ...",
            file=sys.stderr,
        )
        merged = apply_filters(merged, filter_uris, cache_dir, use_cache=use_cache, store_cache=store_cache)

    # --- Percentile thresholds TSV (if requested) ---------------------------
    if percentile_thresholds:
        score_cols = drop_cols if aggregate_genes else all_score_cols
        thresholds_df = merged_pred.select([
            pl.col(c).quantile(q, interpolation="nearest").alias(f"{c}_p{q}")
            for c in score_cols
            for q in percentile_thresholds
        ]).collect()

        rows: list[dict[str, object]] = []
        for c in score_cols:
            row: dict[str, object] = {"score": c}
            for q in percentile_thresholds:
                row[f"p{q}"] = thresholds_df[f"{c}_p{q}"][0]
            rows.append(row)
        result = pl.DataFrame(rows)

        tsv_path = re.sub(r"\.parquet$", "", output_uri) + ".percentile_thresholds.tsv"
        result.write_csv(tsv_path, separator="\t")
        print(f"  Percentile thresholds written to {tsv_path}", file=sys.stderr)

    # --- Write output -----------------------------------------------------
    print(f"  Writing output to {output_uri} ...", file=sys.stderr)
    merged.sink_parquet(
        output_uri,
        compression="zstd"
    )

    if collector is not None and row_counts_output is not None:
        try:
            collector.record("output", "rows", count_parquet_rows(output_uri))
        except Exception:
            pass
        write_report(collector, row_counts_output)

    elapsed = time.perf_counter() - start
    print(f"Done in {elapsed:.1f}s.", file=sys.stderr)


def main() -> None:
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
            rc_path = os.path.join(rc_path, f"{_derive_stem(args.output)}_counts.md")
        row_counts_output = rc_path
      
    pct_thresholds = None
    if args.percentile_thresholds is not None:
        pct_thresholds = [
            float(v.strip()) for v in args.percentile_thresholds.split(",") if v.strip()
        ]

    run_pipeline(
        prediction_uris=_parse_uri_list(args.prediction_tables) if args.prediction_tables else [],
        evaluation_uris=_parse_uri_list(args.evaluation_tables) if args.evaluation_tables else [],
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
        aggregate_genes=args.aggregate_genes,
        collapse_genes=args.collapse_genes,
        use_cache=args.use_cache,
        store_cache=args.store_cache,
        row_counts_output=row_counts_output,
        percentile_thresholds=pct_thresholds,
        subtable=args.subtable,
        precomputed_prediction=args.precomputed_prediction,
        precomputed_evaluation=args.precomputed_evaluation,
    )


if __name__ == "__main__":
    main()
