from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import polars as pl

from biostat_cli.cli import (
    _compute_rows_for_prepared,
    _effective_stats_for_eval,
    _ensure_rate_ratio_denominators,
    _resolve_eval_totals,
    _schema_has_obs_exp,
)
from biostat_cli.config import (
    PAIRWISE_STATS,
    PairwiseColumns,
    detect_pairwise_columns,
    get_table_config,
    load_resources,
    obs_exp_column_names,
    parse_csv_arg,
    parse_eval_totals,
    parse_stats,
    parse_thresholds,
)
from biostat_cli.evaluators.base import BaseEvaluator
from biostat_cli.evaluators.gene import GeneEvaluator, SUM_VARIANTS_SENTINEL
from biostat_cli.evaluators.variant import VariantEvaluator
from biostat_cli.io import scan_table, write_json, write_tsv
from biostat_cli.stats.binary import (
    DEFAULT_PVALUE_METHOD,
    DEFAULT_VSM_COMPARISON_METHOD,
    PVALUE_METHODS,
    VSM_COMPARISON_METHODS,
)

from biostat_cli.utils import WITHIN_GENE_COL, apply_chromosome_filter, parse_chromosomes_arg

ERROR_INVALID_THRESHOLD = 22


@dataclass(frozen=True)
class RunArgs:
    resources_json: str
    table_name: str
    eval_level: str
    stat: str
    eval_set: str | None
    filters: str | None
    thresholds: str | None
    case_total_by_eval: str | None
    ctrl_total_by_eval: str | None
    out_fname: str
    write_missing: str
    within_gene_percentile: bool = False
    pvalue_method: str = DEFAULT_PVALUE_METHOD
    vsm_comparison_method: str = DEFAULT_VSM_COMPARISON_METHOD
    chromosomes: str | None = None


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="BioStat CLI")
    parser.add_argument("--resources-json", default="resources.json")
    parser.add_argument("--table-name", required=True)
    parser.add_argument("--eval-level", required=True, choices=["variant", "gene"])
    parser.add_argument("--stat", default="all")
    parser.add_argument("--eval-set", default=None, help="Comma-separated eval columns")
    parser.add_argument("--filters", default=None, help="Comma-separated logical filter names")
    parser.add_argument("--thresholds", default=None, help="Comma-separated percentile thresholds")
    parser.add_argument(
        "--case-total-by-eval",
        default=None,
        help='Comma-separated per-eval case totals in format "eval_name:value"',
    )
    parser.add_argument(
        "--ctrl-total-by-eval",
        default=None,
        help='Comma-separated per-eval control totals in format "eval_name:value"',
    )
    parser.add_argument(
        "--within-gene-percentile",
        action="store_true",
        default=False,
        help="Transform each score into its within-gene percentile rank (requires 'ensg' column)",
    )
    parser.add_argument(
        "--pvalue-method",
        choices=list(PVALUE_METHODS),
        default=DEFAULT_PVALUE_METHOD,
        help=f"P-value calculation method (default: {DEFAULT_PVALUE_METHOD})",
    )
    parser.add_argument(
        "--vsm-comparison-method",
        choices=list(VSM_COMPARISON_METHODS),
        default=DEFAULT_VSM_COMPARISON_METHOD,
        help="VSM pairwise comparison method for vsm_comparison output (default: fisher)",
    )
    parser.add_argument("--out-fname", required=True)
    parser.add_argument("--write-missing", choices=["none", "all", "any"], default="none")
    parser.add_argument(
        "--chromosomes",
        default=None,
        help="Comma-separated chromosome filter (supports 1-22,X,Y,MT with optional chr prefix)",
    )
    return parser


def _choose_evaluator(eval_level: str, source: pl.LazyFrame) -> BaseEvaluator:
    if eval_level == "variant":
        return VariantEvaluator(source)
    if eval_level == "gene":
        return GeneEvaluator(source)
    raise ValueError(f"Unsupported eval-level: {eval_level}")


def _resolve_eval_cols(raw: str | None, metadata_cols: list[str], eval_level: str) -> list[str]:
    explicit = parse_csv_arg(raw)
    if explicit:
        return explicit
    if metadata_cols:
        return metadata_cols
    if eval_level == "gene":
        # Allow gene weighted mode even when metadata doesn't provide a boolean eval column.
        return [SUM_VARIANTS_SENTINEL]
    raise ValueError("No evaluation columns available from --eval-set or metadata")


def _resolve_filter_cols(raw: str | None, metadata_filters: dict[str, str]) -> list[tuple[str, str | None]]:
    filter_names = parse_csv_arg(raw)
    if not filter_names:
        pairs = list(metadata_filters.items())
    else:
        pairs = []
        for name in filter_names:
            if name.lower() == "none":
                continue
            if name not in metadata_filters:
                raise KeyError(f"Unknown filter name: {name}")
            pairs.append((name, metadata_filters[name]))

    return [("none", None), *pairs]


def _resolve_output_paths(out_fname: str) -> dict[str, str]:
    base = Path(out_fname)
    if base.suffix:
        base = base.with_suffix("")
    prefix = str(base)
    return {
        "tsv": f"{prefix}.tsv",
        "log": f"{prefix}_log.json",
        "missing_tsv": f"{prefix}_missing.tsv",
        "vsm_comparison_tsv": f"{prefix}_vsm_comparison.tsv",
        "curves_tsv": f"{prefix}_curves.tsv",
    }


def _sort_missing_df(df: pl.DataFrame) -> pl.DataFrame:
    cols = set(df.columns)
    group_sort_cols = [col for col in ["eval_name", "filter_name"] if col in cols]
    cat_sort_expr = (
        pl.when(pl.col("missing_category") == "all_methods")
        .then(pl.lit(0))
        .when(pl.col("missing_category") == "partial_methods")
        .then(pl.lit(1))
        .otherwise(pl.lit(2))
        .alias("__cat_sort")
        if "missing_category" in cols
        else pl.lit(2).alias("__cat_sort")
    )
    if "chrom" in cols and "pos" in cols:
        chrom_norm = pl.col("chrom").cast(pl.Utf8).str.replace(r"(?i)^chr", "").str.to_uppercase()
        chr_num = chrom_norm.cast(pl.Int64, strict=False)
        sorted_df = (
            df.with_columns(
                [
                    cat_sort_expr,
                    pl.when(chr_num.is_not_null())
                    .then(chr_num)
                    .when(chrom_norm == "X")
                    .then(pl.lit(23))
                    .when(chrom_norm == "Y")
                    .then(pl.lit(24))
                    .otherwise(pl.lit(10_000))
                    .alias("__chr_sort"),
                    pl.col("pos").cast(pl.Int64, strict=False).fill_null(9_999_999_999).alias("__pos_sort"),
                ]
            )
            .sort([*group_sort_cols, "__cat_sort", "__chr_sort", "__pos_sort"])
            .drop(["__cat_sort", "__chr_sort", "__pos_sort"])
        )
        return sorted_df

    if "CHROM" in cols and "POS" in cols:
        chrom_norm = pl.col("CHROM").cast(pl.Utf8).str.replace(r"(?i)^chr", "").str.to_uppercase()
        chr_num = chrom_norm.cast(pl.Int64, strict=False)
        sorted_df = (
            df.with_columns(
                [
                    cat_sort_expr,
                    pl.when(chr_num.is_not_null())
                    .then(chr_num)
                    .when(chrom_norm == "X")
                    .then(pl.lit(23))
                    .when(chrom_norm == "Y")
                    .then(pl.lit(24))
                    .otherwise(pl.lit(10_000))
                    .alias("__chr_sort"),
                    pl.col("POS").cast(pl.Int64, strict=False).fill_null(9_999_999_999).alias("__pos_sort"),
                ]
            )
            .sort([*group_sort_cols, "__cat_sort", "__chr_sort", "__pos_sort"])
            .drop(["__cat_sort", "__chr_sort", "__pos_sort"])
        )
        return sorted_df

    for gene_col in ["gene_symbol", "GENE_ID", "gene_id", "ensg"]:
        if gene_col in cols:
            return df.with_columns(cat_sort_expr).sort([*group_sort_cols, "__cat_sort", gene_col]).drop("__cat_sort")
    return df


def _resolve_entity_id_cols(frame: pl.LazyFrame) -> list[str]:
    cols = frame.collect_schema().names()
    candidates = [
        ["chrom", "pos", "ref", "alt"],
        ["CHROM", "POS", "REF", "ALT"],
        ["locus", "alleles"],
        ["GENE_ID"],
        ["gene_id"],
        ["ensg"],
        ["gene_symbol"],
    ]
    for candidate in candidates:
        if all(col in cols for col in candidate):
            return candidate
    fallback = [
        c
        for c in ["chrom", "pos", "ref", "alt", "locus", "alleles", "GENE_ID", "gene_id", "ensg", "gene_symbol"]
        if c in cols
    ]
    if fallback:
        return fallback
    raise ValueError("Unable to infer identifier columns for missing-output report.")


def _build_missing_variant_rows(
    prepared_frame: pl.LazyFrame,
    score_cols: list[str],
    eval_name: str,
    filter_name: str,
    mode: str,
) -> list[dict[str, Any]]:
    id_cols = _resolve_entity_id_cols(prepared_frame)
    missing_aliases = [f"__missing__{score}" for score in score_cols]

    lf = prepared_frame.select(
        [pl.col(c) for c in id_cols]
        + [pl.col(score).is_null().alias(alias) for score, alias in zip(score_cols, missing_aliases)]
    )
    grouped = lf.group_by(id_cols).agg([pl.col(alias).any().alias(alias) for alias in missing_aliases])

    missing_name_exprs = [
        pl.when(pl.col(alias)).then(pl.lit(score)).otherwise(pl.lit(None))
        for score, alias in zip(score_cols, missing_aliases)
    ]
    missing_count_expr = pl.sum_horizontal([pl.col(alias).cast(pl.Int64) for alias in missing_aliases]).alias(
        "missing_score_count"
    )
    all_missing_expr = pl.all_horizontal([pl.col(alias) for alias in missing_aliases])
    any_missing_expr = pl.any_horizontal([pl.col(alias) for alias in missing_aliases])

    report = grouped.with_columns(
        [
            pl.lit(eval_name).alias("eval_name"),
            pl.lit(filter_name).alias("filter_name"),
            missing_count_expr,
            pl.concat_str(missing_name_exprs, separator=",", ignore_nulls=True).alias("missing_score_names"),
            pl.when(all_missing_expr)
            .then(pl.lit("all_methods"))
            .otherwise(pl.lit("partial_methods"))
            .alias("missing_category"),
        ]
    )
    if mode == "all":
        report = report.filter(all_missing_expr)
    else:
        report = report.filter(any_missing_expr)

    out_cols = [
        "eval_name",
        "filter_name",
        *id_cols,
        "missing_category",
        "missing_score_count",
        "missing_score_names",
    ]
    return report.select(out_cols).collect(streaming=True).to_dicts()






def _run_eval_filter_combo(
    args: RunArgs,
    source: pl.LazyFrame,
    score_cols: list[str],
    requested_stats: set[str],
    thresholds: list[float],
    eval_col: str,
    filter_name: str,
    filter_col: str | None,
    missing_mode: str | None,
    table_schema: set[str],
    pairwise_cols: PairwiseColumns | None,
    eval_case_total: float | None,
    eval_ctrl_total: float | None,
) -> tuple[list[dict[str, Any]], dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    combo_start = time.perf_counter()
    evaluator = _choose_evaluator(args.eval_level, source)
    obs_exp_columns: tuple[str, str] | None = None
    if _schema_has_obs_exp(table_schema, eval_col):
        if args.eval_level == "gene":
            raise ValueError(
                "Obs/expected evals ({stem}_observed / {stem}_expected) require --eval-level variant."
            )
        obs_exp_columns = obs_exp_column_names(eval_col)
    eff_stats = _effective_stats_for_eval(requested_stats, obs_exp_columns is not None)
    prepared = evaluator.prepare_eval_frame(
        eval_col=eval_col, filter_col=filter_col, obs_exp_columns=obs_exp_columns,
    )
    combo_missing_rows: list[dict[str, Any]] = []
    if missing_mode is not None:
        combo_missing_rows = _build_missing_variant_rows(
            prepared_frame=prepared.frame,
            score_cols=score_cols,
            eval_name=eval_col,
            filter_name=filter_name,
            mode=missing_mode,
        )

    combo_rows, vsm_cmp_rows, curve_rows = _compute_rows_for_prepared(
        evaluator=evaluator,
        prepared=prepared,
        eval_col=eval_col,
        filter_name=filter_name,
        filter_col=filter_col,
        score_cols=score_cols,
        requested_stats=eff_stats,
        thresholds=thresholds,
        eval_case_total=eval_case_total,
        eval_ctrl_total=eval_ctrl_total,
        pairwise_cols=pairwise_cols,
        within_gene_percentile=args.within_gene_percentile,
        pvalue_method=args.pvalue_method,
        vsm_comparison_method=args.vsm_comparison_method,
        gene_col=None,
        obs_exp_columns=obs_exp_columns,
    )

    timing = {
        "eval_name": eval_col,
        "filter_name": filter_name,
        "elapsed_seconds": time.perf_counter() - combo_start,
    }
    return combo_rows, timing, combo_missing_rows, vsm_cmp_rows, curve_rows


def run(
    args: RunArgs,
) -> tuple[pl.DataFrame, list[dict[str, Any]], pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    resources = load_resources(args.resources_json)
    table = get_table_config(resources, args.table_name)
    thresholds = parse_thresholds(args.thresholds)
    requested_stats = parse_stats(args.stat)
    chromosomes = parse_chromosomes_arg(args.chromosomes)

    # Share a single LazyFrame across all workers so parquet metadata is read once.
    source = scan_table(table.path)
    source = apply_chromosome_filter(source, chromosomes)
    table_column_names = source.collect_schema().names()
    table_schema = set(table_column_names)

    need_pairwise = bool(requested_stats & PAIRWISE_STATS)
    pairwise_cols: PairwiseColumns | None = None
    if need_pairwise:
        pairwise_cols = detect_pairwise_columns(table_column_names)
        if pairwise_cols is None:
            raise ValueError(
                "Pairwise stats requested but pairwise column structure not detected. "
                "Expected columns: {anchor}_anchor_percentile, {vsm}_percentile_with_anchor, "
                "{anchor}_anchor_percentile_with_{vsm}"
            )

    if args.within_gene_percentile:
        if args.eval_level == "gene":
            raise ValueError("--within-gene-percentile is not compatible with --eval-level gene.")
        if WITHIN_GENE_COL not in source.collect_schema().names():
            raise ValueError(
                f"--within-gene-percentile requires column '{WITHIN_GENE_COL}' "
                f"but it is not present in {table.path}."
            )

    eval_cols = _resolve_eval_cols(args.eval_set, table.evals, args.eval_level)
    filter_pairs = _resolve_filter_cols(args.filters, table.filters)

    case_totals_by_eval = parse_eval_totals(args.case_total_by_eval, "--case-total-by-eval")
    ctrl_totals_by_eval = parse_eval_totals(args.ctrl_total_by_eval, "--ctrl-total-by-eval")

    totals_by_eval: dict[str, tuple[float | None, float | None]] = {}
    for eval_col in eval_cols:
        eff_stats = _effective_stats_for_eval(requested_stats, _schema_has_obs_exp(table_schema, eval_col))
        eval_case_total, eval_ctrl_total = _resolve_eval_totals(
            eval_col=eval_col,
            table_case_totals=table.case_totals,
            table_ctrl_totals=table.ctrl_totals,
            cli_case_totals=case_totals_by_eval,
            cli_ctrl_totals=ctrl_totals_by_eval,
        )
        _ensure_rate_ratio_denominators(eff_stats, eval_col, eval_case_total, eval_ctrl_total)
        totals_by_eval[eval_col] = (eval_case_total, eval_ctrl_total)

    # Parallel work unit is one eval/filter pair, so multiple eval sets also run concurrently.
    combos = [
        (idx, eval_col, filter_name, filter_col)
        for idx, (eval_col, filter_name, filter_col) in enumerate(
            (ec, fn, fc) for ec in eval_cols for (fn, fc) in filter_pairs
        )
    ]

    rows_by_idx: dict[int, list[dict[str, Any]]] = {}
    timings_by_idx: dict[int, dict[str, Any]] = {}
    missing_by_idx: dict[int, list[dict[str, Any]]] = {}
    vsm_cmp_by_idx: dict[int, list[dict[str, Any]]] = {}
    curve_by_idx: dict[int, list[dict[str, Any]]] = {}
    max_workers = min(len(combos), os.cpu_count() or 1)

    if max_workers <= 1:
        for idx, eval_col, filter_name, filter_col in combos:
            ect, ecf = totals_by_eval[eval_col]
            combo_rows, combo_timing, combo_missing_rows, combo_vsm_cmp, combo_curves = _run_eval_filter_combo(
                args=args,
                source=source,
                score_cols=table.score_cols,
                requested_stats=requested_stats,
                thresholds=thresholds,
                eval_col=eval_col,
                filter_name=filter_name,
                filter_col=filter_col,
                missing_mode=args.write_missing if args.write_missing != "none" else None,
                table_schema=table_schema,
                pairwise_cols=pairwise_cols,
                eval_case_total=ect,
                eval_ctrl_total=ecf,
            )
            rows_by_idx[idx] = combo_rows
            timings_by_idx[idx] = combo_timing
            missing_by_idx[idx] = combo_missing_rows
            vsm_cmp_by_idx[idx] = combo_vsm_cmp
            curve_by_idx[idx] = combo_curves
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(
                    _run_eval_filter_combo,
                    args,
                    source,
                    table.score_cols,
                    requested_stats,
                    thresholds,
                    eval_col,
                    filter_name,
                    filter_col,
                    args.write_missing if args.write_missing != "none" else None,
                    table_schema,
                    pairwise_cols,
                    totals_by_eval[eval_col][0],
                    totals_by_eval[eval_col][1],
                ): idx
                for idx, eval_col, filter_name, filter_col in combos
            }
            for future in as_completed(future_map):
                idx = future_map[future]
                combo_rows, combo_timing, combo_missing_rows, combo_vsm_cmp, combo_curves = future.result()
                rows_by_idx[idx] = combo_rows
                timings_by_idx[idx] = combo_timing
                missing_by_idx[idx] = combo_missing_rows
                vsm_cmp_by_idx[idx] = combo_vsm_cmp
                curve_by_idx[idx] = combo_curves

    rows: list[dict[str, Any]] = []
    for idx in sorted(rows_by_idx.keys()):
        rows.extend(rows_by_idx[idx])
    timings = [timings_by_idx[idx] for idx in sorted(timings_by_idx.keys())]
    missing_rows: list[dict[str, Any]] = []
    for idx in sorted(missing_by_idx.keys()):
        missing_rows.extend(missing_by_idx[idx])
    all_vsm_cmp_rows: list[dict[str, Any]] = []
    for idx in sorted(vsm_cmp_by_idx.keys()):
        all_vsm_cmp_rows.extend(vsm_cmp_by_idx[idx])
    all_curve_rows: list[dict[str, Any]] = []
    for idx in sorted(curve_by_idx.keys()):
        all_curve_rows.extend(curve_by_idx[idx])
    if missing_rows:
        missing_df = _sort_missing_df(pl.DataFrame(missing_rows))
    else:
        missing_df = pl.DataFrame(
            schema={
                "eval_name": pl.String,
                "filter_name": pl.String,
                "missing_category": pl.String,
                "missing_score_count": pl.Int64,
                "missing_score_names": pl.String,
            }
        )
    vsm_comparison_df = pl.DataFrame(all_vsm_cmp_rows) if all_vsm_cmp_rows else pl.DataFrame(
        schema={
            "eval_name": pl.String,
            "filter_name": pl.String,
            "vsm_i": pl.String,
            "vsm_j": pl.String,
            "threshold": pl.Float64,
            "odds_ratio": pl.Float64,
            "p_greater": pl.Float64,
            "p_less": pl.Float64,
            "log_odds_ratio": pl.Float64,
            "standard_error": pl.Float64,
            "log_ci_lower": pl.Float64,
            "log_ci_upper": pl.Float64,
            "conf_interval_lower": pl.Float64,
            "conf_interval_upper": pl.Float64,
            "rows_used_i": pl.Int64,
            "rows_used_j": pl.Int64,
            "rows_used_pair": pl.Int64,
        }
    )
    empty_coverage_df = pl.DataFrame(
        schema={
            "eval_name": pl.String,
            "filter_name": pl.String,
            "score_name": pl.String,
            "gene": pl.String,
            "n_variants_used": pl.Int64,
            "n_variants_excluded": pl.Int64,
            "n_variants_total": pl.Int64,
        }
    )
    curves_df = pl.DataFrame(all_curve_rows) if all_curve_rows else pl.DataFrame(
        schema={
            "eval_name": pl.String,
            "filter_name": pl.String,
            "score_name": pl.String,
            "curve_type": pl.String,
            "point_idx": pl.Int64,
            "score_threshold": pl.Float64,
            "fpr": pl.Float64,
            "tpr": pl.Float64,
            "precision": pl.Float64,
            "recall": pl.Float64,
            "rows_used": pl.Int64,
            "total_eval_rows": pl.Int64,
        }
    )
    return pl.DataFrame(rows), timings, missing_df, vsm_comparison_df, empty_coverage_df, curves_df


def main() -> None:
    parser = _build_parser()
    ns = parser.parse_args()
    args = RunArgs(
        resources_json=ns.resources_json,
        table_name=ns.table_name,
        eval_level=ns.eval_level,
        stat=ns.stat,
        eval_set=ns.eval_set,
        filters=ns.filters,
        thresholds=ns.thresholds,
        case_total_by_eval=ns.case_total_by_eval,
        ctrl_total_by_eval=ns.ctrl_total_by_eval,
        within_gene_percentile=ns.within_gene_percentile,
        out_fname=ns.out_fname,
        write_missing=ns.write_missing,
        pvalue_method=ns.pvalue_method,
        vsm_comparison_method=ns.vsm_comparison_method,
        chromosomes=ns.chromosomes,
    )
    try:
        output_paths = _resolve_output_paths(args.out_fname)
        start = time.perf_counter()
        out, eval_filter_timings, missing_df, vsm_comparison_df, _, curves_df = run(args)
        write_tsv(out, output_paths["tsv"])
        if args.write_missing != "none":
            write_tsv(missing_df, output_paths["missing_tsv"])
        if vsm_comparison_df.height > 0:
            write_tsv(vsm_comparison_df, output_paths["vsm_comparison_tsv"])
        if curves_df.height > 0:
            write_tsv(curves_df, output_paths["curves_tsv"])
        elapsed_seconds = time.perf_counter() - start
        write_json(
            {
                "run_args": asdict(args),
                "table_path": get_table_config(load_resources(args.resources_json), args.table_name).path,
                "output_files": output_paths,
                "elapsed_seconds": elapsed_seconds,
                "eval_filter_elapsed_seconds": eval_filter_timings,
            },
            output_paths["log"],
        )
        print(f"Elapsed time (s): {elapsed_seconds:.3f}")
        for item in eval_filter_timings:
            print(
                f"  eval={item['eval_name']} filter={item['filter_name']} "
                f"time_s={item['elapsed_seconds']:.3f}"
            )
    except ValueError as exc:
        if "threshold" in str(exc).lower():
            print(f"Error [{ERROR_INVALID_THRESHOLD}]: {exc}", file=sys.stderr)
            raise SystemExit(ERROR_INVALID_THRESHOLD) from exc
        raise


if __name__ == "__main__":
    main()
