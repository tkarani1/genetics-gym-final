from __future__ import annotations

import argparse
import math
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

from biostat_cli.config import (
    GENE_AVG_STATS,
    OE_STATS,
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
from biostat_cli.stats.binary import (
    DEFAULT_PVALUE_METHOD,
    DEFAULT_VSM_COMPARISON_METHOD,
    PVALUE_METHODS,
    VSM_COMPARISON_METHODS,
)
from biostat_cli.stats.continuous import compute_auc, compute_auprc, delong_two_auc_p_value
from biostat_cli.utils import WITHIN_GENE_COL, missing_category_sort_expr, normalize_chromosome_sort_expr
from biostat_cli.evaluators.base import BaseEvaluator, Contingency, PreparedFrame, slice_prepared_for_score
from biostat_cli.evaluators.gene import GeneEvaluator, SUM_VARIANTS_SENTINEL
from biostat_cli.evaluators.variant import VariantEvaluator
from biostat_cli.io import scan_table, write_json, write_tsv
from biostat_cli.stats.binary import _compute_enrichment_value, _compute_rate_ratio_value
from biostat_cli.stats.continuous import compute_auc as _raw_auc, compute_auprc as _raw_auprc
from biostat_cli.stats.factory import GeneAvgStatOutput, PairwiseStatOutput, StatFactory
from biostat_cli.stats.obs_exp import obs_exp_ratio_value

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
    case_total: float | None
    ctrl_total: float | None
    case_total_by_eval: str | None
    ctrl_total_by_eval: str | None
    bootstrap_samples: int | None
    out_fname: str
    write_missing: str
    within_gene_percentile: bool = False
    pvalue_method: str = DEFAULT_PVALUE_METHOD
    vsm_comparison_method: str = DEFAULT_VSM_COMPARISON_METHOD
    gene_col: str = "ensg"
    write_gene_variant_coverage: bool = False


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="BioStat CLI")
    parser.add_argument("--resources-json", default="resources.json")
    parser.add_argument("--table-name", required=True)
    parser.add_argument("--eval-level", required=True, choices=["variant", "gene"])
    parser.add_argument("--stat", default="all")
    parser.add_argument("--eval-set", default=None, help="Comma-separated eval columns")
    parser.add_argument("--filters", default=None, help="Comma-separated logical filter names")
    parser.add_argument("--thresholds", default=None, help="Comma-separated percentile thresholds")
    parser.add_argument("--case-total", type=float, default=None)
    parser.add_argument("--ctrl-total", type=float, default=None)
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
        "--bootstrap",
        nargs="?",
        const=100,
        type=int,
        default=None,
        metavar="N",
        help="Enable nonparametric row-bootstrap and optionally set samples, e.g. --bootstrap 50 (default: 100)",
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
    parser.add_argument(
        "--gene-col",
        default="ensg",
        help="Gene identifier column for gene-averaged stats (default: ensg)",
    )
    parser.add_argument(
        "--write-gene-variant-coverage",
        action="store_true",
        default=False,
        help="Write per-gene variant coverage report as a separate TSV",
    )
    parser.add_argument("--out-fname", required=True)
    parser.add_argument("--write-missing", choices=["none", "all", "any"], default="none")
    return parser


def _schema_has_obs_exp(schema_names: set[str], eval_stem: str) -> bool:
    obs_c, exp_c = obs_exp_column_names(eval_stem)
    return obs_c in schema_names and exp_c in schema_names


def _effective_stats_for_eval(requested_stats: set[str], is_obs_exp_eval: bool) -> set[str]:
    """Restrict stats to those valid for a boolean-labeled eval vs. an obs/expected eval."""
    if is_obs_exp_eval:
        eff = requested_stats & OE_STATS
        if not eff:
            raise ValueError(
                f"Table has {OE_STATS!r} columns for this eval stem, but --stat does not request any of "
                f"those. Requested: {sorted(requested_stats)}."
            )
        return eff
    eff = requested_stats - OE_STATS
    if not eff and (requested_stats & OE_STATS):
        raise ValueError(
            f"obs_exp_ratio / pairwise_obs_exp_ratio need {{eval_stem}}_observed and {{eval_stem}}_expected. "
            f"Requested: {sorted(requested_stats)}."
        )
    return eff


def _oe_ratio_tail(
    table: pl.DataFrame, score: str, obs: str, exp: str, threshold: float
) -> float:
    s = table.filter(pl.col(score) >= threshold)
    if s.height == 0:
        return float("nan")
    return obs_exp_ratio_value(
        float(s[obs].sum()), float(s[exp].sum()),
    )


def _oe_ratio_for_thresholds(
    score_frame: Any, score_col: str, obs_col: str, exp_col: str, thresholds: list[float]
) -> list[float]:
    df = score_frame.frame.select([pl.col(score_col), pl.col(obs_col), pl.col(exp_col)]).collect(streaming=True)
    if df.height == 0:
        return [float("nan")] * len(thresholds)
    return [_oe_ratio_tail(df, score_col, obs_col, exp_col, t) for t in thresholds]


def _append_obs_exp_row(
    rows: list[dict[str, Any]],
    eval_name: str,
    filter_name: str,
    score_name: str,
    threshold: float,
    value: float,
    rows_used: int,
    total_eval_rows: int,
) -> None:
    _append_binary_row(
        rows=rows,
        eval_name=eval_name,
        filter_name=filter_name,
        score_name=score_name,
        threshold=threshold,
        stat_name="obs_exp_ratio",
        value=value,
        p_value=float("nan"),
        std_error=float("nan"),
        tp=float("nan"),
        fp=float("nan"),
        tn=float("nan"),
        fn=float("nan"),
        rows_used=rows_used,
        total_eval_rows=total_eval_rows,
    )


def _compute_obs_exp_ratio_stats(
    rows: list[dict[str, Any]],
    score_frame: Any,
    eval_stem: str,
    filter_name: str,
    score_col: str,
    obs_col: str,
    exp_col: str,
    thresholds: list[float],
    total_eval_rows: int,
) -> None:
    if not thresholds:
        return
    ratios = _oe_ratio_for_thresholds(score_frame, score_col, obs_col, exp_col, thresholds)
    for t, v in zip(thresholds, ratios):
        _append_obs_exp_row(
            rows, eval_stem, filter_name, score_col, t, v, score_frame.rows_used, total_eval_rows,
        )


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


def _append_binary_row(
    rows: list[dict[str, Any]],
    eval_name: str,
    filter_name: str,
    score_name: str,
    threshold: float,
    stat_name: str,
    value: float,
    p_value: float,
    std_error: float,
    tp: float,
    fp: float,
    tn: float,
    fn: float,
    rows_used: int,
    total_eval_rows: int,
    *,
    enrichment_ci_lower: float = float("nan"),
    enrichment_ci_upper: float = float("nan"),
    rate_ratio_ci_lower: float = float("nan"),
    rate_ratio_ci_upper: float = float("nan"),
) -> None:
    rows.append(
        {
            "eval_name": eval_name,
            "filter_name": filter_name,
            "score_name": score_name,
            "threshold": threshold,
            "stat": stat_name,
            "value": value,
            "p_value": p_value,
            "std_error": std_error,
            "enrichment_ci_lower": enrichment_ci_lower,
            "enrichment_ci_upper": enrichment_ci_upper,
            "rate_ratio_ci_lower": rate_ratio_ci_lower,
            "rate_ratio_ci_upper": rate_ratio_ci_upper,
            "tp": tp,
            "fp": fp,
            "tn": tn,
            "fn": fn,
            "rows_used": rows_used,
            "total_eval_rows": total_eval_rows,
        }
    )


def _append_pairwise_row(
    rows: list[dict[str, Any]],
    eval_name: str,
    filter_name: str,
    score_name: str,
    threshold: float,
    out: PairwiseStatOutput,
    cont: Contingency,
    rows_used: int,
    total_eval_rows: int,
) -> None:
    rows.append(
        {
            "eval_name": eval_name,
            "filter_name": filter_name,
            "score_name": score_name,
            "threshold": threshold,
            "stat": out.stat,
            "value": out.value,
            "p_value": out.p_value,
            "anchor_value": out.anchor_value,
            "adjustment_ratio": out.adjustment_ratio,
            "enrichment_ci_lower": float("nan"),
            "enrichment_ci_upper": float("nan"),
            "rate_ratio_ci_lower": float("nan"),
            "rate_ratio_ci_upper": float("nan"),
            "tp": cont.tp,
            "fp": cont.fp,
            "tn": cont.tn,
            "fn": cont.fn,
            "rows_used": rows_used,
            "total_eval_rows": total_eval_rows,
        }
    )


def _resolve_eval_totals(
    eval_col: str,
    table_case_totals: dict[str, float],
    table_ctrl_totals: dict[str, float],
    cli_case_totals: dict[str, float],
    cli_ctrl_totals: dict[str, float],
    global_case_total: float | None,
    global_ctrl_total: float | None,
) -> tuple[float | None, float | None]:
    """
    Resolve denominators for a specific eval column.

    Priority (high -> low):
    1) per-eval CLI override
    2) per-eval values from resources JSON
    3) global CLI --case-total/--ctrl-total
    """
    case_total = cli_case_totals.get(eval_col, table_case_totals.get(eval_col, global_case_total))
    ctrl_total = cli_ctrl_totals.get(eval_col, table_ctrl_totals.get(eval_col, global_ctrl_total))
    return case_total, ctrl_total


def _threshold_key(value: float) -> str:
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return f"{value:.12g}"


def _row_identity_key(row: dict[str, Any]) -> tuple[str, str, str, str, str]:
    return (
        str(row["eval_name"]),
        str(row["filter_name"]),
        str(row["score_name"]),
        _threshold_key(float(row["threshold"])),
        str(row["stat"]),
    )


def _compute_std_error(values: list[float]) -> float:
    clean = np.asarray([v for v in values if not (isinstance(v, float) and math.isnan(v))], dtype=float)
    if clean.size < 2:
        return math.nan
    return float(np.std(clean, ddof=1))


def _validate_bootstrap_args(args: RunArgs) -> None:
    if args.bootstrap_samples is not None and args.bootstrap_samples < 2:
        raise ValueError("--bootstrap N must be >= 2 when bootstrap is enabled.")


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
        "gene_variant_coverage_tsv": f"{prefix}_gene_variant_coverage.tsv",
    }


def _sort_missing_df(df: pl.DataFrame) -> pl.DataFrame:
    """Sort missing variant DataFrame by eval, filter, category, and genomic position."""
    cols = set(df.columns)
    group_sort_cols = [col for col in ["eval_name", "filter_name"] if col in cols]

    cat_sort_expr = (
        missing_category_sort_expr().alias("__cat_sort")
        if "missing_category" in cols
        else pl.lit(2).alias("__cat_sort")
    )

    # Try lowercase chrom/pos first
    if "chrom" in cols and "pos" in cols:
        return (
            df.with_columns([
                cat_sort_expr,
                normalize_chromosome_sort_expr("chrom").alias("__chr_sort"),
                pl.col("pos").cast(pl.Int64, strict=False).fill_null(9_999_999_999).alias("__pos_sort"),
            ])
            .sort([*group_sort_cols, "__cat_sort", "__chr_sort", "__pos_sort"])
            .drop(["__cat_sort", "__chr_sort", "__pos_sort"])
        )

    # Try uppercase CHROM/POS
    if "CHROM" in cols and "POS" in cols:
        return (
            df.with_columns([
                cat_sort_expr,
                normalize_chromosome_sort_expr("CHROM").alias("__chr_sort"),
                pl.col("POS").cast(pl.Int64, strict=False).fill_null(9_999_999_999).alias("__pos_sort"),
            ])
            .sort([*group_sort_cols, "__cat_sort", "__chr_sort", "__pos_sort"])
            .drop(["__cat_sort", "__chr_sort", "__pos_sort"])
        )

    # Fallback to gene-level sorting
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


def _compute_continuous_stats(
    rows: list[dict[str, Any]],
    evaluator: BaseEvaluator,
    score_frame: Any,
    eval_col: str,
    filter_name: str,
    score_col: str,
    requested_stats: set[str],
    total_eval_rows: int,
) -> None:
    """Compute AUC and AUPRC statistics."""
    labels_scores = evaluator.labels_and_scores(score_frame, eval_col=eval_col, score_col=score_col)
    labels = labels_scores[0] if labels_scores else None
    scores = labels_scores[1] if labels_scores else None

    for stat_name in ["auc", "auprc"]:
        if stat_name not in requested_stats:
            continue
        out = StatFactory.auc(labels, scores) if stat_name == "auc" else StatFactory.auprc(labels, scores)
        _append_binary_row(
            rows=rows,
            eval_name=eval_col,
            filter_name=filter_name,
            score_name=score_col,
            threshold=float("nan"),
            stat_name=out.stat,
            value=out.value,
            p_value=out.p_value,
            std_error=out.std_error,
            tp=float("nan"),
            fp=float("nan"),
            tn=float("nan"),
            fn=float("nan"),
            rows_used=score_frame.rows_used,
            total_eval_rows=total_eval_rows,
            enrichment_ci_lower=float("nan"),
            enrichment_ci_upper=float("nan"),
            rate_ratio_ci_lower=float("nan"),
            rate_ratio_ci_upper=float("nan"),
        )


def _compute_binary_stats(
    rows: list[dict[str, Any]],
    evaluator: BaseEvaluator,
    score_frame: Any,
    eval_col: str,
    filter_name: str,
    score_col: str,
    requested_stats: set[str],
    thresholds: list[float],
    eval_case_total: float | None,
    eval_ctrl_total: float | None,
    total_eval_rows: int,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> list[Contingency]:
    """Compute enrichment and rate_ratio statistics at each threshold.

    Returns the list of Contingency objects (one per threshold) so callers
    can reuse them for cross-VSM comparisons.
    """
    conts = evaluator.contingency_batch(score_frame, eval_col=eval_col, score_col=score_col, thresholds=thresholds)

    if "enrichment" in requested_stats:
        enr_results = StatFactory.enrichment_batch(conts, pvalue_method=pvalue_method)
        for threshold, cont, out in zip(thresholds, conts, enr_results):
            _append_binary_row(
                rows=rows,
                eval_name=eval_col,
                filter_name=filter_name,
                score_name=score_col,
                threshold=threshold,
                stat_name=out.stat,
                value=out.value,
                p_value=out.p_value,
                std_error=out.std_error,
                tp=cont.tp,
                fp=cont.fp,
                tn=cont.tn,
                fn=cont.fn,
                rows_used=score_frame.rows_used,
                total_eval_rows=total_eval_rows,
                enrichment_ci_lower=out.enrichment_ci_lower,
                enrichment_ci_upper=out.enrichment_ci_upper,
                rate_ratio_ci_lower=float("nan"),
                rate_ratio_ci_upper=float("nan"),
            )

    if "rate_ratio" in requested_stats:
        rr_results = StatFactory.rate_ratio_batch(
            conts, case_total=eval_case_total, ctrl_total=eval_ctrl_total, pvalue_method=pvalue_method,
        )
        for threshold, cont, out in zip(thresholds, conts, rr_results):
            _append_binary_row(
                rows=rows,
                eval_name=eval_col,
                filter_name=filter_name,
                score_name=score_col,
                threshold=threshold,
                stat_name=out.stat,
                value=out.value,
                p_value=out.p_value,
                std_error=out.std_error,
                tp=cont.tp,
                fp=cont.fp,
                tn=cont.tn,
                fn=cont.fn,
                rows_used=score_frame.rows_used,
                total_eval_rows=total_eval_rows,
                enrichment_ci_lower=float("nan"),
                enrichment_ci_upper=float("nan"),
                rate_ratio_ci_lower=out.rate_ratio_ci_lower,
                rate_ratio_ci_upper=out.rate_ratio_ci_upper,
            )

    return conts


def _append_gene_avg_row(
    rows: list[dict[str, Any]],
    eval_name: str,
    filter_name: str,
    score_name: str,
    threshold: float,
    out: GeneAvgStatOutput,
    rows_used: int,
    total_eval_rows: int,
) -> None:
    rows.append(
        {
            "eval_name": eval_name,
            "filter_name": filter_name,
            "score_name": score_name,
            "threshold": threshold,
            "stat": out.stat,
            "value": out.value,
            "p_value": out.p_value,
            "std_error": out.std_error,
            "enrichment_ci_lower": float("nan"),
            "enrichment_ci_upper": float("nan"),
            "rate_ratio_ci_lower": float("nan"),
            "rate_ratio_ci_upper": float("nan"),
            "tp": float("nan"),
            "fp": float("nan"),
            "tn": float("nan"),
            "fn": float("nan"),
            "rows_used": rows_used,
            "total_eval_rows": total_eval_rows,
            "n_genes_used": out.n_genes_used,
            "n_genes_excluded": out.n_genes_excluded,
        }
    )


def _compute_gene_averaged_stats(
    rows: list[dict[str, Any]],
    evaluator: BaseEvaluator,
    score_frame: Any,
    eval_col: str,
    filter_name: str,
    score_col: str,
    requested_stats: set[str],
    thresholds: list[float],
    eval_case_total: float | None,
    eval_ctrl_total: float | None,
    total_eval_rows: int,
    gene_col: str,
) -> None:
    need_binary_ga = "gene_avg_enrichment" in requested_stats or "gene_avg_rate_ratio" in requested_stats
    need_continuous_ga = "gene_avg_auc" in requested_stats or "gene_avg_auprc" in requested_stats

    if need_binary_ga and thresholds:
        gene_conts = evaluator.contingency_by_gene_batch(
            score_frame, eval_col=eval_col, score_col=score_col,
            thresholds=thresholds, gene_col=gene_col,
        )
        n_total_genes = len(gene_conts)
        for t_idx, threshold in enumerate(thresholds):
            if "gene_avg_enrichment" in requested_stats:
                per_gene = [_compute_enrichment_value(conts[t_idx]) for conts in gene_conts.values()]
                out = StatFactory.gene_avg_enrichment_stat(per_gene, n_total_genes)
                _append_gene_avg_row(
                    rows, eval_col, filter_name, score_col, threshold, out,
                    score_frame.rows_used, total_eval_rows,
                )
            if "gene_avg_rate_ratio" in requested_stats:
                if eval_case_total is not None and eval_ctrl_total is not None:
                    per_gene = [
                        _compute_rate_ratio_value(conts[t_idx], eval_case_total, eval_ctrl_total)
                        for conts in gene_conts.values()
                    ]
                else:
                    per_gene = [math.nan for _ in gene_conts]
                out = StatFactory.gene_avg_rate_ratio_stat(per_gene, n_total_genes)
                _append_gene_avg_row(
                    rows, eval_col, filter_name, score_col, threshold, out,
                    score_frame.rows_used, total_eval_rows,
                )

    if need_continuous_ga:
        gene_ls = evaluator.labels_and_scores_by_gene(
            score_frame, eval_col=eval_col, score_col=score_col, gene_col=gene_col,
        )
        n_total_genes = len(gene_ls)
        if "gene_avg_auc" in requested_stats:
            per_gene = [_raw_auc(labels, scores) for labels, scores in gene_ls.values()]
            out = StatFactory.gene_avg_auc_stat(per_gene, n_total_genes)
            _append_gene_avg_row(
                rows, eval_col, filter_name, score_col, float("nan"), out,
                score_frame.rows_used, total_eval_rows,
            )
        if "gene_avg_auprc" in requested_stats:
            per_gene = [_raw_auprc(labels, scores) for labels, scores in gene_ls.values()]
            out = StatFactory.gene_avg_auprc_stat(per_gene, n_total_genes)
            _append_gene_avg_row(
                rows, eval_col, filter_name, score_col, float("nan"), out,
                score_frame.rows_used, total_eval_rows,
            )


_NAN_CONTINGENCY = Contingency(tp=float("nan"), fp=float("nan"), tn=float("nan"), fn=float("nan"))


def _compute_pairwise_stats(
    rows: list[dict[str, Any]],
    evaluator: BaseEvaluator,
    prepared: PreparedFrame,
    eval_col: str,
    filter_name: str,
    requested_stats: set[str],
    thresholds: list[float],
    eval_case_total: float | None,
    eval_ctrl_total: float | None,
    pairwise_cols: PairwiseColumns,
    *,
    within_gene_percentile: bool = False,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
    obs_exp_columns: tuple[str, str] | None = None,
) -> None:
    """Compute pairwise statistics (enrichment, rate_ratio, AUC, AUPRC, O/E)."""
    need_pw_binary = bool(requested_stats & {"pairwise_enrichment", "pairwise_rate_ratio"})
    need_pw_continuous = bool(requested_stats & {"pairwise_auc", "pairwise_auprc"})
    need_pairwise_oe = "pairwise_obs_exp_ratio" in requested_stats and obs_exp_columns is not None

    anchor_score_frame = evaluator.prepare_score_frame(
        prepared, score_col=pairwise_cols.anchor_full_col, within_gene_percentile=within_gene_percentile,
    )

    anchor_conts_full: list[Contingency] = []
    if need_pw_binary and thresholds:
        anchor_conts_full = evaluator.contingency_batch(
            anchor_score_frame, eval_col=eval_col, score_col=pairwise_cols.anchor_full_col, thresholds=thresholds
        )

    anchor_full_auc = math.nan
    anchor_full_auprc = math.nan
    if need_pw_continuous:
        anchor_full_ls = evaluator.labels_and_scores(
            anchor_score_frame, eval_col=eval_col, score_col=pairwise_cols.anchor_full_col,
        )
        if anchor_full_ls:
            if "pairwise_auc" in requested_stats:
                anchor_full_auc = compute_auc(anchor_full_ls[0], anchor_full_ls[1])
            if "pairwise_auprc" in requested_stats:
                anchor_full_auprc = compute_auprc(anchor_full_ls[0], anchor_full_ls[1])

    r_full_by_t: list[float] = []
    if need_pairwise_oe and obs_exp_columns is not None and thresholds:
        o_c, e_c = obs_exp_columns
        r_full_by_t = _oe_ratio_for_thresholds(
            anchor_score_frame, pairwise_cols.anchor_full_col, o_c, e_c, thresholds,
        )

    # Process each VSM pair
    for vsm_base, vsm_col, anchor_pairwise_col in pairwise_cols.vsm_pairs:
        pairwise_lf = prepared.frame.filter(pl.col(vsm_col).is_not_null() & pl.col(anchor_pairwise_col).is_not_null())
        pairwise_df = pairwise_lf.collect(streaming=True)
        pairwise_rows_used = pairwise_df.height
        if pairwise_rows_used == 0:
            continue

        pairwise_prepared = PreparedFrame(frame=pairwise_df.lazy(), total_eval_rows=pairwise_rows_used)

        vsm_pw_sf = evaluator.prepare_score_frame(
            pairwise_prepared, score_col=vsm_col, within_gene_percentile=within_gene_percentile,
        )
        anchor_pw_sf = evaluator.prepare_score_frame(
            pairwise_prepared, score_col=anchor_pairwise_col, within_gene_percentile=within_gene_percentile,
        )

        if need_pw_binary and thresholds:
            vsm_conts = evaluator.contingency_batch(
                vsm_pw_sf, eval_col=eval_col, score_col=vsm_col, thresholds=thresholds,
            )
            anchor_conts_pairwise = evaluator.contingency_batch(
                anchor_pw_sf, eval_col=eval_col, score_col=anchor_pairwise_col, thresholds=thresholds,
            )
            for threshold, anchor_cont_full, anchor_cont_pw, vsm_cont in zip(
                thresholds, anchor_conts_full, anchor_conts_pairwise, vsm_conts
            ):
                if "pairwise_enrichment" in requested_stats:
                    out = StatFactory.pairwise_enrichment(
                        anchor_cont_full, anchor_cont_pw, vsm_cont, pvalue_method=pvalue_method,
                    )
                    _append_pairwise_row(
                        rows=rows, eval_name=eval_col, filter_name=filter_name, score_name=vsm_base,
                        threshold=threshold, out=out, cont=vsm_cont,
                        rows_used=pairwise_rows_used, total_eval_rows=prepared.total_eval_rows,
                    )
                if "pairwise_rate_ratio" in requested_stats:
                    out = StatFactory.pairwise_rate_ratio(
                        anchor_cont_full, anchor_cont_pw, vsm_cont, eval_case_total, eval_ctrl_total,
                        pvalue_method=pvalue_method,
                    )
                    _append_pairwise_row(
                        rows=rows, eval_name=eval_col, filter_name=filter_name, score_name=vsm_base,
                        threshold=threshold, out=out, cont=vsm_cont,
                        rows_used=pairwise_rows_used, total_eval_rows=prepared.total_eval_rows,
                    )

        if need_pw_continuous:
            vsm_pw_ls = evaluator.labels_and_scores(vsm_pw_sf, eval_col=eval_col, score_col=vsm_col)
            anchor_pw_ls = evaluator.labels_and_scores(anchor_pw_sf, eval_col=eval_col, score_col=anchor_pairwise_col)

            if "pairwise_auc" in requested_stats:
                vsm_pw_val = compute_auc(vsm_pw_ls[0], vsm_pw_ls[1]) if vsm_pw_ls else math.nan
                anchor_pw_val = compute_auc(anchor_pw_ls[0], anchor_pw_ls[1]) if anchor_pw_ls else math.nan
                p_delong = math.nan
                if anchor_pw_ls and vsm_pw_ls:
                    la, sa = anchor_pw_ls
                    lb, sb = vsm_pw_ls
                    if len(la) == len(lb) and all(int(a) == int(b) for a, b in zip(la, lb)):
                        p_delong = delong_two_auc_p_value(la, sa, sb)
                out = StatFactory.pairwise_auc(
                    anchor_full_auc, anchor_pw_val, vsm_pw_val, delong_p_value=p_delong,
                )
                _append_pairwise_row(
                    rows=rows, eval_name=eval_col, filter_name=filter_name, score_name=vsm_base,
                    threshold=float("nan"), out=out, cont=_NAN_CONTINGENCY,
                    rows_used=pairwise_rows_used, total_eval_rows=prepared.total_eval_rows,
                )

            if "pairwise_auprc" in requested_stats:
                vsm_pw_val = compute_auprc(vsm_pw_ls[0], vsm_pw_ls[1]) if vsm_pw_ls else math.nan
                anchor_pw_val = compute_auprc(anchor_pw_ls[0], anchor_pw_ls[1]) if anchor_pw_ls else math.nan
                out = StatFactory.pairwise_auprc(anchor_full_auprc, anchor_pw_val, vsm_pw_val)
                _append_pairwise_row(
                    rows=rows, eval_name=eval_col, filter_name=filter_name, score_name=vsm_base,
                    threshold=float("nan"), out=out, cont=_NAN_CONTINGENCY,
                    rows_used=pairwise_rows_used, total_eval_rows=prepared.total_eval_rows,
                )

        if need_pairwise_oe and obs_exp_columns is not None and thresholds:
            o_c, e_c = obs_exp_columns
            for t_i, threshold in enumerate(thresholds):
                r_f = r_full_by_t[t_i] if t_i < len(r_full_by_t) else float("nan")
                r_v = _oe_ratio_tail(pairwise_df, vsm_col, o_c, e_c, threshold)
                r_a = _oe_ratio_tail(pairwise_df, anchor_pairwise_col, o_c, e_c, threshold)
                out = StatFactory.pairwise_obs_exp_ratio(r_f, r_v, r_a)
                _append_pairwise_row(
                    rows=rows, eval_name=eval_col, filter_name=filter_name, score_name=vsm_base,
                    threshold=threshold, out=out, cont=_NAN_CONTINGENCY,
                    rows_used=pairwise_rows_used, total_eval_rows=prepared.total_eval_rows,
                )

    # Add anchor baseline rows
    if need_pw_binary and thresholds:
        for threshold, anchor_cont in zip(thresholds, anchor_conts_full):
            if "pairwise_enrichment" in requested_stats:
                out = StatFactory.pairwise_enrichment(
                    anchor_cont, anchor_cont, anchor_cont, pvalue_method=pvalue_method,
                )
                _append_pairwise_row(
                    rows=rows, eval_name=eval_col, filter_name=filter_name, score_name=pairwise_cols.anchor_base,
                    threshold=threshold, out=out, cont=anchor_cont,
                    rows_used=anchor_score_frame.rows_used, total_eval_rows=prepared.total_eval_rows,
                )
            if "pairwise_rate_ratio" in requested_stats:
                out = StatFactory.pairwise_rate_ratio(
                    anchor_cont, anchor_cont, anchor_cont, eval_case_total, eval_ctrl_total,
                    pvalue_method=pvalue_method,
                )
                _append_pairwise_row(
                    rows=rows, eval_name=eval_col, filter_name=filter_name, score_name=pairwise_cols.anchor_base,
                    threshold=threshold, out=out, cont=anchor_cont,
                    rows_used=anchor_score_frame.rows_used, total_eval_rows=prepared.total_eval_rows,
                )

    if need_pw_continuous:
        for stat_name, anchor_full_val, factory_fn in [
            ("pairwise_auc", anchor_full_auc, StatFactory.pairwise_auc),
            ("pairwise_auprc", anchor_full_auprc, StatFactory.pairwise_auprc),
        ]:
            if stat_name not in requested_stats:
                continue
            out = factory_fn(anchor_full_val, anchor_full_val, anchor_full_val)
            _append_pairwise_row(
                rows=rows, eval_name=eval_col, filter_name=filter_name, score_name=pairwise_cols.anchor_base,
                threshold=float("nan"), out=out, cont=_NAN_CONTINGENCY,
                rows_used=anchor_score_frame.rows_used, total_eval_rows=prepared.total_eval_rows,
            )

    if need_pairwise_oe and obs_exp_columns is not None and thresholds:
        for t_i, threshold in enumerate(thresholds):
            r = r_full_by_t[t_i] if t_i < len(r_full_by_t) else float("nan")
            out = StatFactory.pairwise_obs_exp_ratio(r, r, r)
            _append_pairwise_row(
                rows=rows, eval_name=eval_col, filter_name=filter_name, score_name=pairwise_cols.anchor_base,
                threshold=threshold, out=out, cont=_NAN_CONTINGENCY,
                rows_used=anchor_score_frame.rows_used, total_eval_rows=prepared.total_eval_rows,
            )


def _compute_vsm_comparison(
    conts_by_score: dict[str, tuple[list[Contingency], int]],
    eval_col: str,
    filter_name: str,
    thresholds: list[float],
    method: str = DEFAULT_VSM_COMPARISON_METHOD,
) -> list[dict[str, Any]]:
    """All-pairs VSM comparison test using TP/FP counts."""
    score_cols = list(conts_by_score.keys())
    rows: list[dict[str, Any]] = []
    for idx_i in range(len(score_cols)):
        for idx_j in range(idx_i + 1, len(score_cols)):
            col_i, col_j = score_cols[idx_i], score_cols[idx_j]
            conts_i, rows_used_i = conts_by_score[col_i]
            conts_j, rows_used_j = conts_by_score[col_j]
            for t_idx, threshold in enumerate(thresholds):
                result = StatFactory.vsm_comparison(conts_i[t_idx], conts_j[t_idx], method=method)
                rows.append({
                    "eval_name": eval_col,
                    "filter_name": filter_name,
                    "vsm_i": col_i,
                    "vsm_j": col_j,
                    "threshold": threshold,
                    "odds_ratio": result.odds_ratio,
                    "p_greater": result.p_greater,
                    "p_less": result.p_less,
                    "log_odds_ratio": result.log_odds_ratio,
                    "standard_error": result.standard_error,
                    "log_ci_lower": result.log_ci_lower,
                    "log_ci_upper": result.log_ci_upper,
                    "conf_interval_lower": result.conf_interval_lower,
                    "conf_interval_upper": result.conf_interval_upper,
                    "rows_used_i": rows_used_i,
                    "rows_used_j": rows_used_j,
                })
    return rows


def _compute_rows_for_prepared(
    evaluator: BaseEvaluator,
    prepared: PreparedFrame,
    eval_col: str,
    filter_name: str,
    filter_col: str | None,
    score_cols: list[str],
    requested_stats: set[str],
    thresholds: list[float],
    eval_case_total: float | None,
    eval_ctrl_total: float | None,
    pairwise_cols: PairwiseColumns | None,
    *,
    within_gene_percentile: bool = False,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
    vsm_comparison_method: str = DEFAULT_VSM_COMPARISON_METHOD,
    gene_col: str | None = None,
    obs_exp_columns: tuple[str, str] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """
    Compute all requested statistics for a prepared frame.

    Delegates to specialized functions for continuous, binary, pairwise, and
    gene-averaged stats.

    Returns (rows, vsm_comparison_rows).
    """
    rows: list[dict[str, Any]] = []
    is_obs_exp = obs_exp_columns is not None
    need_oe = "obs_exp_ratio" in requested_stats and is_obs_exp
    need_continuous = (not is_obs_exp) and ("auc" in requested_stats or "auprc" in requested_stats)
    need_binary = (not is_obs_exp) and (
        "enrichment" in requested_stats or "rate_ratio" in requested_stats
    )
    need_pairwise = bool(requested_stats & PAIRWISE_STATS)
    need_vsm_comparison = (not is_obs_exp) and "vsm_comparison" in requested_stats
    need_gene_avg = (not is_obs_exp) and bool(requested_stats & GENE_AVG_STATS)

    conts_by_score: dict[str, tuple[list[Contingency], int]] = {}
    prepared_schema = set(prepared.frame.collect_schema().names())

    for score_col in score_cols:
        sub_prepared = slice_prepared_for_score(
            prepared,
            eval_col=eval_col,
            filter_col=filter_col,
            score_col=score_col,
            within_gene_percentile=within_gene_percentile,
            gene_col=gene_col,
            schema_names=prepared_schema,
            obs_exp_columns=obs_exp_columns,
        )
        score_frame = evaluator.prepare_score_frame(
            sub_prepared, score_col=score_col, within_gene_percentile=within_gene_percentile,
        )

        if need_oe and obs_exp_columns is not None and thresholds:
            oc, ec = obs_exp_columns
            _compute_obs_exp_ratio_stats(
                rows, score_frame, eval_col, filter_name, score_col, oc, ec, thresholds, prepared.total_eval_rows,
            )

        if need_continuous:
            _compute_continuous_stats(
                rows, evaluator, score_frame, eval_col, filter_name, score_col,
                requested_stats, prepared.total_eval_rows,
            )

        if (need_binary or need_vsm_comparison) and thresholds:
            conts = _compute_binary_stats(
                rows, evaluator, score_frame, eval_col, filter_name, score_col,
                requested_stats, thresholds, eval_case_total, eval_ctrl_total, prepared.total_eval_rows,
                pvalue_method=pvalue_method,
            )
            if need_vsm_comparison:
                conts_by_score[score_col] = (conts, score_frame.rows_used)

        if need_gene_avg and gene_col is not None:
            _compute_gene_averaged_stats(
                rows, evaluator, score_frame, eval_col, filter_name, score_col,
                requested_stats, thresholds, eval_case_total, eval_ctrl_total,
                prepared.total_eval_rows, gene_col=gene_col,
            )

    if need_pairwise and pairwise_cols is not None:
        _compute_pairwise_stats(
            rows, evaluator, prepared, eval_col, filter_name,
            requested_stats, thresholds, eval_case_total, eval_ctrl_total, pairwise_cols,
            within_gene_percentile=within_gene_percentile,
            pvalue_method=pvalue_method,
            obs_exp_columns=obs_exp_columns,
        )

    vsm_comparison_rows: list[dict[str, Any]] = []
    if need_vsm_comparison and len(conts_by_score) >= 2 and thresholds:
        vsm_comparison_rows = _compute_vsm_comparison(
            conts_by_score, eval_col, filter_name, thresholds, method=vsm_comparison_method,
        )

    return rows, vsm_comparison_rows


def run(args: RunArgs) -> tuple[pl.DataFrame, list[dict[str, Any]], pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    _validate_bootstrap_args(args)
    resources = load_resources(args.resources_json)
    table = get_table_config(resources, args.table_name)
    thresholds = parse_thresholds(args.thresholds)
    requested_stats = parse_stats(args.stat)

    source = scan_table(table.path)
    table_schema = set(source.collect_schema().names())

    if args.within_gene_percentile:
        if args.eval_level == "gene":
            raise ValueError("--within-gene-percentile is not compatible with --eval-level gene.")
        if WITHIN_GENE_COL not in source.collect_schema().names():
            raise ValueError(
                f"--within-gene-percentile requires column '{WITHIN_GENE_COL}' "
                f"but it is not present in {table.path}."
            )

    evaluator = _choose_evaluator(args.eval_level, source)

    eval_cols = _resolve_eval_cols(args.eval_set, table.evals, args.eval_level)
    filter_pairs = _resolve_filter_cols(args.filters, table.filters)
    case_totals_by_eval = parse_eval_totals(args.case_total_by_eval, "--case-total-by-eval")
    ctrl_totals_by_eval = parse_eval_totals(args.ctrl_total_by_eval, "--ctrl-total-by-eval")

    # Validate gene-averaged stats (only when at least one eval is not obs/expected-only)
    need_gene_avg = bool(requested_stats & GENE_AVG_STATS)
    any_non_oe_eval = any(not _schema_has_obs_exp(table_schema, e) for e in eval_cols)
    gene_col: str | None = None
    if need_gene_avg and any_non_oe_eval:
        if args.eval_level == "gene":
            raise ValueError("Gene-averaged stats are not compatible with --eval-level gene.")
        if args.gene_col not in table_schema:
            raise ValueError(
                f"Gene-averaged stats require column '{args.gene_col}' "
                f"but it is not present in {table.path}. "
                f"Use --gene-col to specify the gene identifier column."
            )
        gene_col = args.gene_col

    # Detect pairwise columns if pairwise stats are requested
    need_pairwise = bool(requested_stats & PAIRWISE_STATS)
    pairwise_cols: PairwiseColumns | None = None
    if need_pairwise:
        table_columns = source.collect_schema().names()
        pairwise_cols = detect_pairwise_columns(table_columns)
        if pairwise_cols is None:
            raise ValueError(
                "Pairwise stats requested but pairwise column structure not detected. "
                "Expected columns: {anchor}_anchor_percentile, {vsm}_percentile_with_anchor, "
                "{anchor}_anchor_percentile_with_{vsm}"
            )

    rows: list[dict[str, Any]] = []
    all_vsm_comparison_rows: list[dict[str, Any]] = []
    eval_filter_timings: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []
    coverage_rows: list[dict[str, Any]] = []

    for eval_col in eval_cols:
        eff_stats = _effective_stats_for_eval(
            requested_stats, _schema_has_obs_exp(table_schema, eval_col),
        )
        obs_exp_columns: tuple[str, str] | None = None
        if _schema_has_obs_exp(table_schema, eval_col):
            if args.eval_level == "gene":
                raise ValueError(
                    "Obs/expected evals require --eval-level variant (use <eval_stem>_observed and <eval_stem>_expected columns)."
                )
            obs_exp_columns = obs_exp_column_names(eval_col)
        eval_case_total, eval_ctrl_total = _resolve_eval_totals(
            eval_col=eval_col,
            table_case_totals=table.case_totals,
            table_ctrl_totals=table.ctrl_totals,
            cli_case_totals=case_totals_by_eval,
            cli_ctrl_totals=ctrl_totals_by_eval,
            global_case_total=args.case_total,
            global_ctrl_total=args.ctrl_total,
        )
        for filter_name, filter_col in filter_pairs:
            combo_start = time.perf_counter()
            prepared = evaluator.prepare_eval_frame(
                eval_col=eval_col, filter_col=filter_col, obs_exp_columns=obs_exp_columns,
            )
            if args.write_missing != "none":
                missing_rows.extend(
                    _build_missing_variant_rows(
                        prepared_frame=prepared.frame,
                        score_cols=table.score_cols,
                        eval_name=eval_col,
                        filter_name=filter_name,
                        mode=args.write_missing,
                    )
                )
            if args.write_gene_variant_coverage and gene_col is not None:
                for sc in table.score_cols:
                    cov = (
                        prepared.frame
                        .group_by(gene_col)
                        .agg([
                            pl.col(sc).is_not_null().sum().cast(pl.Int64).alias("n_variants_used"),
                            pl.col(sc).is_null().sum().cast(pl.Int64).alias("n_variants_excluded"),
                            pl.len().cast(pl.Int64).alias("n_variants_total"),
                        ])
                        .collect(streaming=True)
                        .with_columns([
                            pl.lit(eval_col).alias("eval_name"),
                            pl.lit(filter_name).alias("filter_name"),
                            pl.lit(sc).alias("score_name"),
                        ])
                        .rename({gene_col: "gene"})
                        .select(["eval_name", "filter_name", "score_name", "gene",
                                 "n_variants_used", "n_variants_excluded", "n_variants_total"])
                    )
                    coverage_rows.extend(cov.to_dicts())

            combo_rows, vsm_cmp_rows = _compute_rows_for_prepared(
                evaluator=evaluator,
                prepared=prepared,
                eval_col=eval_col,
                filter_name=filter_name,
                filter_col=filter_col,
                score_cols=table.score_cols,
                requested_stats=eff_stats,
                thresholds=thresholds,
                eval_case_total=eval_case_total,
                eval_ctrl_total=eval_ctrl_total,
                pairwise_cols=pairwise_cols,
                within_gene_percentile=args.within_gene_percentile,
                pvalue_method=args.pvalue_method,
                vsm_comparison_method=args.vsm_comparison_method,
                gene_col=gene_col,
                obs_exp_columns=obs_exp_columns,
            )
            all_vsm_comparison_rows.extend(vsm_cmp_rows)

            # Bootstrap std_error is computed from replicate values only; point value/p_value stay from combo_rows.
            if args.bootstrap_samples is not None and combo_rows:
                base_df = prepared.frame.collect(streaming=True)
                n_rows = base_df.height
                bootstrap_values_by_key: dict[tuple[str, str, str, str, str], list[float]] = {}
                if n_rows > 0:
                    rng = np.random.default_rng(42)
                    for _ in range(args.bootstrap_samples):
                        sampled_df = base_df.sample(n=n_rows, with_replacement=True, shuffle=False, seed=int(rng.integers(0, 2**31 - 1)))
                        sample_prepared = PreparedFrame(frame=sampled_df.lazy(), total_eval_rows=sampled_df.height)
                        sample_evaluator = _choose_evaluator(args.eval_level, sample_prepared.frame)
                        sample_rows, _ = _compute_rows_for_prepared(
                            evaluator=sample_evaluator,
                            prepared=sample_prepared,
                            eval_col=eval_col,
                            filter_name=filter_name,
                            filter_col=filter_col,
                            score_cols=table.score_cols,
                            requested_stats=eff_stats,
                            thresholds=thresholds,
                            eval_case_total=eval_case_total,
                            eval_ctrl_total=eval_ctrl_total,
                            pairwise_cols=pairwise_cols,
                            within_gene_percentile=args.within_gene_percentile,
                            pvalue_method=args.pvalue_method,
                            vsm_comparison_method=args.vsm_comparison_method,
                            gene_col=gene_col,
                            obs_exp_columns=obs_exp_columns,
                        )
                        for row in sample_rows:
                            key = _row_identity_key(row)
                            bootstrap_values_by_key.setdefault(key, []).append(float(row["value"]))
                for row in combo_rows:
                    row["std_error"] = _compute_std_error(bootstrap_values_by_key.get(_row_identity_key(row), []))
                    if row.get("stat") == "enrichment":
                        row["enrichment_ci_lower"] = float("nan")
                        row["enrichment_ci_upper"] = float("nan")
                    if row.get("stat") == "rate_ratio":
                        row["rate_ratio_ci_lower"] = float("nan")
                        row["rate_ratio_ci_upper"] = float("nan")
            else:
                for row in combo_rows:
                    row["std_error"] = float(row.get("std_error", math.nan))
            rows.extend(combo_rows)

            eval_filter_timings.append(
                {
                    "eval_name": eval_col,
                    "filter_name": filter_name,
                    "elapsed_seconds": time.perf_counter() - combo_start,
                }
            )
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
    vsm_comparison_df = pl.DataFrame(all_vsm_comparison_rows) if all_vsm_comparison_rows else pl.DataFrame(
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
        }
    )
    coverage_df = pl.DataFrame(coverage_rows) if coverage_rows else pl.DataFrame(
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
    return pl.DataFrame(rows), eval_filter_timings, missing_df, vsm_comparison_df, coverage_df


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
        case_total=ns.case_total,
        ctrl_total=ns.ctrl_total,
        case_total_by_eval=ns.case_total_by_eval,
        ctrl_total_by_eval=ns.ctrl_total_by_eval,
        bootstrap_samples=ns.bootstrap,
        within_gene_percentile=ns.within_gene_percentile,
        out_fname=ns.out_fname,
        write_missing=ns.write_missing,
        pvalue_method=ns.pvalue_method,
        vsm_comparison_method=ns.vsm_comparison_method,
        gene_col=ns.gene_col,
        write_gene_variant_coverage=ns.write_gene_variant_coverage,
    )
    try:
        output_paths = _resolve_output_paths(args.out_fname)
        start = time.perf_counter()
        out, eval_filter_timings, missing_df, vsm_comparison_df, coverage_df = run(args)
        write_tsv(out, output_paths["tsv"])
        if args.write_missing != "none":
            write_tsv(missing_df, output_paths["missing_tsv"])
        if vsm_comparison_df.height > 0:
            write_tsv(vsm_comparison_df, output_paths["vsm_comparison_tsv"])
        if coverage_df.height > 0:
            write_tsv(coverage_df, output_paths["gene_variant_coverage_tsv"])
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
