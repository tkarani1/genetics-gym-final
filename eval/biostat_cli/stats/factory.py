from __future__ import annotations

import math
from dataclasses import dataclass

from biostat_cli.evaluators.base import Contingency
from biostat_cli.stats.binary import (
    DEFAULT_PVALUE_METHOD,
    DEFAULT_VSM_COMPARISON_METHOD,
    VsmComparisonResult,
    enrichment,
    enrichment_batch,
    pairwise_enrichment,
    pairwise_rate_ratio,
    rate_ratio,
    rate_ratio_batch,
    vsm_comparison,
)
from biostat_cli.stats.continuous import compute_auc, compute_auc_p_value, compute_auprc, pairwise_continuous_adjust
from biostat_cli.stats.gene_averaged import (
    GeneAvgResult,
    gene_avg_auc,
    gene_avg_auprc,
    gene_avg_enrichment,
    gene_avg_rate_ratio,
)


@dataclass(frozen=True)
class StatOutput:
    stat: str
    value: float
    p_value: float
    std_error: float = math.nan
    enrichment_ci_lower: float = math.nan
    enrichment_ci_upper: float = math.nan


@dataclass(frozen=True)
class PairwiseStatOutput:
    """Output for pairwise-adjusted statistics with additional anchor/ratio fields."""

    stat: str
    value: float
    p_value: float
    anchor_value: float
    adjustment_ratio: float


@dataclass(frozen=True)
class GeneAvgStatOutput:
    stat: str
    value: float
    p_value: float
    std_error: float
    n_genes_used: int
    n_genes_excluded: int


class StatFactory:
    @staticmethod
    def auc(labels: list[int] | None, scores: list[float] | None) -> StatOutput:
        if labels is None or scores is None:
            return StatOutput(stat="auc", value=math.nan, p_value=math.nan, std_error=math.nan)
        return StatOutput(
            stat="auc",
            value=compute_auc(labels, scores),
            p_value=compute_auc_p_value(labels, scores),
            std_error=math.nan,
        )

    @staticmethod
    def auprc(labels: list[int] | None, scores: list[float] | None) -> StatOutput:
        if labels is None or scores is None:
            return StatOutput(stat="auprc", value=math.nan, p_value=math.nan, std_error=math.nan)
        return StatOutput(stat="auprc", value=compute_auprc(labels, scores), p_value=math.nan, std_error=math.nan)

    @staticmethod
    def enrichment(cont: Contingency, pvalue_method: str = DEFAULT_PVALUE_METHOD) -> StatOutput:
        out = enrichment(cont, pvalue_method=pvalue_method)
        return StatOutput(
            stat="enrichment",
            value=out.value,
            p_value=out.p_value,
            std_error=out.std_error,
            enrichment_ci_lower=out.enrichment_ci_lower,
            enrichment_ci_upper=out.enrichment_ci_upper,
        )

    @staticmethod
    def rate_ratio(
        cont: Contingency, case_total: float | None, ctrl_total: float | None,
        pvalue_method: str = DEFAULT_PVALUE_METHOD,
    ) -> StatOutput:
        out = rate_ratio(cont, case_total=case_total, ctrl_total=ctrl_total, pvalue_method=pvalue_method)
        return StatOutput(stat="rate_ratio", value=out.value, p_value=out.p_value, std_error=out.std_error)

    @staticmethod
    def enrichment_batch(
        conts: list[Contingency], pvalue_method: str = DEFAULT_PVALUE_METHOD,
    ) -> list[StatOutput]:
        results = enrichment_batch(conts, pvalue_method=pvalue_method)
        return [
            StatOutput(
                stat="enrichment",
                value=r.value,
                p_value=r.p_value,
                std_error=r.std_error,
                enrichment_ci_lower=r.enrichment_ci_lower,
                enrichment_ci_upper=r.enrichment_ci_upper,
            )
            for r in results
        ]

    @staticmethod
    def rate_ratio_batch(
        conts: list[Contingency], case_total: float | None, ctrl_total: float | None,
        pvalue_method: str = DEFAULT_PVALUE_METHOD,
    ) -> list[StatOutput]:
        results = rate_ratio_batch(conts, case_total=case_total, ctrl_total=ctrl_total, pvalue_method=pvalue_method)
        return [StatOutput(stat="rate_ratio", value=r.value, p_value=r.p_value, std_error=r.std_error) for r in results]

    @staticmethod
    def pairwise_enrichment(
        anchor_cont_full: Contingency,
        anchor_cont_pairwise: Contingency,
        vsm_cont_pairwise: Contingency,
        pvalue_method: str = DEFAULT_PVALUE_METHOD,
    ) -> PairwiseStatOutput:
        out = pairwise_enrichment(
            anchor_cont_full, anchor_cont_pairwise, vsm_cont_pairwise, pvalue_method=pvalue_method,
        )
        return PairwiseStatOutput(
            stat="pairwise_enrichment",
            value=out.value,
            p_value=out.p_value,
            anchor_value=out.anchor_value,
            adjustment_ratio=out.adjustment_ratio,
        )

    @staticmethod
    def pairwise_rate_ratio(
        anchor_cont_full: Contingency,
        anchor_cont_pairwise: Contingency,
        vsm_cont_pairwise: Contingency,
        case_total: float | None,
        ctrl_total: float | None,
        pvalue_method: str = DEFAULT_PVALUE_METHOD,
    ) -> PairwiseStatOutput:
        out = pairwise_rate_ratio(
            anchor_cont_full, anchor_cont_pairwise, vsm_cont_pairwise, case_total, ctrl_total,
            pvalue_method=pvalue_method,
        )
        return PairwiseStatOutput(
            stat="pairwise_rate_ratio",
            value=out.value,
            p_value=out.p_value,
            anchor_value=out.anchor_value,
            adjustment_ratio=out.adjustment_ratio,
        )

    @staticmethod
    def pairwise_auc(
        anchor_full_auc: float,
        anchor_pairwise_auc: float,
        vsm_pairwise_auc: float,
        *,
        delong_p_value: float = math.nan,
    ) -> PairwiseStatOutput:
        out = pairwise_continuous_adjust(anchor_full_auc, anchor_pairwise_auc, vsm_pairwise_auc)
        return PairwiseStatOutput(
            stat="pairwise_auc",
            value=out.value,
            p_value=delong_p_value,
            anchor_value=out.anchor_value,
            adjustment_ratio=out.adjustment_ratio,
        )

    @staticmethod
    def vsm_comparison(
        cont_a: Contingency,
        cont_b: Contingency,
        method: str = DEFAULT_VSM_COMPARISON_METHOD,
    ) -> VsmComparisonResult:
        return vsm_comparison(cont_a, cont_b, method=method)

    @staticmethod
    def gene_avg_enrichment_stat(per_gene_values: list[float], n_total_genes: int) -> GeneAvgStatOutput:
        r = gene_avg_enrichment(per_gene_values, n_total_genes)
        return GeneAvgStatOutput(
            stat="gene_avg_enrichment", value=r.value, p_value=r.p_value,
            std_error=r.std_error, n_genes_used=r.n_genes_used, n_genes_excluded=r.n_genes_excluded,
        )

    @staticmethod
    def gene_avg_rate_ratio_stat(per_gene_values: list[float], n_total_genes: int) -> GeneAvgStatOutput:
        r = gene_avg_rate_ratio(per_gene_values, n_total_genes)
        return GeneAvgStatOutput(
            stat="gene_avg_rate_ratio", value=r.value, p_value=r.p_value,
            std_error=r.std_error, n_genes_used=r.n_genes_used, n_genes_excluded=r.n_genes_excluded,
        )

    @staticmethod
    def gene_avg_auc_stat(per_gene_values: list[float], n_total_genes: int) -> GeneAvgStatOutput:
        r = gene_avg_auc(per_gene_values, n_total_genes)
        return GeneAvgStatOutput(
            stat="gene_avg_auc", value=r.value, p_value=r.p_value,
            std_error=r.std_error, n_genes_used=r.n_genes_used, n_genes_excluded=r.n_genes_excluded,
        )

    @staticmethod
    def gene_avg_auprc_stat(per_gene_values: list[float], n_total_genes: int) -> GeneAvgStatOutput:
        r = gene_avg_auprc(per_gene_values, n_total_genes)
        return GeneAvgStatOutput(
            stat="gene_avg_auprc", value=r.value, p_value=r.p_value,
            std_error=r.std_error, n_genes_used=r.n_genes_used, n_genes_excluded=r.n_genes_excluded,
        )

    @staticmethod
    def pairwise_auprc(
        anchor_full_auprc: float,
        anchor_pairwise_auprc: float,
        vsm_pairwise_auprc: float,
    ) -> PairwiseStatOutput:
        out = pairwise_continuous_adjust(anchor_full_auprc, anchor_pairwise_auprc, vsm_pairwise_auprc)
        return PairwiseStatOutput(
            stat="pairwise_auprc",
            value=out.value,
            p_value=math.nan,
            anchor_value=out.anchor_value,
            adjustment_ratio=out.adjustment_ratio,
        )
