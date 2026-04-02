from __future__ import annotations

import math
from dataclasses import dataclass

from biostat_cli.evaluators.base import Contingency
from biostat_cli.stats.binary import (
    DEFAULT_PVALUE_METHOD,
    enrichment,
    enrichment_batch,
    pairwise_enrichment,
    pairwise_rate_ratio,
    rate_ratio,
    rate_ratio_batch,
)
from biostat_cli.stats.continuous import compute_auc, compute_auprc


@dataclass(frozen=True)
class StatOutput:
    stat: str
    value: float
    p_value: float


@dataclass(frozen=True)
class PairwiseStatOutput:
    """Output for pairwise-adjusted statistics with additional anchor/ratio fields."""

    stat: str
    value: float
    p_value: float
    anchor_value: float
    adjustment_ratio: float


class StatFactory:
    @staticmethod
    def auc(labels: list[int] | None, scores: list[float] | None) -> StatOutput:
        if labels is None or scores is None:
            return StatOutput(stat="auc", value=math.nan, p_value=math.nan)
        return StatOutput(stat="auc", value=compute_auc(labels, scores), p_value=math.nan)

    @staticmethod
    def auprc(labels: list[int] | None, scores: list[float] | None) -> StatOutput:
        if labels is None or scores is None:
            return StatOutput(stat="auprc", value=math.nan, p_value=math.nan)
        return StatOutput(stat="auprc", value=compute_auprc(labels, scores), p_value=math.nan)

    @staticmethod
    def enrichment(cont: Contingency, pvalue_method: str = DEFAULT_PVALUE_METHOD) -> StatOutput:
        out = enrichment(cont, pvalue_method=pvalue_method)
        return StatOutput(stat="enrichment", value=out.value, p_value=out.p_value)

    @staticmethod
    def rate_ratio(
        cont: Contingency, case_total: float | None, ctrl_total: float | None,
        pvalue_method: str = DEFAULT_PVALUE_METHOD,
    ) -> StatOutput:
        out = rate_ratio(cont, case_total=case_total, ctrl_total=ctrl_total, pvalue_method=pvalue_method)
        return StatOutput(stat="rate_ratio", value=out.value, p_value=out.p_value)

    @staticmethod
    def enrichment_batch(
        conts: list[Contingency], pvalue_method: str = DEFAULT_PVALUE_METHOD,
    ) -> list[StatOutput]:
        results = enrichment_batch(conts, pvalue_method=pvalue_method)
        return [StatOutput(stat="enrichment", value=r.value, p_value=r.p_value) for r in results]

    @staticmethod
    def rate_ratio_batch(
        conts: list[Contingency], case_total: float | None, ctrl_total: float | None,
        pvalue_method: str = DEFAULT_PVALUE_METHOD,
    ) -> list[StatOutput]:
        results = rate_ratio_batch(conts, case_total=case_total, ctrl_total=ctrl_total, pvalue_method=pvalue_method)
        return [StatOutput(stat="rate_ratio", value=r.value, p_value=r.p_value) for r in results]

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
