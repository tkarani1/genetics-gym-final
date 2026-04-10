from __future__ import annotations

import math
from dataclasses import dataclass

from sklearn.metrics import average_precision_score, roc_auc_score


@dataclass(frozen=True)
class PairwiseContinuousResult:
    value: float
    anchor_value: float
    adjustment_ratio: float


def pairwise_continuous_adjust(
    anchor_full_value: float,
    anchor_pairwise_value: float,
    vsm_pairwise_value: float,
) -> PairwiseContinuousResult:
    """Pairwise adjustment: anchor_full * (vsm_pairwise / anchor_pairwise)."""
    if anchor_pairwise_value == 0 or math.isnan(anchor_pairwise_value):
        adjustment_ratio = math.nan
    else:
        adjustment_ratio = vsm_pairwise_value / anchor_pairwise_value

    if math.isnan(anchor_full_value) or math.isnan(adjustment_ratio):
        value = math.nan
    else:
        value = anchor_full_value * adjustment_ratio

    return PairwiseContinuousResult(
        value=value,
        anchor_value=anchor_full_value,
        adjustment_ratio=adjustment_ratio,
    )


def compute_auc(labels: list[int], scores: list[float]) -> float:
    if len(labels) == 0 or len(set(labels)) < 2:
        return math.nan
    return float(roc_auc_score(labels, scores))


def compute_auprc(labels: list[int], scores: list[float]) -> float:
    if len(labels) == 0 or len(set(labels)) < 2:
        return math.nan
    return float(average_precision_score(labels, scores))
