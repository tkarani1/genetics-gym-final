from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.stats import norm
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


def auc_variance_hanley_mcneil(auc: float, n_positive: int, n_negative: int) -> float:
    """Hanley–McNeil (1982) variance of the trapezoidal AUC estimator."""
    if n_positive < 1 or n_negative < 1 or math.isnan(auc):
        return math.nan
    a = float(auc)
    # Keep AUC in (0,1) so Q1, Q2 are finite; boundaries handled in p-value path.
    eps = 1e-15
    a = min(max(a, eps), 1.0 - eps)
    q1 = a / (2.0 - a)
    q2 = 2.0 * a * a / (1.0 + a)
    num = (
        a * (1.0 - a)
        + (n_positive - 1) * (q1 - a * a)
        + (n_negative - 1) * (q2 - a * a)
    )
    return max(0.0, num / (n_positive * n_negative))


def compute_auc_p_value(labels: list[int], scores: list[float]) -> float:
    """Two-sided p-value for H0: AUC = 0.5 (Hanley–McNeil SE + normal approximation).

    Uses the Hanley–McNeil variance for the AUC point estimate, then
    z = (AUC - 0.5) / SE with a two-sided normal test.
    """
    auc = compute_auc(labels, scores)
    if math.isnan(auc):
        return math.nan
    n_pos = sum(1 for y in labels if int(y) == 1)
    n_neg = sum(1 for y in labels if int(y) == 0)
    if n_pos < 1 or n_neg < 1:
        return math.nan
    var = auc_variance_hanley_mcneil(auc, n_pos, n_neg)
    if math.isnan(var):
        return math.nan
    se = math.sqrt(var)
    if se < 1e-15:
        if abs(auc - 0.5) < 1e-12:
            return math.nan
        return 0.0
    z = (auc - 0.5) / se
    return float(2.0 * norm.sf(abs(z)))


def compute_auprc(labels: list[int], scores: list[float]) -> float:
    if len(labels) == 0 or len(set(labels)) < 2:
        return math.nan
    return float(average_precision_score(labels, scores))
