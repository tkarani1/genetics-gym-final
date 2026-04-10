from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import fisher_exact, poisson

from biostat_cli.evaluators.base import Contingency

PVALUE_METHODS = ("fisher", "poisson")
DEFAULT_PVALUE_METHOD = "fisher"


@dataclass(frozen=True)
class BinaryStatResult:
    value: float
    p_value: float


@dataclass(frozen=True)
class PairwiseStatResult:
    """Result for pairwise-adjusted enrichment/rate_ratio calculations."""

    value: float  # Final adjusted value: anchor_value * adjustment_ratio
    p_value: float  # p-value from vsm contingency on pairwise intersection
    anchor_value: float  # enr(VSM*, S* ∩ S_e) - the baseline
    adjustment_ratio: float  # enr(VSM_i, S_i ∩ S* ∩ S_e) / enr(VSM*, S_i ∩ S* ∩ S_e)


def _safe_div(num: float, den: float) -> float:
    if den == 0:
        return math.nan
    return num / den


# ---------------------------------------------------------------------------
# Single-contingency helpers (kept for backward compatibility)
# ---------------------------------------------------------------------------

def enrichment(cont: Contingency, pvalue_method: str = DEFAULT_PVALUE_METHOD) -> BinaryStatResult:
    case_rate = _safe_div(cont.tp, cont.tp + cont.fn)
    ctrl_rate = _safe_div(cont.fp, cont.fp + cont.tn)
    value = _safe_div(case_rate, ctrl_rate) if not math.isnan(case_rate) and not math.isnan(ctrl_rate) else math.nan
    return BinaryStatResult(value=value, p_value=compute_p_value(cont, pvalue_method))


def rate_ratio(
    cont: Contingency, case_total: float | None, ctrl_total: float | None,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> BinaryStatResult:
    if case_total is None or ctrl_total is None:
        return BinaryStatResult(value=math.nan, p_value=compute_p_value(cont, pvalue_method))
    case_rate = _safe_div(cont.tp, case_total)
    ctrl_rate = _safe_div(cont.fp, ctrl_total)
    value = _safe_div(case_rate, ctrl_rate) if not math.isnan(case_rate) and not math.isnan(ctrl_rate) else math.nan
    return BinaryStatResult(value=value, p_value=compute_p_value(cont, pvalue_method))


def poisson_p_value(cont: Contingency) -> float:
    above_pos = cont.tp
    above_total = cont.tp + cont.fp
    below_pos = cont.fn
    below_total = cont.fn + cont.tn
    if above_total <= 0 or below_total <= 0:
        return math.nan

    below_rate = below_pos / below_total
    expected = below_rate * above_total
    return float(poisson.sf(above_pos - 1, expected))


def fisher_p_value(cont: Contingency) -> float:
    """One-sided Fisher's exact test (alternative='greater') on the 2×2 table."""
    above_total = cont.tp + cont.fp
    below_total = cont.fn + cont.tn
    if above_total <= 0 or below_total <= 0:
        return math.nan
    table = [[int(round(cont.tp)), int(round(cont.fp))],
             [int(round(cont.fn)), int(round(cont.tn))]]
    _, p = fisher_exact(table, alternative="greater")
    return float(p)


def compute_p_value(cont: Contingency, method: str = DEFAULT_PVALUE_METHOD) -> float:
    if method == "poisson":
        return poisson_p_value(cont)
    if method == "fisher":
        return fisher_p_value(cont)
    raise ValueError(f"Unknown pvalue_method: {method!r}. Must be one of {PVALUE_METHODS}.")


# ---------------------------------------------------------------------------
# Vectorised batch helpers – operate on a list of Contingency objects at once
# ---------------------------------------------------------------------------

def _conts_to_arrays(conts: list[Contingency]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    tp = np.array([c.tp for c in conts])
    fp = np.array([c.fp for c in conts])
    tn = np.array([c.tn for c in conts])
    fn = np.array([c.fn for c in conts])
    return tp, fp, tn, fn


def poisson_p_values_batch(conts: list[Contingency]) -> np.ndarray:
    """Vectorised version of :func:`poisson_p_value`."""
    tp, fp, tn, fn = _conts_to_arrays(conts)
    above_total = tp + fp
    below_pos = fn
    below_total = fn + tn

    valid = (above_total > 0) & (below_total > 0)
    p_values = np.full(len(conts), np.nan)
    if valid.any():
        bt = below_total[valid]
        with np.errstate(divide="ignore", invalid="ignore"):
            below_rate = below_pos[valid] / bt
        expected = below_rate * above_total[valid]
        p_values[valid] = poisson.sf(tp[valid] - 1, expected)
    return p_values


def fisher_p_values_batch(conts: list[Contingency]) -> np.ndarray:
    """Compute Fisher's exact p-values for a list of contingency tables."""
    return np.array([fisher_p_value(c) for c in conts])


def compute_p_values_batch(
    conts: list[Contingency], method: str = DEFAULT_PVALUE_METHOD,
) -> np.ndarray:
    if method == "poisson":
        return poisson_p_values_batch(conts)
    if method == "fisher":
        return fisher_p_values_batch(conts)
    raise ValueError(f"Unknown pvalue_method: {method!r}. Must be one of {PVALUE_METHODS}.")


def enrichment_batch(
    conts: list[Contingency], pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> list[BinaryStatResult]:
    """Vectorised enrichment over many contingency tables."""
    if not conts:
        return []
    tp, fp, tn, fn = _conts_to_arrays(conts)
    case_denom = tp + fn
    ctrl_denom = fp + tn
    with np.errstate(divide="ignore", invalid="ignore"):
        case_rate = np.where(case_denom == 0, np.nan, tp / case_denom)
        ctrl_rate = np.where(ctrl_denom == 0, np.nan, fp / ctrl_denom)
        values = np.where(
            np.isnan(case_rate) | np.isnan(ctrl_rate) | (ctrl_rate == 0),
            np.nan,
            case_rate / ctrl_rate,
        )
    p_values = compute_p_values_batch(conts, pvalue_method)
    return [BinaryStatResult(value=float(values[i]), p_value=float(p_values[i])) for i in range(len(conts))]


def rate_ratio_batch(
    conts: list[Contingency], case_total: float | None, ctrl_total: float | None,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> list[BinaryStatResult]:
    """Vectorised rate-ratio over many contingency tables."""
    if not conts:
        return []
    p_values = compute_p_values_batch(conts, pvalue_method)
    if case_total is None or ctrl_total is None:
        return [BinaryStatResult(value=math.nan, p_value=float(p_values[i])) for i in range(len(conts))]
    tp, fp, _tn, _fn = _conts_to_arrays(conts)
    with np.errstate(divide="ignore", invalid="ignore"):
        case_rate = np.where(case_total == 0, np.nan, tp / case_total)
        ctrl_rate = np.where(ctrl_total == 0, np.nan, fp / ctrl_total)
        values = np.where(
            np.isnan(case_rate) | np.isnan(ctrl_rate) | (ctrl_rate == 0),
            np.nan,
            case_rate / ctrl_rate,
        )
    return [BinaryStatResult(value=float(values[i]), p_value=float(p_values[i])) for i in range(len(conts))]


# ---------------------------------------------------------------------------
# Pairwise-adjusted statistics
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class VsmComparisonResult:
    """Result of a Fisher exact test comparing two VSMs' contingency tables."""

    odds_ratio: float
    p_greater: float
    p_less: float


def vsm_comparison_fisher(cont_a: Contingency, cont_b: Contingency) -> VsmComparisonResult:
    """Fisher exact test on [[TP_a, TP_b], [FP_a, FP_b]].

    Tests whether the TP/FP odds differ between two VSMs above the same
    score threshold.  Returns odds_ratio (a vs b), one-sided p-values for
    'greater' and 'less' alternatives.
    """
    tp_a = int(round(cont_a.tp))
    tp_b = int(round(cont_b.tp))
    fp_a = int(round(cont_a.fp))
    fp_b = int(round(cont_b.fp))

    if tp_a + fp_a == 0 or tp_b + fp_b == 0:
        return VsmComparisonResult(odds_ratio=math.nan, p_greater=math.nan, p_less=math.nan)

    table = [[tp_a, tp_b], [fp_a, fp_b]]
    odds_ratio, p_greater = fisher_exact(table, alternative="greater")
    _, p_less = fisher_exact(table, alternative="less")

    return VsmComparisonResult(
        odds_ratio=float(odds_ratio),
        p_greater=float(p_greater),
        p_less=float(p_less),
    )


def _compute_enrichment_value(cont: Contingency) -> float:
    """Compute raw enrichment value from contingency table."""
    case_rate = _safe_div(cont.tp, cont.tp + cont.fn)
    ctrl_rate = _safe_div(cont.fp, cont.fp + cont.tn)
    if math.isnan(case_rate) or math.isnan(ctrl_rate):
        return math.nan
    return _safe_div(case_rate, ctrl_rate)


def _compute_rate_ratio_value(cont: Contingency, case_total: float, ctrl_total: float) -> float:
    """Compute raw rate ratio value from contingency table."""
    case_rate = _safe_div(cont.tp, case_total)
    ctrl_rate = _safe_div(cont.fp, ctrl_total)
    if math.isnan(case_rate) or math.isnan(ctrl_rate):
        return math.nan
    return _safe_div(case_rate, ctrl_rate)


def pairwise_enrichment(
    anchor_cont_full: Contingency,
    anchor_cont_pairwise: Contingency,
    vsm_cont_pairwise: Contingency,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> PairwiseStatResult:
    """
    Compute pairwise-adjusted enrichment.

    Formula: enr(VSM_i) = enr(VSM*, S*) × [enr(VSM_i, S_i ∩ S*) / enr(VSM*, S_i ∩ S*)]

    Args:
        anchor_cont_full: Contingency for anchor VSM on full set S* ∩ S_e
        anchor_cont_pairwise: Contingency for anchor VSM on pairwise intersection S_i ∩ S* ∩ S_e
        vsm_cont_pairwise: Contingency for VSM_i on pairwise intersection S_i ∩ S* ∩ S_e
        pvalue_method: "fisher" or "poisson"

    Returns:
        PairwiseStatResult with adjusted value, anchor baseline, and adjustment ratio
    """
    anchor_value = _compute_enrichment_value(anchor_cont_full)
    anchor_pairwise_value = _compute_enrichment_value(anchor_cont_pairwise)
    vsm_pairwise_value = _compute_enrichment_value(vsm_cont_pairwise)

    adjustment_ratio = _safe_div(vsm_pairwise_value, anchor_pairwise_value)
    value = anchor_value * adjustment_ratio if not math.isnan(adjustment_ratio) else math.nan

    return PairwiseStatResult(
        value=value,
        p_value=compute_p_value(vsm_cont_pairwise, pvalue_method),
        anchor_value=anchor_value,
        adjustment_ratio=adjustment_ratio,
    )


def pairwise_rate_ratio(
    anchor_cont_full: Contingency,
    anchor_cont_pairwise: Contingency,
    vsm_cont_pairwise: Contingency,
    case_total: float | None,
    ctrl_total: float | None,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> PairwiseStatResult:
    """
    Compute pairwise-adjusted rate ratio.

    Formula: rr(VSM_i) = rr(VSM*, S*) × [rr(VSM_i, S_i ∩ S*) / rr(VSM*, S_i ∩ S*)]

    Args:
        anchor_cont_full: Contingency for anchor VSM on full set S* ∩ S_e
        anchor_cont_pairwise: Contingency for anchor VSM on pairwise intersection S_i ∩ S* ∩ S_e
        vsm_cont_pairwise: Contingency for VSM_i on pairwise intersection S_i ∩ S* ∩ S_e
        case_total: Total number of cases (N1)
        ctrl_total: Total number of controls (N2)
        pvalue_method: "fisher" or "poisson"

    Returns:
        PairwiseStatResult with adjusted value, anchor baseline, and adjustment ratio
    """
    if case_total is None or ctrl_total is None:
        return PairwiseStatResult(
            value=math.nan,
            p_value=compute_p_value(vsm_cont_pairwise, pvalue_method),
            anchor_value=math.nan,
            adjustment_ratio=math.nan,
        )

    anchor_value = _compute_rate_ratio_value(anchor_cont_full, case_total, ctrl_total)
    anchor_pairwise_value = _compute_rate_ratio_value(anchor_cont_pairwise, case_total, ctrl_total)
    vsm_pairwise_value = _compute_rate_ratio_value(vsm_cont_pairwise, case_total, ctrl_total)

    adjustment_ratio = _safe_div(vsm_pairwise_value, anchor_pairwise_value)
    value = anchor_value * adjustment_ratio if not math.isnan(adjustment_ratio) else math.nan

    return PairwiseStatResult(
        value=value,
        p_value=compute_p_value(vsm_cont_pairwise, pvalue_method),
        anchor_value=anchor_value,
        adjustment_ratio=adjustment_ratio,
    )
