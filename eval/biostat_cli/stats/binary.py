from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import poisson

from biostat_cli.evaluators.base import Contingency


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

def enrichment(cont: Contingency) -> BinaryStatResult:
    case_rate = _safe_div(cont.tp, cont.tp + cont.fn)
    ctrl_rate = _safe_div(cont.fp, cont.fp + cont.tn)
    value = _safe_div(case_rate, ctrl_rate) if not math.isnan(case_rate) and not math.isnan(ctrl_rate) else math.nan
    return BinaryStatResult(value=value, p_value=poisson_p_value(cont))


def rate_ratio(cont: Contingency, case_total: float | None, ctrl_total: float | None) -> BinaryStatResult:
    if case_total is None or ctrl_total is None:
        return BinaryStatResult(value=math.nan, p_value=poisson_p_value(cont))
    case_rate = _safe_div(cont.tp, case_total)
    ctrl_rate = _safe_div(cont.fp, ctrl_total)
    value = _safe_div(case_rate, ctrl_rate) if not math.isnan(case_rate) and not math.isnan(ctrl_rate) else math.nan
    return BinaryStatResult(value=value, p_value=poisson_p_value(cont))


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


def enrichment_batch(conts: list[Contingency]) -> list[BinaryStatResult]:
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
    p_values = poisson_p_values_batch(conts)
    return [BinaryStatResult(value=float(values[i]), p_value=float(p_values[i])) for i in range(len(conts))]


def rate_ratio_batch(
    conts: list[Contingency], case_total: float | None, ctrl_total: float | None
) -> list[BinaryStatResult]:
    """Vectorised rate-ratio over many contingency tables."""
    if not conts:
        return []
    p_values = poisson_p_values_batch(conts)
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
) -> PairwiseStatResult:
    """
    Compute pairwise-adjusted enrichment.

    Formula: enr(VSM_i) = enr(VSM*, S*) × [enr(VSM_i, S_i ∩ S*) / enr(VSM*, S_i ∩ S*)]

    Args:
        anchor_cont_full: Contingency for anchor VSM on full set S* ∩ S_e
        anchor_cont_pairwise: Contingency for anchor VSM on pairwise intersection S_i ∩ S* ∩ S_e
        vsm_cont_pairwise: Contingency for VSM_i on pairwise intersection S_i ∩ S* ∩ S_e

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
        p_value=poisson_p_value(vsm_cont_pairwise),
        anchor_value=anchor_value,
        adjustment_ratio=adjustment_ratio,
    )


def pairwise_rate_ratio(
    anchor_cont_full: Contingency,
    anchor_cont_pairwise: Contingency,
    vsm_cont_pairwise: Contingency,
    case_total: float | None,
    ctrl_total: float | None,
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

    Returns:
        PairwiseStatResult with adjusted value, anchor baseline, and adjustment ratio
    """
    if case_total is None or ctrl_total is None:
        return PairwiseStatResult(
            value=math.nan,
            p_value=poisson_p_value(vsm_cont_pairwise),
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
        p_value=poisson_p_value(vsm_cont_pairwise),
        anchor_value=anchor_value,
        adjustment_ratio=adjustment_ratio,
    )
