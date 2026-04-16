from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import binom, fisher_exact, poisson

from biostat_cli.evaluators.base import Contingency

PVALUE_METHODS = ("fisher", "poisson")
DEFAULT_PVALUE_METHOD = "fisher"
VSM_COMPARISON_METHODS = ("fisher", "poisson")
DEFAULT_VSM_COMPARISON_METHOD = "fisher"

# 95% Wald intervals: exp(log(metric) ± z * SE(log metric)).
_LOG_RATIO_CI_Z = 1.96


@dataclass(frozen=True)
class BinaryStatResult:
    value: float
    p_value: float
    std_error: float = math.nan
    enrichment_ci_lower: float = math.nan
    enrichment_ci_upper: float = math.nan
    rate_ratio_ci_lower: float = math.nan
    rate_ratio_ci_upper: float = math.nan


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


def _enrichment_stderr_ci_from_cells(tp: float, fp: float, fn: float, tn: float) -> tuple[float, float, float]:
    """SE(ln LR+) and 95% CI bounds on LR+ scale for one 2×2 table.

    LR+ = (TP/(TP+FN)) / (FP/(FP+TN));
    SE(ln LR+) = sqrt((1/TP - 1/(TP+FN)) + (1/FP - 1/(FP+TN))).
    Returns (nan, nan, nan) if undefined.
    """
    pos_d = tp + fn
    neg_d = fp + tn
    if tp <= 0 or fp <= 0 or pos_d <= 0 or neg_d <= 0:
        return math.nan, math.nan, math.nan
    term_pos = (1.0 / tp) - (1.0 / pos_d)
    term_neg = (1.0 / fp) - (1.0 / neg_d)
    if term_pos < 0 or term_neg < 0:
        return math.nan, math.nan, math.nan
    rad = term_pos + term_neg
    if rad < 0:
        return math.nan, math.nan, math.nan
    se = math.sqrt(rad)
    case_r = tp / pos_d
    ctrl_r = fp / neg_d
    if ctrl_r <= 0:
        return math.nan, math.nan, math.nan
    lr = case_r / ctrl_r
    if lr <= 0 or math.isnan(lr):
        return math.nan, math.nan, math.nan
    log_lr = math.log(lr)
    margin = _LOG_RATIO_CI_Z * se
    return se, math.exp(log_lr - margin), math.exp(log_lr + margin)


def enrichment_log_lr_stderr_ci(cont: Contingency) -> tuple[float, float, float]:
    """Analytic SE(ln LR+) and 95% CI on LR+; raw cells first, else +0.5 to all four cells."""
    se, lo, hi = _enrichment_stderr_ci_from_cells(cont.tp, cont.fp, cont.fn, cont.tn)
    if not math.isnan(se):
        return se, lo, hi
    return _enrichment_stderr_ci_from_cells(
        cont.tp + 0.5, cont.fp + 0.5, cont.fn + 0.5, cont.tn + 0.5,
    )


# ---------------------------------------------------------------------------
# Single-contingency helpers (kept for backward compatibility)
# ---------------------------------------------------------------------------

def enrichment(cont: Contingency, pvalue_method: str = DEFAULT_PVALUE_METHOD) -> BinaryStatResult:
    case_rate = _safe_div(cont.tp, cont.tp + cont.fn)
    ctrl_rate = _safe_div(cont.fp, cont.fp + cont.tn)
    value = _safe_div(case_rate, ctrl_rate) if not math.isnan(case_rate) and not math.isnan(ctrl_rate) else math.nan
    se, ci_lo, ci_hi = enrichment_log_lr_stderr_ci(cont)
    return BinaryStatResult(
        value=value,
        p_value=compute_p_value(cont, pvalue_method),
        std_error=se,
        enrichment_ci_lower=ci_lo,
        enrichment_ci_upper=ci_hi,
    )


def enrichment_fisher_std_error(cont: Contingency) -> float:
    """Analytic stderr proxy for Fisher enrichment using log(OR) Wald SE.

    Uses Haldane-Anscombe correction and returns the standard error on the
    log-odds-ratio scale:
      SE(log(OR)) = sqrt(1/a + 1/b + 1/c + 1/d)
    where a,b,c,d are TP,FP,FN,TN cells with +0.5 correction.
    """
    above_total = cont.tp + cont.fp
    below_total = cont.fn + cont.tn
    if above_total <= 0 or below_total <= 0:
        return math.nan
    a = float(cont.tp) + 0.5
    b = float(cont.fp) + 0.5
    c = float(cont.fn) + 0.5
    d = float(cont.tn) + 0.5
    return math.sqrt((1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d))


def rate_ratio(
    cont: Contingency, case_total: float | None, ctrl_total: float | None,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> BinaryStatResult:
    if case_total is None or ctrl_total is None:
        return BinaryStatResult(value=math.nan, p_value=compute_p_value(cont, pvalue_method), std_error=math.nan)
    case_rate = _safe_div(cont.tp, case_total)
    ctrl_rate = _safe_div(cont.fp, ctrl_total)
    value = _safe_div(case_rate, ctrl_rate) if not math.isnan(case_rate) and not math.isnan(ctrl_rate) else math.nan
    std_error = rate_ratio_poisson_std_error(cont, value) if pvalue_method == "poisson" else math.nan
    rr_lo, rr_hi = rate_ratio_log_rr_ci_bounds(cont, value)
    return BinaryStatResult(
        value=value,
        p_value=compute_p_value(cont, pvalue_method),
        std_error=std_error,
        rate_ratio_ci_lower=rr_lo,
        rate_ratio_ci_upper=rr_hi,
    )


def rate_ratio_poisson_std_error(cont: Contingency, value: float) -> float:
    """Analytic stderr of RR on RR scale under Poisson approximation.

    Uses SE(log(RR)) = sqrt(1/TP + 1/FP), then converts via RR * SE(log(RR)).
    Returns NaN when TP or FP is non-positive.
    """
    tp = float(cont.tp)
    fp = float(cont.fp)
    if tp <= 0 or fp <= 0 or math.isnan(value):
        return math.nan
    se_log_rr = math.sqrt((1.0 / tp) + (1.0 / fp))
    return value * se_log_rr


def rate_ratio_log_rr_ci_bounds(cont: Contingency, value: float) -> tuple[float, float]:
    """95% Wald CI on the RR scale: exp(log(RR) ± z * SE(log RR)) with SE = sqrt(1/TP + 1/FP)."""
    tp = float(cont.tp)
    fp = float(cont.fp)
    if tp <= 0 or fp <= 0 or math.isnan(value) or value <= 0:
        return math.nan, math.nan
    se_log_rr = math.sqrt((1.0 / tp) + (1.0 / fp))
    log_rr = math.log(value)
    margin = _LOG_RATIO_CI_Z * se_log_rr
    return math.exp(log_rr - margin), math.exp(log_rr + margin)


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
    return [
        BinaryStatResult(
            value=float(values[i]),
            p_value=float(p_values[i]),
            std_error=se,
            enrichment_ci_lower=lo,
            enrichment_ci_upper=hi,
        )
        for i, (se, lo, hi) in enumerate(enrichment_log_lr_stderr_ci(c) for c in conts)
    ]


def rate_ratio_batch(
    conts: list[Contingency], case_total: float | None, ctrl_total: float | None,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> list[BinaryStatResult]:
    """Vectorised rate-ratio over many contingency tables."""
    if not conts:
        return []
    p_values = compute_p_values_batch(conts, pvalue_method)
    if case_total is None or ctrl_total is None:
        return [
            BinaryStatResult(value=math.nan, p_value=float(p_values[i]), std_error=math.nan)
            for i in range(len(conts))
        ]
    tp, fp, _tn, _fn = _conts_to_arrays(conts)
    with np.errstate(divide="ignore", invalid="ignore"):
        case_rate = np.where(case_total == 0, np.nan, tp / case_total)
        ctrl_rate = np.where(ctrl_total == 0, np.nan, fp / ctrl_total)
        values = np.where(
            np.isnan(case_rate) | np.isnan(ctrl_rate) | (ctrl_rate == 0),
            np.nan,
            case_rate / ctrl_rate,
        )
        se_log_rr = np.sqrt((1.0 / tp) + (1.0 / fp))
        poisson_std_errors = np.where((tp > 0) & (fp > 0) & ~np.isnan(values), values * se_log_rr, np.nan)
        valid_ci = (tp > 0) & (fp > 0) & ~np.isnan(values) & (values > 0)
        log_rr = np.log(np.where(valid_ci, values, np.nan))
        margin = _LOG_RATIO_CI_Z * se_log_rr
        rr_ci_lo = np.exp(np.where(valid_ci, log_rr - margin, np.nan))
        rr_ci_hi = np.exp(np.where(valid_ci, log_rr + margin, np.nan))
    std_errors = poisson_std_errors if pvalue_method == "poisson" else np.full(len(conts), np.nan)
    return [
        BinaryStatResult(
            value=float(values[i]),
            p_value=float(p_values[i]),
            std_error=float(std_errors[i]),
            rate_ratio_ci_lower=float(rr_ci_lo[i]),
            rate_ratio_ci_upper=float(rr_ci_hi[i]),
        )
        for i in range(len(conts))
    ]


# ---------------------------------------------------------------------------
# Pairwise-adjusted statistics
# ---------------------------------------------------------------------------


_LOG_OR_CI_Z = 1.96  # 95% two-sided normal multiplier for log(OR) interval


@dataclass(frozen=True)
class VsmComparisonResult:
    """Result of a Fisher exact test comparing two VSMs' contingency tables."""

    odds_ratio: float
    p_greater: float
    p_less: float
    log_odds_ratio: float
    standard_error: float
    log_ci_lower: float
    log_ci_upper: float
    conf_interval_lower: float
    conf_interval_upper: float


def vsm_comparison_fisher(cont_a: Contingency, cont_b: Contingency) -> VsmComparisonResult:
    """Fisher exact test on [[TP_a, TP_b], [FP_a, FP_b]].

    Tests whether the TP/FP odds differ between two VSMs above the same
    score threshold.  Returns odds_ratio (a vs b), one-sided p-values for
    'greater' and 'less' alternatives.

    log_odds_ratio and standard_error use the Haldane–Anscombe +0.5 cell
    correction on the same integer table; SE is for ln(OR).  95% CIs use
    log_OR ± 1.96 * SE, then exp for odds-ratio scale.
    """
    nan_full = VsmComparisonResult(
        odds_ratio=math.nan,
        p_greater=math.nan,
        p_less=math.nan,
        log_odds_ratio=math.nan,
        standard_error=math.nan,
        log_ci_lower=math.nan,
        log_ci_upper=math.nan,
        conf_interval_lower=math.nan,
        conf_interval_upper=math.nan,
    )

    tp_a = int(round(cont_a.tp))
    tp_b = int(round(cont_b.tp))
    fp_a = int(round(cont_a.fp))
    fp_b = int(round(cont_b.fp))

    if tp_a + fp_a == 0 or tp_b + fp_b == 0:
        return nan_full

    table = [[tp_a, tp_b], [fp_a, fp_b]]
    odds_ratio, p_greater = fisher_exact(table, alternative="greater")
    _, p_less = fisher_exact(table, alternative="less")

    a, b, c, d = tp_a + 0.5, tp_b + 0.5, fp_a + 0.5, fp_b + 0.5
    log_odds_ratio = math.log((a * d) / (b * c))
    standard_error = math.sqrt((1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d))
    margin = _LOG_OR_CI_Z * standard_error
    log_ci_lower = log_odds_ratio - margin
    log_ci_upper = log_odds_ratio + margin

    return VsmComparisonResult(
        odds_ratio=float(odds_ratio),
        p_greater=float(p_greater),
        p_less=float(p_less),
        log_odds_ratio=float(log_odds_ratio),
        standard_error=float(standard_error),
        log_ci_lower=float(log_ci_lower),
        log_ci_upper=float(log_ci_upper),
        conf_interval_lower=float(math.exp(log_ci_lower)),
        conf_interval_upper=float(math.exp(log_ci_upper)),
    )


def vsm_comparison_poisson_exact(cont_a: Contingency, cont_b: Contingency) -> VsmComparisonResult:
    """Exact conditional Poisson test using full 2x2 TP/FP counts.

    The test compares model-specific TP rates with FP counts as exposures:
      rate_a = TP_a / FP_a and rate_b = TP_b / FP_b.
    Under H0 (rate_a == rate_b), conditioning on TP_a + TP_b gives:
      TP_a ~ Binomial(TP_a + TP_b, p0 = FP_a / (FP_a + FP_b)).
    """
    nan_full = VsmComparisonResult(
        odds_ratio=math.nan,
        p_greater=math.nan,
        p_less=math.nan,
        log_odds_ratio=math.nan,
        standard_error=math.nan,
        log_ci_lower=math.nan,
        log_ci_upper=math.nan,
        conf_interval_lower=math.nan,
        conf_interval_upper=math.nan,
    )

    tp_a = int(round(cont_a.tp))
    tp_b = int(round(cont_b.tp))
    fp_a = int(round(cont_a.fp))
    fp_b = int(round(cont_b.fp))

    if tp_a + fp_a == 0 or tp_b + fp_b == 0:
        return nan_full

    total_tp = tp_a + tp_b
    total_fp = fp_a + fp_b
    if total_tp <= 0 or total_fp <= 0:
        return nan_full

    p0 = fp_a / total_fp
    p_greater = float(binom.sf(tp_a - 1, total_tp, p0))
    p_less = float(binom.cdf(tp_a, total_tp, p0))

    # Haldane-Anscombe correction for stable log-ratio summaries.
    a, b, c, d = tp_a + 0.5, tp_b + 0.5, fp_a + 0.5, fp_b + 0.5
    ratio = (a * d) / (b * c)
    log_ratio = math.log(ratio)
    standard_error = math.sqrt((1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d))
    margin = _LOG_OR_CI_Z * standard_error
    log_ci_lower = log_ratio - margin
    log_ci_upper = log_ratio + margin

    return VsmComparisonResult(
        odds_ratio=float(ratio),
        p_greater=p_greater,
        p_less=p_less,
        log_odds_ratio=float(log_ratio),
        standard_error=float(standard_error),
        log_ci_lower=float(log_ci_lower),
        log_ci_upper=float(log_ci_upper),
        conf_interval_lower=float(math.exp(log_ci_lower)),
        conf_interval_upper=float(math.exp(log_ci_upper)),
    )


def vsm_comparison(
    cont_a: Contingency,
    cont_b: Contingency,
    method: str = DEFAULT_VSM_COMPARISON_METHOD,
) -> VsmComparisonResult:
    if method == "fisher":
        return vsm_comparison_fisher(cont_a, cont_b)
    if method == "poisson":
        return vsm_comparison_poisson_exact(cont_a, cont_b)
    raise ValueError(
        f"Unknown vsm comparison method: {method!r}. Must be one of {VSM_COMPARISON_METHODS}."
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
