from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import norm
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve


@dataclass(frozen=True)
class PairwiseContinuousResult:
    value: float
    anchor_value: float
    adjustment_ratio: float


@dataclass(frozen=True)
class ThresholdPointResult:
    tpr: float
    fpr: float
    precision: float
    recall: float
    rows_retained: int
    n_pos_retained: int
    n_neg_retained: int


@dataclass(frozen=True)
class CurvePoint:
    curve_type: str  # "roc" or "pr"
    point_idx: int
    score_threshold: float
    fpr: float
    tpr: float
    precision: float
    recall: float


def _validate_binary(labels: list[int], scores: list[float]) -> bool:
    return len(labels) > 0 and len(set(labels)) >= 2 and len(labels) == len(scores)


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
    if not _validate_binary(labels, scores):
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
    if not _validate_binary(labels, scores):
        return math.nan
    return float(average_precision_score(labels, scores))


def compute_threshold_point_metrics(
    labels: list[int], scores: list[float], threshold: float,
) -> ThresholdPointResult:
    if len(labels) == 0 or len(labels) != len(scores):
        return ThresholdPointResult(
            tpr=math.nan,
            fpr=math.nan,
            precision=math.nan,
            recall=math.nan,
            rows_retained=0,
            n_pos_retained=0,
            n_neg_retained=0,
        )
    preds = [1 if float(s) >= threshold else 0 for s in scores]
    tp = sum(1 for y, p in zip(labels, preds) if int(y) == 1 and p == 1)
    fp = sum(1 for y, p in zip(labels, preds) if int(y) == 0 and p == 1)
    tn = sum(1 for y, p in zip(labels, preds) if int(y) == 0 and p == 0)
    fn = sum(1 for y, p in zip(labels, preds) if int(y) == 1 and p == 0)
    n_pos = tp + fn
    n_neg = fp + tn
    pred_pos = tp + fp
    tpr = float(tp / n_pos) if n_pos > 0 else math.nan
    fpr = float(fp / n_neg) if n_neg > 0 else math.nan
    precision = float(tp / pred_pos) if pred_pos > 0 else math.nan
    return ThresholdPointResult(
        tpr=tpr,
        fpr=fpr,
        precision=precision,
        recall=tpr,
        rows_retained=int(pred_pos),
        n_pos_retained=int(tp),
        n_neg_retained=int(fp),
    )


def _truncate_labels_scores(
    labels: list[int], scores: list[float], threshold: float,
) -> tuple[list[int], list[float]]:
    kept = [(int(y), float(s)) for y, s in zip(labels, scores) if float(s) >= threshold]
    if not kept:
        return [], []
    out_labels = [y for y, _s in kept]
    out_scores = [s for _y, s in kept]
    return out_labels, out_scores


def compute_auc_trunc(labels: list[int], scores: list[float], threshold: float) -> float:
    t_labels, t_scores = _truncate_labels_scores(labels, scores, threshold)
    return compute_auc(t_labels, t_scores)


def compute_auprc_trunc(labels: list[int], scores: list[float], threshold: float) -> float:
    t_labels, t_scores = _truncate_labels_scores(labels, scores, threshold)
    return compute_auprc(t_labels, t_scores)


def truncate_counts(labels: list[int], scores: list[float], threshold: float) -> tuple[int, int, int]:
    t_labels, _t_scores = _truncate_labels_scores(labels, scores, threshold)
    n_total = len(t_labels)
    n_pos = sum(1 for y in t_labels if int(y) == 1)
    n_neg = sum(1 for y in t_labels if int(y) == 0)
    return n_total, n_pos, n_neg


def compute_curve_points(labels: list[int], scores: list[float]) -> list[CurvePoint]:
    if not _validate_binary(labels, scores):
        return []
    out: list[CurvePoint] = []
    roc_fpr, roc_tpr, roc_thresholds = roc_curve(labels, scores)
    for idx, (fpr, tpr, thr) in enumerate(zip(roc_fpr, roc_tpr, roc_thresholds)):
        out.append(
            CurvePoint(
                curve_type="roc",
                point_idx=idx,
                score_threshold=float(thr),
                fpr=float(fpr),
                tpr=float(tpr),
                precision=math.nan,
                recall=math.nan,
            )
        )
    pr_precision, pr_recall, pr_thresholds = precision_recall_curve(labels, scores)
    for idx, (precision, recall) in enumerate(zip(pr_precision, pr_recall)):
        thr = float(pr_thresholds[idx]) if idx < len(pr_thresholds) else math.nan
        out.append(
            CurvePoint(
                curve_type="pr",
                point_idx=idx,
                score_threshold=thr,
                fpr=math.nan,
                tpr=math.nan,
                precision=float(precision),
                recall=float(recall),
            )
        )
    return out


# ---------------------------------------------------------------------------
# Paired DeLong test for two correlated ROC AUCs (Sun & Xu 2014 fast method)
# ---------------------------------------------------------------------------


def _compute_midrank(x: np.ndarray) -> np.ndarray:
    """Midranks for 1D array x (1-based ranks averaged within ties)."""
    x = np.asarray(x, dtype=np.float64)
    j = np.argsort(x)
    z = x[j]
    n = len(x)
    t = np.zeros(n, dtype=np.float64)
    i = 0
    while i < n:
        j_end = i
        while j_end < n and z[j_end] == z[i]:
            j_end += 1
        t[i:j_end] = 0.5 * (i + j_end - 1)
        i = j_end
    t2 = np.empty(n, dtype=np.float64)
    t2[j] = t + 1.0
    return t2


def _cov_kk(v: np.ndarray) -> np.ndarray:
    """Return k×k covariance for v with shape (k, m); zeros if m < 2."""
    k = v.shape[0]
    if v.shape[1] < 2:
        return np.zeros((k, k), dtype=np.float64)
    c = np.cov(v)
    if np.ndim(c) == 0:
        out = np.zeros((k, k), dtype=np.float64)
        out[0, 0] = float(c)
        return out
    return np.asarray(c, dtype=np.float64)


def _fast_delong(
    predictions_sorted_transposed: np.ndarray, label_1_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (aucs, delong_cov) for k classifiers; rows sorted positives first."""
    m = int(label_1_count)
    n = int(predictions_sorted_transposed.shape[1] - m)
    k = predictions_sorted_transposed.shape[0]
    pos = predictions_sorted_transposed[:, :m]
    neg = predictions_sorted_transposed[:, m:]
    tx = np.empty((k, m), dtype=np.float64)
    ty = np.empty((k, n), dtype=np.float64)
    tz = np.empty((k, m + n), dtype=np.float64)
    for r in range(k):
        tx[r, :] = _compute_midrank(pos[r, :])
        ty[r, :] = _compute_midrank(neg[r, :])
        tz[r, :] = _compute_midrank(predictions_sorted_transposed[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - float(m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = _cov_kk(v01)
    sy = _cov_kk(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def _delong_ground_truth_order(labels: np.ndarray) -> tuple[np.ndarray, int]:
    y = np.asarray(labels, dtype=int)
    uniq = np.unique(y)
    if uniq.size != 2 or not np.array_equal(np.sort(uniq), np.array([0, 1])):
        raise ValueError("DeLong AUC comparison requires binary labels in {0, 1}.")
    order = (-y).argsort(kind="mergesort")
    label_1_count = int(y.sum())
    return order, label_1_count


def delong_two_auc_p_value(
    labels: list[int], scores_a: list[float], scores_b: list[float],
) -> float:
    """Two-sided DeLong p-value for H0: AUC(scores_a) = AUC(scores_b) on the same labels.

    Uses Sun & Xu (2014) fast DeLong covariance. Requires identical-length paired
    scores and binary {0, 1} labels.
    """
    if len(labels) != len(scores_a) or len(labels) != len(scores_b):
        return math.nan
    if len(labels) == 0 or len(set(int(y) for y in labels)) < 2:
        return math.nan
    y = np.asarray(labels, dtype=int)
    s1 = np.asarray(scores_a, dtype=np.float64)
    s2 = np.asarray(scores_b, dtype=np.float64)
    try:
        order, m = _delong_ground_truth_order(y)
    except ValueError:
        return math.nan
    preds = np.vstack((s1, s2))[:, order]
    aucs, sigma = _fast_delong(preds, m)
    diff = float(np.abs(aucs[0] - aucs[1]))
    l = np.array([[1.0, -1.0]])
    var_diff = float((l @ sigma @ l.T).squeeze())
    if var_diff <= 0 or math.isnan(var_diff):
        if diff < 1e-15:
            return 1.0
        return math.nan
    z = diff / math.sqrt(var_diff)
    return float(2.0 * norm.sf(z))
