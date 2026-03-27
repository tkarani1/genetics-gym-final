from __future__ import annotations

import math

from sklearn.metrics import average_precision_score, roc_auc_score


def compute_auc(labels: list[int], scores: list[float]) -> float:
    if len(labels) == 0 or len(set(labels)) < 2:
        return math.nan
    return float(roc_auc_score(labels, scores))


def compute_auprc(labels: list[int], scores: list[float]) -> float:
    if len(labels) == 0 or len(set(labels)) < 2:
        return math.nan
    return float(average_precision_score(labels, scores))
