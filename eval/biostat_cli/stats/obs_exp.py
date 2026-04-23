"""Observed/expected sum ratio helpers for biostat_cli."""

from __future__ import annotations

import math


def obs_exp_ratio_value(sum_observed: float, sum_expected: float) -> float:
    """
    O/E = sum(observed) / sum(expected) for a set of rows.

    Returns NaN if the denominator is zero, invalid, or non-finite.
    """
    if not math.isfinite(sum_observed) or not math.isfinite(sum_expected) or sum_expected == 0:
        return float("nan")
    return float(sum_observed) / float(sum_expected)
