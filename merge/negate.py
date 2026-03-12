"""
Align score columns to a reference scoring method by negating those whose
Pearson correlation with the reference is negative.
"""
from __future__ import annotations

import sys

import polars as pl


def compute_negations(
    lf: pl.LazyFrame,
    score_cols: list[str],
    reference: str,
) -> list[str]:
    """
    Return score column names that have negative Pearson correlation with
    *reference* and therefore need to be negated for directional alignment.

    Materializes a single row of correlation scalars; the reference column
    itself is never included in the result.
    """
    candidates = [c for c in score_cols if c != reference]
    if not candidates:
        return []

    corr_row = lf.select(
        pl.corr(reference, c).alias(c) for c in candidates
    ).collect()

    to_negate: list[str] = []
    for c in candidates:
        r = corr_row[c][0]
        if r is None:
            print(
                f"  WARNING: correlation between '{reference}' and '{c}' "
                f"is null (insufficient overlap); skipping negation for '{c}'",
                file=sys.stderr,
            )
            continue
        direction = "positive" if r >= 0 else "negative"
        print(
            f"  corr({reference}, {c}) = {r:+.4f} ({direction})",
            file=sys.stderr,
        )
        if r < 0:
            to_negate.append(c)

    return to_negate


def negate_scores(
    lf: pl.LazyFrame,
    cols_to_negate: list[str],
) -> pl.LazyFrame:
    """Multiply *cols_to_negate* by -1, keeping the original column names."""
    if not cols_to_negate:
        return lf
    return lf.with_columns((-pl.col(c)).alias(c) for c in cols_to_negate)
