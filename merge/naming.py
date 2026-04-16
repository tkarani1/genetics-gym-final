"""
Column naming conventions for the merge pipeline.

Provides a single source of truth for how derived column names are
constructed (percentiles, pairwise pairs, raw anchor, etc.).  This
module is intended to be used by the merge pipeline itself, by the
eval pipeline's ``detect_pairwise_columns``, and by a future UI layer.
"""
from __future__ import annotations


def percentile_col(score: str) -> str:
    """Name of the percentile column for a given score."""
    return f"{score}_percentile"


def pairwise_anchor_col(anchor: str, score: str) -> str:
    """Name of the pairwise anchor-masked column for a given score pair."""
    return f"{anchor}_pairwise_{score}"


def raw_anchor_col(anchor: str) -> str:
    """Name of the raw (non-percentiled) anchor column retained for add-one."""
    return f"_raw_anchor_{anchor}"


def pairwise_percentile_rename_map(
    anchor: str, non_anchor_cols: list[str]
) -> dict[str, str]:
    """Build the rename dict for pairwise percentile columns (variant-level).

    Maps the default percentile naming to the pairwise naming convention:
      - ``{anchor}_percentile`` -> ``{anchor}_anchor_percentile``
      - ``{score}_percentile`` -> ``{score}_percentile_with_anchor``
      - ``{anchor}_pairwise_{score}_percentile`` -> ``{anchor}_anchor_percentile_with_{score}``
    """
    renames: dict[str, str] = {
        percentile_col(anchor): f"{anchor}_anchor_percentile",
    }
    for c in non_anchor_cols:
        renames[percentile_col(c)] = f"{c}_percentile_with_anchor"
        renames[percentile_col(pairwise_anchor_col(anchor, c))] = (
            f"{anchor}_anchor_percentile_with_{c}"
        )
    return renames


def gene_pairwise_percentile_rename_map(
    anchor: str, non_anchor_cols: list[str]
) -> dict[str, str]:
    """Build the rename dict for pairwise percentile columns (gene-level).

    Gene-level aggregation produces ``{score}_mean`` and ``{score}_max``
    columns, each getting their own percentile.  This maps them to:
      - ``{anchor}_{suffix}_percentile`` -> ``{anchor}_anchor_{suffix}_percentile``
      - ``{score}_{suffix}_percentile`` -> ``{score}_{suffix}_percentile_with_anchor``
      - ``{anchor}_pairwise_{score}_{suffix}_percentile`` -> ``{anchor}_anchor_{suffix}_percentile_with_{score}``
    """
    renames: dict[str, str] = {}
    for suffix in ("mean", "max"):
        renames[f"{anchor}_{suffix}_percentile"] = (
            f"{anchor}_anchor_{suffix}_percentile"
        )
        for c in non_anchor_cols:
            renames[f"{c}_{suffix}_percentile"] = (
                f"{c}_{suffix}_percentile_with_anchor"
            )
            renames[f"{anchor}_pairwise_{c}_{suffix}_percentile"] = (
                f"{anchor}_anchor_{suffix}_percentile_with_{c}"
            )
    return renames
