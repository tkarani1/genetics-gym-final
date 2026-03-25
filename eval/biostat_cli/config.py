"""Configuration loading and parsing for biostat_cli."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


DEFAULT_THRESHOLDS = [0.90, 0.95, 0.98, 0.99]
ALL_STATS = {"auc", "auprc", "enrichment", "rate_ratio", "pairwise_enrichment", "pairwise_rate_ratio"}
PAIRWISE_STATS = {"pairwise_enrichment", "pairwise_rate_ratio"}

__all__ = [
    "DEFAULT_THRESHOLDS",
    "ALL_STATS",
    "PAIRWISE_STATS",
    "TableConfig",
    "PairwiseColumns",
    "load_resources",
    "get_table_config",
    "parse_csv_arg",
    "parse_thresholds",
    "parse_stats",
    "parse_eval_totals",
    "detect_pairwise_columns",
]


@dataclass(frozen=True)
class TableConfig:
    name: str
    path: str
    level: str
    score_cols: list[str]
    filters: dict[str, str]
    evals: list[str]
    case_totals: dict[str, float]
    ctrl_totals: dict[str, float]


def load_resources(resources_json: str) -> dict[str, Any]:
    """Load resources configuration from JSON file."""
    with Path(resources_json).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def get_table_config(resources: dict[str, Any], table_name: str) -> TableConfig:
    table_info = resources.get("Table_info", {})
    if table_name not in table_info:
        raise KeyError(f"Unknown table-name: {table_name}")

    item = table_info[table_name]
    case_totals_raw = item.get("Case_totals", item.get("case_totals", {}))
    ctrl_totals_raw = item.get("Ctrl_totals", item.get("ctrl_totals", {}))
    return TableConfig(
        name=table_name,
        path=item["Path"],
        level=str(item["Level"]).lower(),
        score_cols=list(item.get("Score_cols", [])),
        filters=dict(item.get("Filters", {})),
        evals=list(item.get("evals", item.get("Evals", []))),
        case_totals={str(k): float(v) for k, v in dict(case_totals_raw).items()},
        ctrl_totals={str(k): float(v) for k, v in dict(ctrl_totals_raw).items()},
    )


def parse_csv_arg(raw: str | None) -> list[str]:
    if raw is None or not raw.strip():
        return []
    return [part.strip() for part in raw.split(",") if part.strip()]


def parse_thresholds(raw: str | None) -> list[float]:
    if raw is None or not raw.strip():
        thresholds = list(DEFAULT_THRESHOLDS)
    else:
        thresholds = [float(part.strip()) for part in raw.split(",") if part.strip()]

    if any(t > 1.0 for t in thresholds):
        raise ValueError(
            f"Invalid threshold(s): {thresholds}. Thresholds are percentile-based fractions in [0, 1] "
            "(e.g., 0.90 for the 90th percentile)."
        )
    return thresholds


def parse_stats(raw: str) -> set[str]:
    value = raw.strip().lower()
    if value == "all":
        return set(ALL_STATS)
    stats = {part.strip().lower() for part in value.split(",") if part.strip()}
    unknown = stats - ALL_STATS
    if unknown:
        raise ValueError(f"Unknown stat(s): {sorted(unknown)}")
    return stats


def parse_eval_totals(raw: str | None, arg_name: str) -> dict[str, float]:
    """
    Parse --*-by-eval CLI argument format: "eval_A:123,eval_B:456".
    """
    if raw is None or not raw.strip():
        return {}

    out: dict[str, float] = {}
    for part in raw.split(","):
        piece = part.strip()
        if not piece:
            continue
        if ":" not in piece:
            raise ValueError(f"Invalid {arg_name} entry '{piece}'. Expected format: eval_name:value")
        eval_name, value_raw = piece.split(":", 1)
        eval_name = eval_name.strip()
        value_raw = value_raw.strip()
        if not eval_name:
            raise ValueError(f"Invalid {arg_name} entry '{piece}'. Eval name cannot be empty.")
        try:
            value = float(value_raw)
        except ValueError as exc:
            raise ValueError(f"Invalid {arg_name} value '{value_raw}' for eval '{eval_name}'.") from exc
        out[eval_name] = value
    return out


@dataclass(frozen=True)
class PairwiseColumns:
    """Detected pairwise column structure for adjusted enrichment/rate_ratio calculations."""

    anchor_base: str  # e.g., "mpc_score"
    anchor_full_col: str  # e.g., "mpc_score_anchor_percentile"
    vsm_pairs: tuple[tuple[str, str, str], ...]  # (vsm_base, vsm_col, anchor_pairwise_col)


def detect_pairwise_columns(columns: list[str]) -> PairwiseColumns | None:
    """
    Detect pairwise column structure from table column names.

    Expected patterns:
    - {anchor}_anchor_percentile: anchor percentile on full set S*
    - {vsm}_percentile_with_anchor: VSM_i percentile on S_i ∩ S*
    - {anchor}_anchor_percentile_with_{vsm}: anchor percentile on S_i ∩ S*

    Returns None if pairwise structure is not detected.
    """
    import re

    # Find anchor column: matches *_anchor_percentile but NOT *_anchor_percentile_with_*
    anchor_pattern = re.compile(r"^(.+)_anchor_percentile$")
    anchor_with_pattern = re.compile(r"^(.+)_anchor_percentile_with_(.+)$")

    anchor_full_col: str | None = None
    anchor_base: str | None = None

    for col in columns:
        # Skip columns that match the "with" pattern
        if anchor_with_pattern.match(col):
            continue
        match = anchor_pattern.match(col)
        if match:
            anchor_full_col = col
            anchor_base = match.group(1)
            break

    if anchor_full_col is None or anchor_base is None:
        return None

    # Find VSM columns: matches *_percentile_with_anchor
    vsm_pattern = re.compile(r"^(.+)_percentile_with_anchor$")
    vsm_cols: dict[str, str] = {}  # vsm_base -> vsm_col

    for col in columns:
        match = vsm_pattern.match(col)
        if match:
            vsm_base = match.group(1)
            vsm_cols[vsm_base] = col

    if not vsm_cols:
        return None

    # Match VSM columns with anchor pairwise columns
    # The anchor pairwise column uses a shortened VSM name (e.g., "esm1b" instead of "esm1b_score")
    vsm_pairs: list[tuple[str, str, str]] = []
    for vsm_base, vsm_col in vsm_cols.items():
        # Try full vsm_base first, then try without common suffixes like "_score"
        vsm_short_names = [vsm_base]
        if vsm_base.endswith("_score"):
            vsm_short_names.append(vsm_base[:-6])  # Remove "_score" suffix

        anchor_pairwise_col = None
        for vsm_short in vsm_short_names:
            candidate = f"{anchor_base}_anchor_percentile_with_{vsm_short}"
            if candidate in columns:
                anchor_pairwise_col = candidate
                break

        if anchor_pairwise_col is not None:
            vsm_pairs.append((vsm_base, vsm_col, anchor_pairwise_col))

    if not vsm_pairs:
        return None

    return PairwiseColumns(
        anchor_base=anchor_base,
        anchor_full_col=anchor_full_col,
        vsm_pairs=tuple(vsm_pairs),
    )
