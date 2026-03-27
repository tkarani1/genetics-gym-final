"""Shared types, enums, and dataclasses for biostat_cli."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Literal


class StatType(str, Enum):
    """Supported statistic types."""

    AUC = "auc"
    AUPRC = "auprc"
    ENRICHMENT = "enrichment"
    RATE_RATIO = "rate_ratio"
    PAIRWISE_ENRICHMENT = "pairwise_enrichment"
    PAIRWISE_RATE_RATIO = "pairwise_rate_ratio"

    @classmethod
    def all(cls) -> set[str]:
        return {s.value for s in cls}

    @classmethod
    def pairwise(cls) -> set[str]:
        return {cls.PAIRWISE_ENRICHMENT.value, cls.PAIRWISE_RATE_RATIO.value}

    @classmethod
    def continuous(cls) -> set[str]:
        return {cls.AUC.value, cls.AUPRC.value}

    @classmethod
    def binary(cls) -> set[str]:
        return {cls.ENRICHMENT.value, cls.RATE_RATIO.value}


class EvalLevel(str, Enum):
    """Evaluation level for statistics."""

    VARIANT = "variant"
    GENE = "gene"


class PipelineMode(str, Enum):
    """Pipeline execution mode."""

    RAW = "raw"
    PAIRWISE = "pairwise"
    BOTH = "both"

    def includes_raw(self) -> bool:
        return self in {PipelineMode.RAW, PipelineMode.BOTH}

    def includes_pairwise(self) -> bool:
        return self in {PipelineMode.PAIRWISE, PipelineMode.BOTH}


class OutputLayout(str, Enum):
    """Output file layout for pipeline results."""

    COMBINED = "combined"
    PER_EVAL = "per_eval"
    BOTH = "both"


class MissingMode(str, Enum):
    """Missing variant report mode."""

    NONE = "none"
    ALL = "all"
    ANY = "any"


# Type aliases
ThresholdList = list[float]
EvalSet = list[str]
FilterPairs = list[tuple[str, str | None]]

# Literal types for strict typing
ProfileType = Literal["paper_figure1", "all_variant"]


@dataclass(frozen=True)
class PanelLayoutConfig:
    """Configuration for panel layout in Figure 1 pipeline."""

    panel_order: list[str]
    panel_eval_map: dict[str, str]
    panel_titles: dict[str, str]
    panel_metrics: dict[str, dict[str, str]]

    def get_eval_for_panel(self, panel_id: str) -> str:
        """Get the eval column name for a panel."""
        return self.panel_eval_map.get(panel_id, panel_id)

    def get_title_for_panel(self, panel_id: str) -> str:
        """Get the display title for a panel."""
        return self.panel_titles.get(panel_id, panel_id)

    def get_stat_for_panel(self, panel_id: str, mode: str) -> str:
        """Get the stat type for a panel in a given mode (raw/pairwise)."""
        return self.panel_metrics.get(panel_id, {}).get(mode, "enrichment")


@dataclass(frozen=True)
class RateRatioDenominators:
    """Per-eval denominators for rate ratio calculation."""

    case_totals: dict[str, float]
    ctrl_totals: dict[str, float]

    def get_totals_for_eval(self, eval_name: str) -> tuple[float | None, float | None]:
        """Get case and control totals for a specific eval."""
        return self.case_totals.get(eval_name), self.ctrl_totals.get(eval_name)

    @classmethod
    def from_dict(cls, data: dict[str, dict[str, float]]) -> RateRatioDenominators:
        """Create from a dict of {eval_name: {case_total: X, ctrl_total: Y}}."""
        case_totals: dict[str, float] = {}
        ctrl_totals: dict[str, float] = {}
        for eval_name, totals in data.items():
            if "case_total" in totals:
                case_totals[eval_name] = totals["case_total"]
            if "ctrl_total" in totals:
                ctrl_totals[eval_name] = totals["ctrl_total"]
        return cls(case_totals=case_totals, ctrl_totals=ctrl_totals)


__all__ = [
    "StatType",
    "EvalLevel",
    "PipelineMode",
    "OutputLayout",
    "MissingMode",
    "ProfileType",
    "ThresholdList",
    "EvalSet",
    "FilterPairs",
    "PanelLayoutConfig",
    "RateRatioDenominators",
]
