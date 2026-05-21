"""Shared types, enums, and dataclasses for biostat_cli."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Literal


class StatType(str, Enum):
    """Supported statistic types."""

    AUC = "auc"
    AUPRC = "auprc"
    TPR_AT_THRESHOLD = "tpr_at_threshold"
    FPR_AT_THRESHOLD = "fpr_at_threshold"
    PRECISION_AT_THRESHOLD = "precision_at_threshold"
    RECALL_AT_THRESHOLD = "recall_at_threshold"
    AUC_TRUNC = "auc_trunc"
    AUPRC_TRUNC = "auprc_trunc"
    ENRICHMENT = "enrichment"
    RATE_RATIO = "rate_ratio"
    PAIRWISE_ENRICHMENT = "pairwise_enrichment"
    PAIRWISE_RATE_RATIO = "pairwise_rate_ratio"
    PAIRWISE_AUC = "pairwise_auc"
    PAIRWISE_AUPRC = "pairwise_auprc"
    PAIRWISE_TPR_AT_THRESHOLD = "pairwise_tpr_at_threshold"
    PAIRWISE_FPR_AT_THRESHOLD = "pairwise_fpr_at_threshold"
    PAIRWISE_PRECISION_AT_THRESHOLD = "pairwise_precision_at_threshold"
    PAIRWISE_RECALL_AT_THRESHOLD = "pairwise_recall_at_threshold"
    PAIRWISE_AUC_TRUNC = "pairwise_auc_trunc"
    PAIRWISE_AUPRC_TRUNC = "pairwise_auprc_trunc"
    VSM_COMPARISON = "vsm_comparison"
    OBS_EXP_RATIO = "obs_exp_ratio"
    PAIRWISE_OBS_EXP_RATIO = "pairwise_obs_exp_ratio"
    GENE_AVG_ENRICHMENT = "gene_avg_enrichment"
    GENE_AVG_RATE_RATIO = "gene_avg_rate_ratio"
    GENE_AVG_AUC = "gene_avg_auc"
    GENE_AVG_AUPRC = "gene_avg_auprc"
    GENE_AVG_TPR_AT_THRESHOLD = "gene_avg_tpr_at_threshold"
    GENE_AVG_FPR_AT_THRESHOLD = "gene_avg_fpr_at_threshold"
    GENE_AVG_PRECISION_AT_THRESHOLD = "gene_avg_precision_at_threshold"
    GENE_AVG_RECALL_AT_THRESHOLD = "gene_avg_recall_at_threshold"
    GENE_AVG_AUC_TRUNC = "gene_avg_auc_trunc"
    GENE_AVG_AUPRC_TRUNC = "gene_avg_auprc_trunc"

    @classmethod
    def all(cls) -> set[str]:
        return {s.value for s in cls}

    @classmethod
    def pairwise(cls) -> set[str]:
        return {
            cls.PAIRWISE_ENRICHMENT.value, cls.PAIRWISE_RATE_RATIO.value,
            cls.PAIRWISE_AUC.value, cls.PAIRWISE_AUPRC.value,
            cls.PAIRWISE_TPR_AT_THRESHOLD.value, cls.PAIRWISE_FPR_AT_THRESHOLD.value,
            cls.PAIRWISE_PRECISION_AT_THRESHOLD.value, cls.PAIRWISE_RECALL_AT_THRESHOLD.value,
            cls.PAIRWISE_AUC_TRUNC.value, cls.PAIRWISE_AUPRC_TRUNC.value,
            cls.PAIRWISE_OBS_EXP_RATIO.value,
        }

    @classmethod
    def continuous(cls) -> set[str]:
        return {
            cls.AUC.value, cls.AUPRC.value,
            cls.TPR_AT_THRESHOLD.value, cls.FPR_AT_THRESHOLD.value,
            cls.PRECISION_AT_THRESHOLD.value, cls.RECALL_AT_THRESHOLD.value,
            cls.AUC_TRUNC.value, cls.AUPRC_TRUNC.value,
        }

    @classmethod
    def binary(cls) -> set[str]:
        return {cls.ENRICHMENT.value, cls.RATE_RATIO.value}

    @classmethod
    def gene_averaged(cls) -> set[str]:
        return {
            cls.GENE_AVG_ENRICHMENT.value, cls.GENE_AVG_RATE_RATIO.value,
            cls.GENE_AVG_AUC.value, cls.GENE_AVG_AUPRC.value,
            cls.GENE_AVG_TPR_AT_THRESHOLD.value, cls.GENE_AVG_FPR_AT_THRESHOLD.value,
            cls.GENE_AVG_PRECISION_AT_THRESHOLD.value, cls.GENE_AVG_RECALL_AT_THRESHOLD.value,
            cls.GENE_AVG_AUC_TRUNC.value, cls.GENE_AVG_AUPRC_TRUNC.value,
        }


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
