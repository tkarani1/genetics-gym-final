"""BioStat CLI package - Memory-efficient genomic statistics using Polars."""

from biostat_cli.cli import RunArgs, run, run_from_frame
from biostat_cli.config import ALL_STATS, DEFAULT_THRESHOLDS, TableConfig
from biostat_cli.types import (
    EvalLevel,
    MissingMode,
    OutputLayout,
    PanelLayoutConfig,
    PipelineMode,
    RateRatioDenominators,
    StatType,
)

__all__ = [
    "__version__",
    "ALL_STATS",
    "DEFAULT_THRESHOLDS",
    "EvalLevel",
    "MissingMode",
    "OutputLayout",
    "PanelLayoutConfig",
    "PipelineMode",
    "RateRatioDenominators",
    "RunArgs",
    "StatType",
    "TableConfig",
    "run",
    "run_from_frame",
]
__version__ = "0.2.0"
