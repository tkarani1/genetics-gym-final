"""BioStat CLI package - Memory-efficient genomic statistics using Polars."""

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
    "StatType",
    "EvalLevel",
    "PipelineMode",
    "OutputLayout",
    "MissingMode",
    "PanelLayoutConfig",
    "RateRatioDenominators",
]
__version__ = "0.2.0"
