"""
Merge pipeline for creating consolidated VSM (Variant Scoring Method) tables.

Reads prediction and evaluation tables, merges them on genomic keys
(chrom, pos, ref, alt), optionally computes percentile ranks, and writes
output parquet files.
"""
from __future__ import annotations

from .merge import JOIN_KEYS
from .paths import CACHE_DIR, KNOWN_EXTENSIONS
from .types import LoadedInputs, MergedPrediction, PipelineResult

__all__ = [
    "JOIN_KEYS",
    "CACHE_DIR",
    "KNOWN_EXTENSIONS",
    "LoadedInputs",
    "MergedPrediction",
    "PipelineResult",
]
