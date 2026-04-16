"""Structured return types for the merge pipeline's phased functions."""
from __future__ import annotations

from dataclasses import dataclass, field

import polars as pl


@dataclass
class LoadedInputs:
    """Output of the input-loading phase."""

    pred_frames: list[pl.LazyFrame]
    all_score_cols: list[str]
    merged_eval: pl.LazyFrame | None
    eval_keys: list[str]
    fields_set: set[str] | None = None
    gene_pred_frames: list[pl.LazyFrame] = field(default_factory=list)
    gene_score_cols: list[str] = field(default_factory=list)


@dataclass
class MergedPrediction:
    """Output of the prediction-merge phase (and subsequent transforms)."""

    frame: pl.LazyFrame
    score_cols: list[str]
    non_anchor_cols: list[str] = field(default_factory=list)
    drop_cols: list[str] = field(default_factory=list)
    negated_cols: list[str] = field(default_factory=list)


@dataclass
class PipelineResult:
    """Metadata about the final output of the pipeline."""

    frame: pl.LazyFrame
    score_cols: list[str]
    eval_cols: list[str]
    join_type: str
    percentile_order: str
    has_gene_aggregation: bool
