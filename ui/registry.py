"""Typed access to ui/config.json.

Reads the config once and exposes structured dataclasses for the UI,
session, and engine layers. No duplication -- this is a thin accessor
over the JSON structure defined in json_config_design.md.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PredictionEntry:
    path: str
    score_columns: list[str]
    level: str  # "variant" or "gene"
    key_mapping: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class EvaluationEntry:
    path: str
    level: str  # "variant" or "gene"
    key_mapping: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class PipelineConfig:
    join_type: str = "inner"
    percentile_order: str = "post"
    reference_score: str | None = "AM_score"
    anchor: str | None = None
    aggregate_genes: bool = False
    collapse_genes: bool = True
    keep_raw_scores: bool = False
    retain_raw_anchor: bool = False
    percentile_thresholds: list[float] = field(
        default_factory=lambda: [0.85, 0.9, 0.95, 0.995]
    )


@dataclass(frozen=True)
class SmoothingConfig:
    order: str = "none"
    reference_dir: str | None = None
    sigma: float = 10.0


def _parse_prediction(raw: dict[str, Any]) -> PredictionEntry:
    return PredictionEntry(
        path=raw["path"],
        score_columns=list(raw["score_columns"]),
        level=raw.get("level", "variant"),
        key_mapping=dict(raw.get("key_mapping", {})),
    )


def _parse_evaluation(raw: dict[str, Any]) -> EvaluationEntry:
    return EvaluationEntry(
        path=raw["path"],
        level=raw.get("level", "variant"),
        key_mapping=dict(raw.get("key_mapping", {})),
    )


def _parse_pipeline(raw: dict[str, Any] | None) -> PipelineConfig:
    if raw is None:
        return PipelineConfig()
    return PipelineConfig(
        join_type=raw.get("join_type", "inner"),
        percentile_order=raw.get("percentile_order", "post"),
        reference_score=raw.get("reference_score", "AM_score"),
        anchor=raw.get("anchor"),
        aggregate_genes=raw.get("aggregate_genes", False),
        collapse_genes=raw.get("collapse_genes", True),
        keep_raw_scores=raw.get("keep_raw_scores", False),
        retain_raw_anchor=raw.get("retain_raw_anchor", False),
        percentile_thresholds=list(
            raw.get("percentile_thresholds", [0.85, 0.9, 0.95, 0.995])
        ),
    )


def _parse_smoothing(raw: dict[str, Any] | None) -> SmoothingConfig:
    if raw is None:
        return SmoothingConfig()
    return SmoothingConfig(
        order=raw.get("order", "none"),
        reference_dir=raw.get("reference_dir"),
        sigma=float(raw.get("sigma", 10.0)),
    )


class AssetRegistry:
    """Structured access to config.json for the UI, session, and engine layers."""

    def __init__(self, config_path: str = "ui/config.json") -> None:
        path = Path(config_path)
        if not path.exists():
            raise FileNotFoundError(f"Config not found: {path}")
        with path.open("r", encoding="utf-8") as f:
            raw: dict[str, Any] = json.load(f)

        self._predictions = [
            _parse_prediction(p) for p in raw.get("prediction_tables", [])
        ]
        self._evaluations = [
            _parse_evaluation(e) for e in raw.get("evaluation_tables", [])
        ]
        self._filter_paths: list[str] = list(raw.get("filter_tables", []))
        self._linker: str | None = raw.get("linker_table")
        self._output: str = raw.get("output", "output.parquet")
        self._pipeline = _parse_pipeline(raw.get("pipeline"))
        self._smoothing = _parse_smoothing(raw.get("smoothing"))
        self._config_path = str(path)

        self._validate_paths()

    def _validate_paths(self) -> None:
        """Warn about missing files; don't fail hard (GCS paths won't exist locally)."""
        for entry in self._predictions:
            if not entry.path.startswith("gs://") and not Path(entry.path).exists():
                logger.warning("Prediction table not found: %s", entry.path)
        for entry in self._evaluations:
            if not entry.path.startswith("gs://") and not Path(entry.path).exists():
                logger.warning("Evaluation table not found: %s", entry.path)
        for p in self._filter_paths:
            if not p.startswith("gs://") and not Path(p).exists():
                logger.warning("Filter table not found: %s", p)
        if self._linker and not self._linker.startswith("gs://"):
            if not Path(self._linker).exists():
                logger.warning("Linker table not found: %s", self._linker)

    def list_predictions(self) -> list[PredictionEntry]:
        return list(self._predictions)

    def list_evaluations(self) -> list[EvaluationEntry]:
        return list(self._evaluations)

    def list_filters(self) -> list[str]:
        return list(self._filter_paths)

    def get_prediction(self, path: str) -> PredictionEntry:
        for entry in self._predictions:
            if entry.path == path:
                return entry
        raise KeyError(f"No prediction table with path: {path}")

    def get_prediction_by_score(self, score_col: str) -> PredictionEntry:
        """Find the prediction entry declaring a given score column."""
        for entry in self._predictions:
            if score_col in entry.score_columns:
                return entry
        raise KeyError(f"No prediction table declares score column: {score_col}")

    def get_evaluation(self, path: str) -> EvaluationEntry:
        for entry in self._evaluations:
            if entry.path == path:
                return entry
        raise KeyError(f"No evaluation table with path: {path}")

    @property
    def linker_table(self) -> str | None:
        return self._linker

    @property
    def output(self) -> str:
        return self._output

    @property
    def pipeline(self) -> PipelineConfig:
        return self._pipeline

    @property
    def smoothing(self) -> SmoothingConfig:
        return self._smoothing

    def all_score_columns(self) -> list[str]:
        """Flattened list of all declared score columns across all prediction tables."""
        cols: list[str] = []
        for entry in self._predictions:
            cols.extend(entry.score_columns)
        return cols

    def variant_predictions(self) -> list[PredictionEntry]:
        return [e for e in self._predictions if e.level == "variant"]

    def gene_predictions(self) -> list[PredictionEntry]:
        return [e for e in self._predictions if e.level == "gene"]

    def variant_evaluations(self) -> list[EvaluationEntry]:
        return [e for e in self._evaluations if e.level == "variant"]

    def gene_evaluations(self) -> list[EvaluationEntry]:
        return [e for e in self._evaluations if e.level == "gene"]
