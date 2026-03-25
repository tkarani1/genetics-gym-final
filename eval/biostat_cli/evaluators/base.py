from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import polars as pl


@dataclass(frozen=True)
class PreparedFrame:
    frame: pl.LazyFrame
    total_eval_rows: int


@dataclass(frozen=True)
class ScoreFrame:
    frame: pl.LazyFrame
    rows_used: int


@dataclass(frozen=True)
class Contingency:
    tp: float
    fp: float
    tn: float
    fn: float


class BaseEvaluator(ABC):
    def __init__(self, source: pl.LazyFrame) -> None:
        self.source = source

    def requires_eval_non_null(self, eval_col: str) -> bool:
        return True

    def prepare_eval_frame(self, eval_col: str, filter_col: str | None) -> PreparedFrame:
        lf = self.source
        conditions: list[pl.Expr] = []
        if self.requires_eval_non_null(eval_col):
            conditions.append(pl.col(eval_col).is_not_null())
        if filter_col:
            conditions.append(pl.col(filter_col) == True)  # noqa: E712
        if conditions:
            lf = lf.filter(pl.all_horizontal(conditions))
        total_eval_rows = int(lf.select(pl.len().alias("n")).collect(streaming=True)["n"][0])
        return PreparedFrame(frame=lf, total_eval_rows=total_eval_rows)

    def prepare_score_frame(self, prepared: PreparedFrame, score_col: str) -> ScoreFrame:
        df = prepared.frame.filter(pl.col(score_col).is_not_null()).collect(streaming=True)
        return ScoreFrame(frame=df.lazy(), rows_used=df.height)

    @abstractmethod
    def contingency(self, score_frame: ScoreFrame, eval_col: str, score_col: str, threshold: float) -> Contingency:
        raise NotImplementedError

    @abstractmethod
    def contingency_batch(
        self, score_frame: ScoreFrame, eval_col: str, score_col: str, thresholds: list[float]
    ) -> list[Contingency]:
        raise NotImplementedError

    @abstractmethod
    def labels_and_scores(
        self, score_frame: ScoreFrame, eval_col: str, score_col: str
    ) -> tuple[list[int], list[float]] | None:
        raise NotImplementedError
