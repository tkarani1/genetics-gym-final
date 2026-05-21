from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import polars as pl

from biostat_cli.utils import WITHIN_GENE_COL, apply_within_gene_percentile


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


def slice_prepared_for_score(
    prepared: PreparedFrame,
    *,
    eval_col: str,
    filter_col: str | None,
    score_col: str,
    within_gene_percentile: bool,
    gene_col: str | None,
    schema_names: set[str],
    obs_exp_columns: tuple[str, str] | None = None,
) -> PreparedFrame:
    """Project ``prepared`` to columns needed for one score (avoids loading wide rows on collect)."""
    if obs_exp_columns is not None:
        obs_c, exp_c = obs_exp_columns
        want: list[str] = [obs_c, exp_c, score_col]
    else:
        want: list[str] = [eval_col, score_col]
    # GeneEvaluator.sum_variants contingencies need per-row burden counts on the slice.
    if not obs_exp_columns and eval_col == "sum_variants":
        for c in ("n_case", "n_ctrl"):
            if c in schema_names:
                want.append(c)
    if filter_col:
        want.append(filter_col)
    if within_gene_percentile:
        want.append(WITHIN_GENE_COL)
    if gene_col:
        want.append(gene_col)
    sel: list[str] = []
    seen: set[str] = set()
    for c in want:
        if c not in schema_names:
            raise KeyError(f"Column {c!r} not found in table (while processing score {score_col!r}).")
        if c not in seen:
            sel.append(c)
            seen.add(c)
    return PreparedFrame(
        frame=prepared.frame.select(sel),
        total_eval_rows=prepared.total_eval_rows,
    )


class BaseEvaluator(ABC):
    def __init__(self, source: pl.LazyFrame) -> None:
        self.source = source

    def requires_eval_non_null(self, eval_col: str) -> bool:
        return True

    def prepare_eval_frame(
        self,
        eval_col: str,
        filter_col: str | None,
        *,
        obs_exp_columns: tuple[str, str] | None = None,
    ) -> PreparedFrame:
        lf = self.source
        conditions: list[pl.Expr] = []
        if obs_exp_columns is not None:
            obs_c, exp_c = obs_exp_columns
            conditions.append(pl.col(obs_c).is_not_null() & pl.col(exp_c).is_not_null())
        else:
            if self.requires_eval_non_null(eval_col):
                conditions.append(pl.col(eval_col).is_not_null())
        if filter_col:
            conditions.append(pl.col(filter_col) == True)  # noqa: E712
        if conditions:
            lf = lf.filter(pl.all_horizontal(conditions))
        total_eval_rows = int(lf.select(pl.len().alias("n")).collect(streaming=True)["n"][0])
        return PreparedFrame(frame=lf, total_eval_rows=total_eval_rows)

    def prepare_score_frame(
        self, prepared: PreparedFrame, score_col: str, *, within_gene_percentile: bool = False
    ) -> ScoreFrame:
        lf = prepared.frame
        if within_gene_percentile:
            lf = apply_within_gene_percentile(lf, score_col)
        df = lf.filter(pl.col(score_col).is_not_null()).collect(streaming=True)
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

    def contingency_by_gene_batch(
        self,
        score_frame: ScoreFrame,
        eval_col: str,
        score_col: str,
        thresholds: list[float],
        gene_col: str,
    ) -> dict[str, list[Contingency]]:
        raise NotImplementedError(f"{type(self).__name__} does not support per-gene contingency.")

    def labels_and_scores_by_gene(
        self,
        score_frame: ScoreFrame,
        eval_col: str,
        score_col: str,
        gene_col: str,
    ) -> dict[str, tuple[list[int], list[float]]]:
        raise NotImplementedError(f"{type(self).__name__} does not support per-gene labels/scores.")
