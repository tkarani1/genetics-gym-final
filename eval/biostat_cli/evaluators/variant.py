from __future__ import annotations

import polars as pl

from biostat_cli.evaluators.base import BaseEvaluator, Contingency, ScoreFrame


class VariantEvaluator(BaseEvaluator):
    def contingency(self, score_frame: ScoreFrame, eval_col: str, score_col: str, threshold: float) -> Contingency:
        return self.contingency_batch(score_frame, eval_col, score_col, [threshold])[0]

    def contingency_batch(
        self, score_frame: ScoreFrame, eval_col: str, score_col: str, thresholds: list[float]
    ) -> list[Contingency]:
        if not thresholds:
            return []
        is_pos = pl.col(eval_col) == True  # noqa: E712
        is_neg = pl.col(eval_col) == False  # noqa: E712
        exprs: list[pl.Expr] = []
        for i, t in enumerate(thresholds):
            above = pl.col(score_col) > t
            exprs.extend([
                pl.when(above & is_pos).then(1).otherwise(0).sum().cast(pl.Float64).alias(f"tp_{i}"),
                pl.when(above & is_neg).then(1).otherwise(0).sum().cast(pl.Float64).alias(f"fp_{i}"),
                pl.when((~above) & is_neg).then(1).otherwise(0).sum().cast(pl.Float64).alias(f"tn_{i}"),
                pl.when((~above) & is_pos).then(1).otherwise(0).sum().cast(pl.Float64).alias(f"fn_{i}"),
            ])
        row = score_frame.frame.select(exprs).collect(streaming=True).to_dicts()[0]
        return [
            Contingency(tp=row[f"tp_{i}"], fp=row[f"fp_{i}"], tn=row[f"tn_{i}"], fn=row[f"fn_{i}"])
            for i in range(len(thresholds))
        ]

    def labels_and_scores(
        self, score_frame: ScoreFrame, eval_col: str, score_col: str
    ) -> tuple[list[int], list[float]] | None:
        out = score_frame.frame.select(
            pl.col(eval_col).cast(pl.Int64).alias("label"),
            pl.col(score_col).cast(pl.Float64).alias("score"),
        ).collect(streaming=True)
        labels = out["label"].to_list()
        scores = out["score"].to_list()
        return labels, scores
