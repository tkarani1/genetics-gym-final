from __future__ import annotations

import polars as pl

from biostat_cli.evaluators.base import BaseEvaluator, Contingency, ScoreFrame


SUM_VARIANTS_SENTINEL = "sum_variants"


class GeneEvaluator(BaseEvaluator):
    @staticmethod
    def _resolve_burden_cols(eval_col: str, schema_names: list[str]) -> tuple[str, str]:
        schema_set = set(schema_names)
        specific_case = f"n_case_{eval_col}"
        specific_ctrl = f"n_ctrl_{eval_col}"
        if specific_case in schema_set and specific_ctrl in schema_set:
            return specific_case, specific_ctrl
        if "n_case" in schema_set and "n_ctrl" in schema_set:
            return "n_case", "n_ctrl"

        if eval_col == SUM_VARIANTS_SENTINEL:
            matched_cases = [c for c in schema_names if c.startswith("n_case_")]
            candidates: list[tuple[str, str]] = []
            for case_col in matched_cases:
                stem = case_col[len("n_case_") :]
                ctrl_col = f"n_ctrl_{stem}"
                if ctrl_col in schema_set:
                    candidates.append((case_col, ctrl_col))
            if len(candidates) == 1:
                return candidates[0]
            if len(candidates) > 1:
                pairs = ", ".join(f"{a}/{b}" for a, b in candidates)
                raise ValueError(
                    "Ambiguous burden columns for sum_variants: found multiple n_case_*/n_ctrl_* pairs: "
                    f"{pairs}. Use an eval-specific run (eval_set=<stem>) or provide n_case/n_ctrl."
                )

        raise ValueError(
            f"Burden mode requires n_case/n_ctrl or n_case_{eval_col}/n_ctrl_{eval_col} columns."
        )

    @classmethod
    def _is_burden_eval(cls, eval_col: str, schema_names: list[str]) -> bool:
        schema_set = set(schema_names)
        if eval_col == SUM_VARIANTS_SENTINEL:
            return True
        return f"n_case_{eval_col}" in schema_set and f"n_ctrl_{eval_col}" in schema_set

    def requires_eval_non_null(self, eval_col: str) -> bool:
        schema_names = self.source.collect_schema().names()
        return not self._is_burden_eval(eval_col, schema_names)

    def contingency(self, score_frame: ScoreFrame, eval_col: str, score_col: str, threshold: float) -> Contingency:
        return self.contingency_batch(score_frame, eval_col, score_col, [threshold])[0]

    def contingency_batch(
        self, score_frame: ScoreFrame, eval_col: str, score_col: str, thresholds: list[float]
    ) -> list[Contingency]:
        if not thresholds:
            return []
        exprs: list[pl.Expr] = []
        schema_names = score_frame.frame.collect_schema().names()
        if self._is_burden_eval(eval_col, schema_names):
            case_col, ctrl_col = self._resolve_burden_cols(eval_col, schema_names)
            for i, t in enumerate(thresholds):
                above = pl.col(score_col) >= t
                exprs.extend([
                    pl.when(above).then(pl.col(case_col)).otherwise(0).sum().cast(pl.Float64).alias(f"tp_{i}"),
                    pl.when(above).then(pl.col(ctrl_col)).otherwise(0).sum().cast(pl.Float64).alias(f"fp_{i}"),
                    pl.when(~above).then(pl.col(ctrl_col)).otherwise(0).sum().cast(pl.Float64).alias(f"tn_{i}"),
                    pl.when(~above).then(pl.col(case_col)).otherwise(0).sum().cast(pl.Float64).alias(f"fn_{i}"),
                ])
        else:
            is_pos = pl.col(eval_col) == True  # noqa: E712
            is_neg = pl.col(eval_col) == False  # noqa: E712
            for i, t in enumerate(thresholds):
                above = pl.col(score_col) >= t
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
        schema_names = score_frame.frame.collect_schema().names()
        if self._is_burden_eval(eval_col, schema_names):
            return None
        out = score_frame.frame.select(
            pl.col(eval_col).cast(pl.Int64).alias("label"),
            pl.col(score_col).cast(pl.Float64).alias("score"),
        ).collect(streaming=True)
        labels = out["label"].to_list()
        scores = out["score"].to_list()
        return labels, scores
