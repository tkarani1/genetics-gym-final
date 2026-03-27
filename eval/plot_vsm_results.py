#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PLOT_STATS = {"enrichment", "rate_ratio", "pairwise_enrichment", "pairwise_rate_ratio"}


def _safe_name(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in text)


def _pretty_score_name(text: str) -> str:
    text = text.replace("_", " ")
    return " ".join(token for token in text.split() if token.lower() != "percentile")


def _y_axis_label(stat: str) -> str:
    if stat.endswith("rate_ratio"):
        return "Rate ratio value"
    return "Enrichment value"


def _rows_used_pct_text(rows_used: object, total_eval_rows: object) -> str:
    rows_used_num = pd.to_numeric(pd.Series([rows_used]), errors="coerce").iloc[0]
    total_eval_num = pd.to_numeric(pd.Series([total_eval_rows]), errors="coerce").iloc[0]
    if pd.isna(rows_used_num) or pd.isna(total_eval_num) or not total_eval_num:
        return "NA"
    return f"{(rows_used_num / total_eval_num) * 100:.1f}%"


def _anchor_row_index_for_pairwise(group: pd.DataFrame) -> int | None:
    """Identify anchor row in pairwise outputs via adjustment_ratio ~= 1."""
    adjustment_vals = pd.to_numeric(group["adjustment_ratio"], errors="coerce")
    for idx, value in enumerate(adjustment_vals.tolist()):
        if pd.notna(value) and abs(float(value) - 1.0) < 1e-9:
            return idx
    return None


def plot_enrichment_by_group(input_tsv: Path, output_dir: Path) -> None:
    df = pd.read_csv(input_tsv, sep="\t")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df["threshold"] = df["threshold"].astype(str)

    # Plot enrichment/rate-ratio and their pairwise-adjusted counterparts as dot plots.
    plot_df = df[df["stat"].isin(PLOT_STATS) & df["value"].notna()].copy()
    group_cols = ["eval_name", "filter_name", "stat", "threshold"]
    output_dir.mkdir(parents=True, exist_ok=True)

    for keys, group in plot_df.groupby(group_cols, dropna=False):
        eval_name, filter_name, stat, threshold = keys
        group = group.sort_values("value", ascending=False)

        fig, ax = plt.subplots(figsize=(10, 5))
        x_labels = [_pretty_score_name(v) for v in group["score_name"].astype(str).tolist()]
        if stat.startswith("pairwise_"):
            anchor_idx = _anchor_row_index_for_pairwise(group)
            pairwise_labels: list[str] = []
            for idx, (method_name, rows_used, total_eval) in enumerate(
                zip(x_labels, group["rows_used"].tolist(), group["total_eval_rows"].tolist())
            ):
                method_label = f"{method_name}*" if anchor_idx is not None and idx == anchor_idx else method_name
                pairwise_labels.append(f"{method_label}\n{_rows_used_pct_text(rows_used, total_eval)}")
            x_labels = pairwise_labels
        x_pos = range(len(x_labels))
        y_vals = group["value"].tolist()
        same_rows_used_for_rate = False
        common_rows_used_pct_text = "NA"

        ax.scatter(x_pos, y_vals, s=70)
        ax.set_xticks(list(x_pos))
        if stat.startswith("pairwise_"):
            ax.set_xticklabels(x_labels, rotation=0, ha="center")
        else:
            ax.set_xticklabels(x_labels, rotation=45, ha="right")
        ax.set_ylabel(_y_axis_label(stat))
        if stat.endswith("enrichment"):
            ax.set_ylim(bottom=0)
        ax.set_xlabel("Score")
        ax.set_title(
            f"eval={eval_name} | filter={filter_name} | stat={stat} | threshold={threshold}"
        )
        ax.grid(axis="y", alpha=0.25)

        if stat.endswith("rate_ratio"):
            row_used_vals = pd.to_numeric(group["rows_used"], errors="coerce")
            total_eval_vals = pd.to_numeric(group["total_eval_rows"], errors="coerce")
            max_total_eval = total_eval_vals.max(skipna=True)
            same_rows_used_for_rate = row_used_vals.nunique(dropna=True) == 1
            y_span = max(y_vals) - min(y_vals) if y_vals else 0.0
            y_offset = (0.03 * y_span) if y_span > 0 else 0.05

            if same_rows_used_for_rate:
                common_rows_used = row_used_vals.dropna()
                common_total_eval = total_eval_vals.dropna()
                if not common_rows_used.empty and not common_total_eval.empty:
                    total_ref = common_total_eval.iloc[0]
                    if total_ref:
                        common_rows_used_pct_text = f"{(common_rows_used.iloc[0] / total_ref) * 100:.1f}%"

            for x, y, rows_used, total_eval in zip(x_pos, y_vals, row_used_vals, total_eval_vals):
                if pd.isna(y):
                    continue

                total_pct_text = "NA"
                if pd.notna(total_eval) and pd.notna(max_total_eval) and max_total_eval:
                    total_pct_text = f"{(total_eval / max_total_eval) * 100:.1f}%"

                if same_rows_used_for_rate:
                    point_label = f"total:{total_pct_text}"
                else:
                    used_pct_text = "NA"
                    if pd.notna(rows_used) and pd.notna(total_eval) and total_eval:
                        used_pct_text = f"{(rows_used / total_eval) * 100:.1f}%"
                    point_label = f"used:{used_pct_text}\ntotal:{total_pct_text}"

                ax.text(
                    x,
                    y + y_offset,
                    point_label,
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )

        rows_used_values = group["rows_used"].dropna().unique().tolist()
        total_eval_rows_values = group["total_eval_rows"].dropna().unique().tolist()
        rows_used_text = ", ".join(str(int(v)) for v in rows_used_values) or "NA"
        total_rows_text = ", ".join(str(int(v)) for v in total_eval_rows_values) or "NA"
        note = f"rows_used: {rows_used_text}\ntotal_eval_rows: {total_rows_text}"
        if stat.endswith("rate_ratio") and same_rows_used_for_rate:
            note = f"{note}\nrows_used_pct: {common_rows_used_pct_text}"
        if stat.startswith("pairwise_"):
            anchor_values = pd.to_numeric(group["anchor_value"], errors="coerce").dropna().unique().tolist()
            adjustment_values = (
                pd.to_numeric(group["adjustment_ratio"], errors="coerce").dropna().unique().tolist()
            )
            anchor_text = ", ".join(f"{v:.4g}" for v in anchor_values) or "NA"
            adjustment_text = ", ".join(f"{v:.4g}" for v in adjustment_values) or "NA"
            note = f"{note}\nanchor_value: {anchor_text}\nadjustment_ratio: {adjustment_text}"

        ax.text(
            0.98,
            0.98,
            note,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.9},
        )

        fig.tight_layout()

        file_name = "_".join(
            _safe_name(str(part)) for part in (eval_name, filter_name, stat, threshold)
        )
        out_path = output_dir / f"{file_name}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"Saved: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plot enrichment/rate-ratio dot plots (including pairwise-adjusted stats) "
            "grouped by eval_name/filter_name/stat/threshold."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to input TSV file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/Users/tk508/Work/new/results/vsm_plots"),
        help="Directory for output PNG plots.",
    )
    args = parser.parse_args()

    plot_enrichment_by_group(args.input, args.output_dir)


if __name__ == "__main__":
    main()
