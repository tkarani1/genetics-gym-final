"""Plotly-based interactive charts for eval output DataFrames.

Each function accepts a Polars DataFrame (as returned by
``session.run_stat()`` or ``session.compare()``) and returns a
Plotly ``Figure`` that Streamlit renders via ``st.plotly_chart()``.
"""

from __future__ import annotations

import math
from typing import Any

import plotly.express as px
import plotly.graph_objects as go
import polars as pl


def enrichment_bar(
    df: pl.DataFrame,
    *,
    eval_name: str | None = None,
    filter_name: str = "none",
) -> go.Figure:
    """Grouped bar chart: enrichment across scores at each threshold.

    Filters the results DataFrame to the specified eval/filter combination,
    selects the ``enrichment`` stat, and produces a grouped bar chart with
    one group per threshold and one bar per score.
    """
    filtered = _filter_results(df, stat="enrichment", eval_name=eval_name, filter_name=filter_name)
    if filtered.is_empty():
        return _empty_figure("No enrichment data to plot")

    pdf = filtered.to_pandas()
    fig = px.bar(
        pdf,
        x="threshold",
        y="value",
        color="score_name",
        barmode="group",
        title="Enrichment by Score and Threshold",
        labels={"value": "Enrichment", "threshold": "Percentile Threshold", "score_name": "Score"},
    )
    fig.update_layout(
        xaxis_title="Percentile Threshold",
        yaxis_title="Enrichment",
        legend_title="Score",
        hovermode="x unified",
    )
    return fig


def rate_ratio_bar(
    df: pl.DataFrame,
    *,
    eval_name: str | None = None,
    filter_name: str = "none",
) -> go.Figure:
    """Grouped bar chart for rate_ratio stat."""
    filtered = _filter_results(df, stat="rate_ratio", eval_name=eval_name, filter_name=filter_name)
    if filtered.is_empty():
        return _empty_figure("No rate_ratio data to plot")

    pdf = filtered.to_pandas()
    fig = px.bar(
        pdf,
        x="threshold",
        y="value",
        color="score_name",
        barmode="group",
        title="Rate Ratio by Score and Threshold",
        labels={"value": "Rate Ratio", "threshold": "Percentile Threshold", "score_name": "Score"},
    )
    fig.update_layout(
        xaxis_title="Percentile Threshold",
        yaxis_title="Rate Ratio",
        legend_title="Score",
        hovermode="x unified",
    )
    return fig


def enrichment_dot_with_ci(
    df: pl.DataFrame,
    *,
    eval_name: str | None = None,
    filter_name: str = "none",
) -> go.Figure:
    """Dot plot with bootstrap confidence intervals for enrichment.

    Requires ``std_error`` column (populated by bootstrap in eval).
    """
    filtered = _filter_results(df, stat="enrichment", eval_name=eval_name, filter_name=filter_name)
    if filtered.is_empty():
        return _empty_figure("No enrichment data for CI plot")

    fig = go.Figure()
    scores = filtered.get_column("score_name").unique().sort().to_list()

    for score in scores:
        score_df = filtered.filter(pl.col("score_name") == score).sort("threshold")
        thresholds = score_df.get_column("threshold").to_list()
        values = score_df.get_column("value").to_list()

        if "std_error" in score_df.columns:
            errors = score_df.get_column("std_error").to_list()
            errors_clean = [e if not (isinstance(e, float) and math.isnan(e)) else 0 for e in errors]
        else:
            errors_clean = [0] * len(values)

        fig.add_trace(go.Scatter(
            x=[str(t) for t in thresholds],
            y=values,
            error_y=dict(type="data", array=errors_clean, visible=True),
            mode="markers+lines",
            name=score,
            marker=dict(size=8),
        ))

    fig.update_layout(
        title="Enrichment with Bootstrap CI",
        xaxis_title="Percentile Threshold",
        yaxis_title="Enrichment",
        hovermode="x unified",
    )
    return fig


def score_eval_heatmap(
    df: pl.DataFrame,
    *,
    stat: str = "enrichment",
    threshold: float | None = None,
    filter_name: str = "none",
) -> go.Figure:
    """Heatmap: score x eval enrichment (or other stat) matrix.

    If ``threshold`` is None, uses the highest threshold in the data.
    """
    filtered = df.filter(
        (pl.col("stat") == stat)
        & (pl.col("filter_name") == filter_name)
    )
    if filtered.is_empty():
        return _empty_figure(f"No {stat} data for heatmap")

    if threshold is not None:
        filtered = filtered.filter(pl.col("threshold") == threshold)
    else:
        max_thresh = filtered.get_column("threshold").max()
        if max_thresh is not None and not (isinstance(max_thresh, float) and math.isnan(max_thresh)):
            filtered = filtered.filter(pl.col("threshold") == max_thresh)

    if filtered.is_empty():
        return _empty_figure(f"No {stat} data at the selected threshold")

    pivot = filtered.pivot(
        on="eval_name",
        index="score_name",
        values="value",
    )

    score_names = pivot.get_column("score_name").to_list()
    eval_names = [c for c in pivot.columns if c != "score_name"]
    z_values = pivot.select(eval_names).to_numpy()

    fig = go.Figure(data=go.Heatmap(
        z=z_values,
        x=eval_names,
        y=score_names,
        colorscale="Viridis",
        texttemplate="%{z:.2f}",
        hovertemplate="Score: %{y}<br>Eval: %{x}<br>Value: %{z:.4f}<extra></extra>",
    ))
    fig.update_layout(
        title=f"{stat.replace('_', ' ').title()} Matrix (Score x Eval)",
        xaxis_title="Evaluation Set",
        yaxis_title="Score",
    )
    return fig


def threshold_line(
    df: pl.DataFrame,
    *,
    stat: str = "enrichment",
    eval_name: str | None = None,
    filter_name: str = "none",
) -> go.Figure:
    """Line chart showing stat value across thresholds for each score."""
    filtered = _filter_results(df, stat=stat, eval_name=eval_name, filter_name=filter_name)
    if filtered.is_empty():
        return _empty_figure(f"No {stat} data for line chart")

    pdf = filtered.to_pandas()
    fig = px.line(
        pdf,
        x="threshold",
        y="value",
        color="score_name",
        markers=True,
        title=f"{stat.replace('_', ' ').title()} Across Thresholds",
        labels={"value": stat.replace("_", " ").title(), "threshold": "Percentile Threshold"},
    )
    fig.update_layout(hovermode="x unified")
    return fig


def contingency_summary(
    df: pl.DataFrame,
    *,
    eval_name: str | None = None,
    threshold: float | None = None,
    filter_name: str = "none",
) -> go.Figure:
    """Stacked bar showing TP/FP/TN/FN composition per score at a threshold."""
    filtered = _filter_results(
        df, stat="enrichment", eval_name=eval_name, filter_name=filter_name
    )
    if threshold is not None:
        filtered = filtered.filter(pl.col("threshold") == threshold)
    elif not filtered.is_empty():
        max_thresh = filtered.get_column("threshold").max()
        if max_thresh is not None:
            filtered = filtered.filter(pl.col("threshold") == max_thresh)

    if filtered.is_empty():
        return _empty_figure("No contingency data")

    pdf = filtered.to_pandas()
    fig = go.Figure()
    for component, color in [("tp", "#2ecc71"), ("fp", "#e74c3c"), ("tn", "#3498db"), ("fn", "#e67e22")]:
        if component in pdf.columns:
            fig.add_trace(go.Bar(
                x=pdf["score_name"],
                y=pdf[component],
                name=component.upper(),
                marker_color=color,
            ))

    fig.update_layout(
        barmode="stack",
        title="Contingency Table Composition",
        xaxis_title="Score",
        yaxis_title="Count",
    )
    return fig


# ------------------------------------------------------------------ #
# Helpers
# ------------------------------------------------------------------ #

def _filter_results(
    df: pl.DataFrame,
    stat: str,
    eval_name: str | None = None,
    filter_name: str = "none",
) -> pl.DataFrame:
    """Filter a results DataFrame to a specific stat/eval/filter combination."""
    mask = (pl.col("stat") == stat) & (pl.col("filter_name") == filter_name)
    if eval_name is not None:
        mask = mask & (pl.col("eval_name") == eval_name)

    available = set(df.columns)
    required = {"stat", "filter_name"}
    if not required.issubset(available):
        return pl.DataFrame()

    return df.filter(mask)


def _empty_figure(message: str) -> go.Figure:
    """Return a blank figure with a centered message."""
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        xref="paper", yref="paper",
        x=0.5, y=0.5,
        showarrow=False,
        font=dict(size=16, color="gray"),
    )
    fig.update_layout(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )
    return fig
