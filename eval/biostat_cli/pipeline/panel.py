"""Panel table building and QC summary generation."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import polars as pl

from biostat_cli.types import PanelLayoutConfig, PipelineMode

if TYPE_CHECKING:
    from biostat_cli.pipeline.config import PipelineConfig


def method_label(method: str, display_names: dict[str, str]) -> str:
    """Get display label for a method."""
    return display_names.get(method, method)


def build_panel_table(
    *,
    raw_df: pl.DataFrame | None,
    pairwise_df: pl.DataFrame | None,
    config: PipelineConfig,
    threshold: float,
    eval_set: list[str],
    panel_layout: PanelLayoutConfig,
) -> pl.DataFrame:
    """
    Build the panel table from raw and/or pairwise metrics.

    The panel table filters metrics to the specified threshold and default filter,
    then adds panel metadata (panel_id, panel_title, method labels, ranks).

    Args:
        raw_df: Raw metrics DataFrame (or None).
        pairwise_df: Pairwise metrics DataFrame (or None).
        config: Pipeline configuration.
        threshold: Threshold value to filter on.
        eval_set: List of eval column names to include.
        panel_layout: Panel layout configuration.

    Returns:
        Combined panel table DataFrame.
    """
    panel_by_eval = {eval_name: panel_id for panel_id, eval_name in panel_layout.panel_eval_map.items()}
    method_rank = {name: idx for idx, name in enumerate(config.method_order)}

    frames: list[pl.DataFrame] = []

    for metric_family, df in [("raw", raw_df), ("pairwise", pairwise_df)]:
        if df is None or df.is_empty():
            continue

        expr = (
            pl.col("eval_name").is_in(eval_set)
            & (pl.col("filter_name") == config.default_filter_name)
            & (pl.col("threshold") == threshold)
        )
        subset = df.filter(expr)
        if subset.is_empty():
            continue

        subset = subset.with_columns(
            [
                pl.col("eval_name")
                .map_elements(lambda v: panel_by_eval.get(str(v), "NA"), return_dtype=pl.String)
                .alias("panel_id"),
                pl.col("eval_name")
                .map_elements(
                    lambda v: panel_layout.panel_titles.get(panel_by_eval.get(str(v), ""), str(v)),
                    return_dtype=pl.String,
                )
                .alias("panel_title"),
                pl.lit(metric_family).alias("metric_family"),
                (pl.col("rows_used") / pl.col("total_eval_rows")).alias("rows_used_frac"),
                pl.col("score_name")
                .map_elements(lambda v: method_label(str(v), config.method_display_names), return_dtype=pl.String)
                .alias("method_label"),
                pl.col("score_name")
                .map_elements(lambda v: method_rank.get(str(v), 10_000), return_dtype=pl.Int64)
                .alias("method_rank"),
            ]
        )
        frames.append(subset)

    if not frames:
        return pl.DataFrame()

    panel_df = pl.concat(frames, how="diagonal_relaxed").sort(
        ["metric_family", "panel_id", "stat", "method_rank", "score_name"]
    )
    return panel_df


def build_qc_summary(
    panel_df: pl.DataFrame,
    config: PipelineConfig,
    mode: PipelineMode,
    panel_order: list[str],
    panel_metrics: dict[str, dict[str, str]],
) -> tuple[pl.DataFrame, list[str]]:
    """
    Build QC summary checking panel row counts and pairwise anchor consistency.

    Args:
        panel_df: Panel table DataFrame.
        config: Pipeline configuration.
        mode: Pipeline mode.
        panel_order: Ordered list of panel IDs.
        panel_metrics: Per-panel metric configuration.

    Returns:
        Tuple of (QC summary DataFrame, list of warning messages).
    """
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    families = ["raw", "pairwise"] if mode == PipelineMode.BOTH else [mode.value]

    for family in families:
        fam_df = panel_df.filter(pl.col("metric_family") == family)
        for panel_id in panel_order:
            panel_stat = panel_metrics[panel_id][family]
            count = fam_df.filter(
                (pl.col("panel_id") == panel_id) & (pl.col("stat") == panel_stat)
            ).height

            rows.append(
                {
                    "check": "panel_stat_row_count",
                    "metric_family": family,
                    "panel_id": panel_id,
                    "stat": panel_stat,
                    "value": float(count),
                    "status": "ok" if count > 0 else "warning",
                }
            )
            if count == 0:
                warnings.append(f"No rows for panel `{panel_id}` with stat `{panel_stat}` in `{family}` mode.")

    # Check pairwise anchor adjustment ratio
    if mode in {PipelineMode.PAIRWISE, PipelineMode.BOTH} and not panel_df.is_empty():
        pw = panel_df.filter(pl.col("metric_family") == "pairwise")
        if "adjustment_ratio" in pw.columns:
            anchor_rows = pw.filter(
                pl.col("score_name").is_in(["mpc_score", "mpc_score_percentile", "mpc_score_anchor_percentile"])
            )
            if anchor_rows.height > 0:
                deviations = (
                    anchor_rows.with_columns((pl.col("adjustment_ratio") - 1.0).abs().alias("abs_dev"))
                    .select(pl.col("abs_dev").max())
                    .item()
                )
                status = "ok" if deviations <= 1e-6 else "warning"
                rows.append(
                    {
                        "check": "pairwise_anchor_adjustment_ratio_max_abs_deviation",
                        "metric_family": "pairwise",
                        "panel_id": "all",
                        "stat": "pairwise_*",
                        "value": float(deviations),
                        "status": status,
                    }
                )
                if status != "ok":
                    warnings.append("Pairwise anchor adjustment ratio is not consistently ~1.")

    qc_df = pl.DataFrame(rows) if rows else pl.DataFrame(schema={"check": pl.String, "status": pl.String})
    return qc_df, warnings


def write_qc_report(
    qc_summary: pl.DataFrame,
    warnings: list[str],
    out_path: str,
    *,
    config_path: str,
    outdir: str,
) -> None:
    """
    Write markdown QC report file.

    Args:
        qc_summary: QC summary DataFrame.
        warnings: List of warning messages.
        out_path: Output file path.
        config_path: Path to pipeline config file.
        outdir: Output directory path.
    """
    lines = [
        "# Figure 1 Pipeline QC Report",
        "",
        f"- config: `{config_path}`",
        f"- outdir: `{outdir}`",
        "",
        "## Warnings",
    ]

    if warnings:
        for warning in warnings:
            lines.append(f"- {warning}")
    else:
        lines.append("- none")

    lines.extend(["", "## Checks", ""])

    if qc_summary.is_empty():
        lines.append("No checks were produced.")
    else:
        for row in qc_summary.to_dicts():
            lines.append(
                f"- [{row.get('status', 'unknown')}] {row.get('check')} | "
                f"family={row.get('metric_family')} panel={row.get('panel_id')} "
                f"stat={row.get('stat')} value={row.get('value')}"
            )

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text("\n".join(lines), encoding="utf-8")


def resolve_panel_df_for_plot(
    args: Any,  # PipelineArgs
    config: PipelineConfig,
    outdir: str,
    eval_set: list[str] | None = None,
    panel_layout: PanelLayoutConfig | None = None,
) -> tuple[pl.DataFrame, str]:
    """
    Resolve panel DataFrame for plotting from various input sources.

    Args:
        args: Pipeline arguments.
        config: Pipeline configuration.
        outdir: Output directory.
        eval_set: Optional eval set override.
        panel_layout: Optional panel layout override.

    Returns:
        Tuple of (panel DataFrame, panel table path).

    Raises:
        ValueError: If no valid input source is provided.
    """
    from biostat_cli.io import write_tsv

    if args.panel_table:
        return pl.read_csv(args.panel_table, separator="\t"), args.panel_table

    raw_df = pl.read_csv(args.raw_metrics_tsv, separator="\t") if args.raw_metrics_tsv else None
    pairwise_df = pl.read_csv(args.pairwise_metrics_tsv, separator="\t") if args.pairwise_metrics_tsv else None

    if raw_df is None and pairwise_df is None:
        raise ValueError(
            "plot requires one of: --panel-table OR --raw-metrics-tsv/--pairwise-metrics-tsv."
        )

    threshold_for_panels = args.threshold if args.threshold is not None else config.default_threshold
    layout = panel_layout or config.panel_layout

    panel_df = build_panel_table(
        raw_df=raw_df,
        pairwise_df=pairwise_df,
        config=config,
        threshold=threshold_for_panels,
        eval_set=eval_set or config.eval_set,
        panel_layout=layout,
    )

    panel_path = str(Path(outdir) / "panel_table.tsv")
    write_tsv(panel_df, panel_path)

    return panel_df, panel_path


__all__ = [
    "method_label",
    "build_panel_table",
    "build_qc_summary",
    "write_qc_report",
    "resolve_panel_df_for_plot",
]
