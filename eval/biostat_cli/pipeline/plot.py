"""Plotting functions for Figure 1 pipeline."""

from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import polars as pl

from biostat_cli.types import PanelLayoutConfig, PipelineMode

if TYPE_CHECKING:
    from biostat_cli.pipeline.config import PipelineConfig


def format_method_tick(method_label: str, rows_used_frac: float | None) -> str:
    """Format method label with optional coverage percentage."""
    if rows_used_frac is None or (isinstance(rows_used_frac, float) and math.isnan(rows_used_frac)):
        return method_label
    return f"{method_label}\n({rows_used_frac * 100:.1f}%)"


def render_mode_figure(
    panel_df: pl.DataFrame,
    *,
    mode: str,
    config: PipelineConfig,
    out_png: str,
    out_pdf: str,
    panel_order: list[str],
    panel_titles: dict[str, str],
    panel_metrics: dict[str, dict[str, str]],
) -> None:
    """
    Render figure for a single mode (raw or pairwise).

    Args:
        panel_df: Panel table DataFrame.
        mode: Mode to render ("raw" or "pairwise").
        config: Pipeline configuration.
        out_png: Output PNG path.
        out_pdf: Output PDF path.
        panel_order: Ordered list of panel IDs.
        panel_titles: Panel title mapping.
        panel_metrics: Per-panel metric configuration.
    """
    df = panel_df.filter(pl.col("metric_family") == mode)
    n_panels = len(panel_order)
    n_cols = 2
    n_rows = max(2, (n_panels + 1) // 2)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 4.5 * n_rows), constrained_layout=True)
    flat_axes = list(axes.ravel()) if n_rows > 1 else [axes[0], axes[1]] if n_cols > 1 else [axes]

    for idx, panel_id in enumerate(panel_order):
        if idx >= len(flat_axes):
            break

        ax = flat_axes[idx]
        panel_stat = panel_metrics[panel_id][mode]
        sub = df.filter(
            (pl.col("panel_id") == panel_id) & (pl.col("stat") == panel_stat)
        ).sort(["method_rank", "score_name"])

        if sub.is_empty():
            ax.set_title(f"{panel_id}: {panel_titles.get(panel_id, panel_id)}")
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            ax.set_xticks([])
            continue

        methods = sub["method_label"].to_list()
        values = sub["value"].to_list()
        stderrs = sub["std_error"].to_list()
        rows_used_frac = sub["rows_used_frac"].to_list()

        x = list(range(len(methods)))
        yerr = [float("nan") if (v is None or (isinstance(v, float) and math.isnan(v))) else v for v in stderrs]

        ax.errorbar(x, values, yerr=yerr, fmt="o", capsize=3)
        ax.axhline(1.0, linestyle="--", linewidth=1)
        ax.set_xticks(x)
        ax.set_xticklabels(
            [format_method_tick(m, r) for m, r in zip(methods, rows_used_frac)],
            rotation=35,
            ha="right",
        )
        ax.set_ylabel(panel_stat)
        ax.set_title(f"{panel_id}: {panel_titles.get(panel_id, panel_id)}")

    # Hide unused axes
    for idx in range(n_panels, len(flat_axes)):
        flat_axes[idx].axis("off")

    fig.suptitle(f"Figure 1-style panels ({mode})")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_pdf)
    plt.close(fig)


def render_combined_figure(
    panel_df: pl.DataFrame,
    config: PipelineConfig,
    out_png: str,
    out_pdf: str,
    panel_order: list[str],
    panel_titles: dict[str, str],
    panel_metrics: dict[str, dict[str, str]],
) -> None:
    """
    Render combined figure with raw and pairwise side-by-side.

    Args:
        panel_df: Panel table DataFrame.
        config: Pipeline configuration.
        out_png: Output PNG path.
        out_pdf: Output PDF path.
        panel_order: Ordered list of panel IDs.
        panel_titles: Panel title mapping.
        panel_metrics: Per-panel metric configuration.
    """
    families = [f for f in ["raw", "pairwise"] if f in panel_df["metric_family"].unique().to_list()]
    n_cols = len(families)
    n_rows = len(panel_order)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(7 * max(n_cols, 1), 2.7 * n_rows),
        squeeze=False,
        constrained_layout=True,
    )

    for r, panel_id in enumerate(panel_order):
        for c, family in enumerate(families):
            ax = axes[r][c]
            panel_stat = panel_metrics[panel_id][family]
            sub = panel_df.filter(
                (pl.col("metric_family") == family)
                & (pl.col("panel_id") == panel_id)
                & (pl.col("stat") == panel_stat)
            ).sort(["method_rank", "score_name"])

            if sub.is_empty():
                ax.text(0.5, 0.5, "No data", ha="center", va="center")
                ax.set_xticks([])
                ax.set_title(f"{panel_id} {family}")
                continue

            methods = sub["method_label"].to_list()
            values = sub["value"].to_list()
            x = list(range(len(methods)))

            ax.scatter(x, values)
            ax.axhline(1.0, linestyle="--", linewidth=1)
            ax.set_xticks(x)
            ax.set_xticklabels(methods, rotation=35, ha="right")
            ax.set_ylabel(panel_stat)
            ax.set_title(f"{panel_id}: {panel_titles.get(panel_id, panel_id)} [{family}]")

    fig.suptitle("Figure 1-style panels (combined)")
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_pdf)
    plt.close(fig)


def render_plots(
    panel_df: pl.DataFrame,
    config: PipelineConfig,
    outdir: str,
    mode: PipelineMode,
    output_layout: str,
    panel_layout: PanelLayoutConfig,
) -> list[str]:
    """
    Render all plots based on mode and output layout.

    Args:
        panel_df: Panel table DataFrame.
        config: Pipeline configuration.
        outdir: Output directory.
        mode: Pipeline mode.
        output_layout: Output layout ("combined", "per_eval", or "both").
        panel_layout: Panel layout configuration.

    Returns:
        List of output file paths.
    """
    outputs: list[str] = []
    panel_order = panel_layout.panel_order
    panel_titles = panel_layout.panel_titles
    panel_metrics = panel_layout.panel_metrics

    if output_layout in {"combined", "both"}:
        outputs.extend(
            _render_combined_outputs(
                panel_df, config, outdir, mode, panel_order, panel_titles, panel_metrics
            )
        )

    if output_layout in {"per_eval", "both"}:
        outputs.extend(
            _render_per_eval_outputs(
                panel_df, config, outdir, mode, panel_order, panel_titles, panel_metrics
            )
        )

    return outputs


def _render_combined_outputs(
    panel_df: pl.DataFrame,
    config: PipelineConfig,
    outdir: str,
    mode: PipelineMode,
    panel_order: list[str],
    panel_titles: dict[str, str],
    panel_metrics: dict[str, dict[str, str]],
) -> list[str]:
    """Render combined output plots."""
    outputs: list[str] = []

    if mode.includes_raw():
        raw_png = str(Path(outdir) / "figure1_raw.png")
        raw_pdf = str(Path(outdir) / "figure1_raw.pdf")
        render_mode_figure(
            panel_df,
            mode="raw",
            config=config,
            out_png=raw_png,
            out_pdf=raw_pdf,
            panel_order=panel_order,
            panel_titles=panel_titles,
            panel_metrics=panel_metrics,
        )
        outputs.extend([raw_png, raw_pdf])

    if mode.includes_pairwise():
        pw_png = str(Path(outdir) / "figure1_pairwise.png")
        pw_pdf = str(Path(outdir) / "figure1_pairwise.pdf")
        render_mode_figure(
            panel_df,
            mode="pairwise",
            config=config,
            out_png=pw_png,
            out_pdf=pw_pdf,
            panel_order=panel_order,
            panel_titles=panel_titles,
            panel_metrics=panel_metrics,
        )
        outputs.extend([pw_png, pw_pdf])

    if mode == PipelineMode.BOTH:
        cmb_png = str(Path(outdir) / "figure1_combined.png")
        cmb_pdf = str(Path(outdir) / "figure1_combined.pdf")
        render_combined_figure(
            panel_df,
            config=config,
            out_png=cmb_png,
            out_pdf=cmb_pdf,
            panel_order=panel_order,
            panel_titles=panel_titles,
            panel_metrics=panel_metrics,
        )
        outputs.extend([cmb_png, cmb_pdf])

    return outputs


def _render_per_eval_outputs(
    panel_df: pl.DataFrame,
    config: PipelineConfig,
    outdir: str,
    mode: PipelineMode,
    panel_order: list[str],
    panel_titles: dict[str, str],
    panel_metrics: dict[str, dict[str, str]],
) -> list[str]:
    """Render per-eval output plots."""
    outputs: list[str] = []
    per_eval_dir = Path(outdir) / "per_eval"
    per_eval_dir.mkdir(parents=True, exist_ok=True)

    for panel_id in panel_order:
        single_panel_order = [panel_id]
        single_panel_titles = {panel_id: panel_titles.get(panel_id, panel_id)}
        single_panel_metrics = {
            panel_id: panel_metrics.get(panel_id, {"raw": "enrichment", "pairwise": "pairwise_enrichment"})
        }
        single_panel_df = panel_df.filter(pl.col("panel_id") == panel_id)

        if single_panel_df.is_empty():
            continue

        safe_panel_id = panel_id.replace("/", "_").replace(" ", "_")

        if mode.includes_raw():
            raw_png = str(per_eval_dir / f"figure1_raw_{safe_panel_id}.png")
            raw_pdf = str(per_eval_dir / f"figure1_raw_{safe_panel_id}.pdf")
            render_mode_figure(
                single_panel_df,
                mode="raw",
                config=config,
                out_png=raw_png,
                out_pdf=raw_pdf,
                panel_order=single_panel_order,
                panel_titles=single_panel_titles,
                panel_metrics=single_panel_metrics,
            )
            outputs.extend([raw_png, raw_pdf])

        if mode.includes_pairwise():
            pw_png = str(per_eval_dir / f"figure1_pairwise_{safe_panel_id}.png")
            pw_pdf = str(per_eval_dir / f"figure1_pairwise_{safe_panel_id}.pdf")
            render_mode_figure(
                single_panel_df,
                mode="pairwise",
                config=config,
                out_png=pw_png,
                out_pdf=pw_pdf,
                panel_order=single_panel_order,
                panel_titles=single_panel_titles,
                panel_metrics=single_panel_metrics,
            )
            outputs.extend([pw_png, pw_pdf])

    return outputs


__all__ = [
    "format_method_tick",
    "render_mode_figure",
    "render_combined_figure",
    "render_plots",
]
