"""
Figure 2 Pipeline — Gene-level evaluation across cohorts and aggregation types.

For each aggregation type (mean / max) and each cohort (dd / asd / chd / ...):
  1. Rename cohort-specific n_case/n_ctrl columns to generic names.
  2. Normalize score column names (strip mean/max infix) so the existing
     pairwise detection regex works unchanged.
  3. Call biostat_cli.cli.run() with eval_level="gene".
  4. Concat results, build panel tables (one per aggregation type), and render plots.

Usage:
    figure2-pipeline run   --config figure2_pipeline_config.json --raw-parquet ... [--pairwise-parquet ...]
    figure2-pipeline compute --config figure2_pipeline_config.json --raw-parquet ... [--pairwise-parquet ...]
    figure2-pipeline plot  --config figure2_pipeline_config.json --panel-table-mean ... --panel-table-max ...
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import polars as pl

from biostat_cli.cli import RunArgs, run as biostat_run
from biostat_cli.io import write_json, write_tsv
from biostat_cli.pipeline.plot import render_mode_figure, render_combined_figure
from biostat_cli.stats.binary import DEFAULT_PVALUE_METHOD, PVALUE_METHODS
from biostat_cli.types import PipelineMode

# ---------------------------------------------------------------------------
# Config dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CohortConfig:
    name: str
    n_case_col: str
    n_ctrl_col: str
    stats: str
    case_total: float | None
    ctrl_total: float | None


@dataclass(frozen=True)
class PairwiseScoreSpec:
    vsm_cols: list[str]
    anchor_full: str
    anchor_pairwise: dict[str, str]


@dataclass(frozen=True)
class GenePipelineConfig:
    cohorts: list[CohortConfig]
    aggregation_types: list[str]
    raw_score_columns: dict[str, list[str]]
    pairwise_score_specs: dict[str, PairwiseScoreSpec]
    normalized_score_columns: list[str]
    panel_order: list[str]
    panel_cohort_map: dict[str, str]
    panel_titles: dict[str, str]
    panel_metrics: dict[str, dict[str, str]]
    method_display_names: dict[str, str]
    method_order: list[str]
    default_threshold: float
    thresholds: list[float]
    default_filter_name: str


def load_gene_pipeline_config(config_path: str) -> tuple[GenePipelineConfig, str]:
    cfg_path = Path(config_path).resolve()
    with cfg_path.open("r", encoding="utf-8") as handle:
        p = json.load(handle)

    cohorts = [
        CohortConfig(
            name=name,
            n_case_col=cfg["n_case_col"],
            n_ctrl_col=cfg["n_ctrl_col"],
            stats=cfg["stats"],
            case_total=cfg.get("case_total"),
            ctrl_total=cfg.get("ctrl_total"),
        )
        for name, cfg in p["cohorts"].items()
    ]

    pw_specs: dict[str, PairwiseScoreSpec] = {}
    for agg, spec in p.get("pairwise_score_columns", {}).items():
        pw_specs[agg] = PairwiseScoreSpec(
            vsm_cols=list(spec["vsm"]),
            anchor_full=spec["anchor_full"],
            anchor_pairwise=dict(spec["anchor_pairwise"]),
        )

    return GenePipelineConfig(
        cohorts=cohorts,
        aggregation_types=list(p["aggregation_types"]),
        raw_score_columns={k: list(v) for k, v in p["raw_score_columns"].items()},
        pairwise_score_specs=pw_specs,
        normalized_score_columns=list(p["normalized_score_columns"]),
        panel_order=list(p["panel_order"]),
        panel_cohort_map=dict(p["panel_cohort_map"]),
        panel_titles=dict(p["panel_titles"]),
        panel_metrics=dict(p["panel_metrics"]),
        method_display_names=dict(p["method_display_names"]),
        method_order=list(p["method_order"]),
        default_threshold=float(p.get("default_threshold", 0.95)),
        thresholds=[float(v) for v in p.get("thresholds", [0.90, 0.95, 0.98, 0.99])],
        default_filter_name=str(p.get("default_filter_name", "none")),
    ), str(cfg_path)


# ---------------------------------------------------------------------------
# Column renaming helpers
# ---------------------------------------------------------------------------

def _build_raw_rename_map(agg_score_cols: list[str], normalized_cols: list[str]) -> dict[str, str]:
    """Map aggregation-specific score names to normalized (agg-stripped) names."""
    if len(agg_score_cols) != len(normalized_cols):
        raise ValueError(
            f"raw_score_columns ({len(agg_score_cols)}) and "
            f"normalized_score_columns ({len(normalized_cols)}) must have the same length."
        )
    return dict(zip(agg_score_cols, normalized_cols))


def _build_pairwise_rename_map(spec: PairwiseScoreSpec) -> dict[str, str]:
    """
    Rename gene-level pairwise columns to the variant-level convention
    expected by detect_pairwise_columns():
      {vsm}_{agg}_percentile_with_anchor  →  {vsm}_percentile_with_anchor
      mpc_score_anchor_{agg}_percentile     →  mpc_score_anchor_percentile
      mpc_score_anchor_{agg}_percentile_with_{vsm}  →  mpc_score_anchor_percentile_with_{vsm}
    """
    rename: dict[str, str] = {}

    for vsm_col in spec.vsm_cols:
        base = vsm_col.replace("_mean_percentile_with_anchor", "_percentile_with_anchor")
        base = base.replace("_max_percentile_with_anchor", "_percentile_with_anchor")
        rename[vsm_col] = base

    anchor_norm = spec.anchor_full.replace("_mean_percentile", "_percentile").replace("_max_percentile", "_percentile")
    rename[spec.anchor_full] = anchor_norm

    for vsm_short, anchor_pw_col in spec.anchor_pairwise.items():
        norm = anchor_pw_col.replace("_mean_percentile_with_", "_percentile_with_")
        norm = norm.replace("_max_percentile_with_", "_percentile_with_")
        rename[anchor_pw_col] = norm

    return rename


# ---------------------------------------------------------------------------
# Cohort parquet preparation
# ---------------------------------------------------------------------------

def _prepare_cohort_parquet(
    source_df: pl.DataFrame,
    cohort: CohortConfig,
    score_cols: list[str],
    rename_map: dict[str, str],
    filter_cols: list[str],
    out_path: str,
) -> str:
    """Select relevant columns, rename n_case/n_ctrl and scores, drop null genes."""
    select_exprs: list[pl.Expr] = [
        pl.col("ensg"),
        pl.col(cohort.n_case_col).alias("n_case"),
        pl.col(cohort.n_ctrl_col).alias("n_ctrl"),
    ]
    for orig in score_cols:
        target = rename_map.get(orig, orig)
        select_exprs.append(pl.col(orig).alias(target))
    for fc in filter_cols:
        if fc in source_df.columns:
            select_exprs.append(pl.col(fc))

    cohort_df = source_df.select(select_exprs).drop_nulls(subset=["n_case", "n_ctrl"])
    cohort_df.write_parquet(out_path)
    return out_path


def _write_ephemeral_resources(
    parquet_path: str,
    table_name: str,
    score_cols: list[str],
    filter_map: dict[str, str],
    out_dir: Path,
) -> str:
    resources = {
        "Table_info": {
            table_name: {
                "Path": str(Path(parquet_path).resolve()),
                "Level": "gene",
                "Score_cols": score_cols,
                "Filters": filter_map,
            }
        }
    }
    resources_path = str(out_dir / f"resources_{table_name}.json")
    Path(resources_path).write_text(json.dumps(resources), encoding="utf-8")
    return resources_path


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------

def _run_cohort(
    resources_path: str,
    table_name: str,
    cohort: CohortConfig,
    thresholds: list[float],
    out_prefix: str,
    filters: str | None,
    bootstrap_samples: int | None = None,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> tuple[pl.DataFrame, list[dict[str, Any]]]:
    args = RunArgs(
        resources_json=resources_path,
        table_name=table_name,
        eval_level="gene",
        stat=cohort.stats,
        eval_set=None,
        filters=filters,
        thresholds=",".join(f"{v:g}" for v in thresholds),
        case_total=cohort.case_total,
        ctrl_total=cohort.ctrl_total,
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=bootstrap_samples,
        out_fname=out_prefix,
        write_missing="none",
        pvalue_method=pvalue_method,
    )
    out_df, timings, _ = biostat_run(args)
    return out_df, timings


def compute_gene_metrics(
    config: GenePipelineConfig,
    raw_parquet: str | None,
    pairwise_parquet: str | None,
    mode: PipelineMode,
    outdir: str,
    filters: str | None = None,
    bootstrap_samples: int | None = None,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
) -> dict[str, dict[str, pl.DataFrame]]:
    """
    Run gene-level stats for every (aggregation_type, cohort) combination.

    Returns:
        {agg_type: {"raw": DataFrame, "pairwise": DataFrame}}
        where each DataFrame is the concat of all cohort results, tagged with
        a ``cohort`` column and ``eval_name`` remapped to the cohort name.
    """
    raw_df_full = pl.read_parquet(raw_parquet) if raw_parquet else None
    pw_df_full = pl.read_parquet(pairwise_parquet) if pairwise_parquet else None

    source_df = raw_df_full if raw_df_full is not None else pw_df_full
    filter_cols = [c for c in source_df.columns if c.startswith("filter_")]
    filter_map = {c: c for c in filter_cols}

    results: dict[str, dict[str, pl.DataFrame]] = {}

    with tempfile.TemporaryDirectory(prefix="gene_fig2_") as tmp_str:
        tmp = Path(tmp_str)

        for agg in config.aggregation_types:
            raw_rename = _build_raw_rename_map(
                config.raw_score_columns[agg], config.normalized_score_columns
            )
            pw_rename: dict[str, str] = {}
            if agg in config.pairwise_score_specs:
                pw_rename = _build_pairwise_rename_map(config.pairwise_score_specs[agg])

            all_raw_score_cols = list(raw_rename.values())
            all_pw_score_cols = list(pw_rename.values())

            agg_raw_dfs: list[pl.DataFrame] = []
            agg_pw_dfs: list[pl.DataFrame] = []

            for cohort in config.cohorts:
                cohort_tag = cohort.name
                print(f"  [{agg}] cohort={cohort_tag} stats={cohort.stats}")

                # --- raw ---
                if mode.includes_raw() and raw_df_full is not None:
                    raw_pq = str(tmp / f"raw_{agg}_{cohort_tag}.parquet")
                    _prepare_cohort_parquet(
                        raw_df_full, cohort,
                        score_cols=config.raw_score_columns[agg],
                        rename_map=raw_rename,
                        filter_cols=filter_cols,
                        out_path=raw_pq,
                    )
                    table_name = f"raw_{agg}_{cohort_tag}"
                    res_path = _write_ephemeral_resources(
                        raw_pq, table_name, all_raw_score_cols, filter_map, tmp,
                    )
                    prefix = str(tmp / f"out_raw_{agg}_{cohort_tag}")
                    df, timings = _run_cohort(
                        res_path, table_name, cohort,
                        config.thresholds, prefix, filters,
                        bootstrap_samples=bootstrap_samples,
                        pvalue_method=pvalue_method,
                    )
                    df = df.with_columns(
                        pl.lit(cohort_tag).alias("cohort"),
                        pl.lit(cohort_tag).alias("eval_name"),
                    )
                    agg_raw_dfs.append(df)
                    for t in timings:
                        print(
                            f"    raw  eval={t['eval_name']} filter={t['filter_name']} "
                            f"time_s={t['elapsed_seconds']:.3f}"
                        )

                # --- pairwise ---
                if mode.includes_pairwise() and pw_df_full is not None and agg in config.pairwise_score_specs:
                    pw_score_cols_orig = (
                        config.pairwise_score_specs[agg].vsm_cols
                        + [config.pairwise_score_specs[agg].anchor_full]
                        + list(config.pairwise_score_specs[agg].anchor_pairwise.values())
                    )
                    pw_pq = str(tmp / f"pw_{agg}_{cohort_tag}.parquet")
                    _prepare_cohort_parquet(
                        pw_df_full, cohort,
                        score_cols=pw_score_cols_orig,
                        rename_map=pw_rename,
                        filter_cols=filter_cols,
                        out_path=pw_pq,
                    )
                    table_name = f"pw_{agg}_{cohort_tag}"
                    pw_stats = cohort.stats.replace("enrichment", "pairwise_enrichment").replace(
                        "rate_ratio", "pairwise_rate_ratio"
                    )
                    res_path = _write_ephemeral_resources(
                        pw_pq, table_name, all_pw_score_cols, filter_map, tmp,
                    )
                    prefix = str(tmp / f"out_pw_{agg}_{cohort_tag}")
                    pw_cohort = CohortConfig(
                        name=cohort.name,
                        n_case_col=cohort.n_case_col,
                        n_ctrl_col=cohort.n_ctrl_col,
                        stats=pw_stats,
                        case_total=cohort.case_total,
                        ctrl_total=cohort.ctrl_total,
                    )
                    df, timings = _run_cohort(
                        res_path, table_name, pw_cohort,
                        config.thresholds, prefix, filters,
                        bootstrap_samples=bootstrap_samples,
                        pvalue_method=pvalue_method,
                    )
                    df = df.with_columns(
                        pl.lit(cohort_tag).alias("cohort"),
                        pl.lit(cohort_tag).alias("eval_name"),
                    )
                    agg_pw_dfs.append(df)
                    for t in timings:
                        print(
                            f"    pw   eval={t['eval_name']} filter={t['filter_name']} "
                            f"time_s={t['elapsed_seconds']:.3f}"
                        )

            raw_combined = pl.concat(agg_raw_dfs, how="diagonal_relaxed") if agg_raw_dfs else None
            pw_combined = pl.concat(agg_pw_dfs, how="diagonal_relaxed") if agg_pw_dfs else None
            results[agg] = {"raw": raw_combined, "pairwise": pw_combined}

    return results


# ---------------------------------------------------------------------------
# Panel table
# ---------------------------------------------------------------------------

def build_gene_panel_table(
    df: pl.DataFrame | None,
    metric_family: str,
    config: GenePipelineConfig,
    threshold: float,
) -> pl.DataFrame:
    if df is None or df.is_empty():
        return pl.DataFrame()

    cohort_to_panel = {coh: pid for pid, coh in config.panel_cohort_map.items()}
    method_rank = {name: idx for idx, name in enumerate(config.method_order)}

    expr = (
        pl.col("filter_name").eq(config.default_filter_name)
        & pl.col("threshold").eq(threshold)
    )
    subset = df.filter(expr)
    if subset.is_empty():
        return pl.DataFrame()

    subset = subset.with_columns([
        pl.col("cohort")
        .map_elements(lambda v: cohort_to_panel.get(str(v), "NA"), return_dtype=pl.String)
        .alias("panel_id"),
        pl.col("cohort")
        .map_elements(
            lambda v: config.panel_titles.get(cohort_to_panel.get(str(v), ""), str(v)),
            return_dtype=pl.String,
        )
        .alias("panel_title"),
        pl.lit(metric_family).alias("metric_family"),
        (pl.col("rows_used") / pl.col("total_eval_rows")).alias("rows_used_frac"),
        pl.col("score_name")
        .map_elements(lambda v: config.method_display_names.get(str(v), str(v)), return_dtype=pl.String)
        .alias("method_label"),
        pl.col("score_name")
        .map_elements(lambda v: method_rank.get(str(v), 10_000), return_dtype=pl.Int64)
        .alias("method_rank"),
    ])
    return subset.sort(["metric_family", "panel_id", "stat", "method_rank", "score_name"])


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------

def build_gene_qc_summary(
    panel_df: pl.DataFrame,
    config: GenePipelineConfig,
    mode: PipelineMode,
) -> tuple[pl.DataFrame, list[str]]:
    rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    families = ["raw", "pairwise"] if mode == PipelineMode.BOTH else [mode.value]

    for family in families:
        fam_df = panel_df.filter(pl.col("metric_family") == family)
        for panel_id in config.panel_order:
            panel_stat = config.panel_metrics[panel_id][family]
            count = fam_df.filter(
                (pl.col("panel_id") == panel_id) & (pl.col("stat") == panel_stat)
            ).height
            rows.append({
                "check": "panel_stat_row_count",
                "metric_family": family,
                "panel_id": panel_id,
                "stat": panel_stat,
                "value": float(count),
                "status": "ok" if count > 0 else "warning",
            })
            if count == 0:
                warnings.append(f"No rows for panel `{panel_id}` stat `{panel_stat}` in `{family}`.")

    qc_df = pl.DataFrame(rows) if rows else pl.DataFrame(schema={"check": pl.String, "status": pl.String})
    return qc_df, warnings


def write_gene_qc_report(
    qc_summary: pl.DataFrame,
    warnings: list[str],
    out_path: str,
    *,
    config_path: str,
    outdir: str,
    agg_type: str,
) -> None:
    lines = [
        "# Figure 2 Pipeline QC Report",
        "",
        f"- config: `{config_path}`",
        f"- outdir: `{outdir}`",
        f"- aggregation: `{agg_type}`",
        "",
        "## Warnings",
    ]
    if warnings:
        for w in warnings:
            lines.append(f"- {w}")
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


# ---------------------------------------------------------------------------
# Plotting (thin wrapper around existing plot.py functions)
# ---------------------------------------------------------------------------

class _GeneConfigAdapter:
    """Minimal adapter so existing plot functions can read method_display_names / method_order."""

    def __init__(self, config: GenePipelineConfig) -> None:
        self.method_display_names = config.method_display_names
        self.method_order = config.method_order


def render_gene_plots(
    panel_df: pl.DataFrame,
    config: GenePipelineConfig,
    outdir: str,
    mode: PipelineMode,
    agg_type: str,
) -> list[str]:
    outputs: list[str] = []
    adapter = _GeneConfigAdapter(config)

    if mode.includes_raw():
        png = str(Path(outdir) / f"figure2_raw_{agg_type}.png")
        pdf = str(Path(outdir) / f"figure2_raw_{agg_type}.pdf")
        render_mode_figure(
            panel_df, mode="raw", config=adapter,
            out_png=png, out_pdf=pdf,
            panel_order=config.panel_order,
            panel_titles=config.panel_titles,
            panel_metrics=config.panel_metrics,
            suptitle=f"Gene-level panels — {agg_type} aggregation (raw)",
        )
        outputs.extend([png, pdf])

    if mode.includes_pairwise():
        png = str(Path(outdir) / f"figure2_pairwise_{agg_type}.png")
        pdf = str(Path(outdir) / f"figure2_pairwise_{agg_type}.pdf")
        render_mode_figure(
            panel_df, mode="pairwise", config=adapter,
            out_png=png, out_pdf=pdf,
            panel_order=config.panel_order,
            panel_titles=config.panel_titles,
            panel_metrics=config.panel_metrics,
            suptitle=f"Gene-level panels — {agg_type} aggregation (pairwise)",
        )
        outputs.extend([png, pdf])

    if mode == PipelineMode.BOTH:
        png = str(Path(outdir) / f"figure2_combined_{agg_type}.png")
        pdf = str(Path(outdir) / f"figure2_combined_{agg_type}.pdf")
        render_combined_figure(
            panel_df, config=adapter,
            out_png=png, out_pdf=pdf,
            panel_order=config.panel_order,
            panel_titles=config.panel_titles,
            panel_metrics=config.panel_metrics,
            suptitle=f"Gene-level panels — {agg_type} aggregation (combined)",
        )
        outputs.extend([png, pdf])

    return outputs


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def execute_gene_pipeline(
    config: GenePipelineConfig,
    cfg_path: str,
    mode: PipelineMode,
    outdir: str,
    raw_parquet: str | None,
    pairwise_parquet: str | None,
    filters: str | None = None,
    bootstrap_samples: int | None = None,
    pvalue_method: str = DEFAULT_PVALUE_METHOD,
    render: bool = True,
) -> dict[str, Any]:
    start = time.perf_counter()
    boot_str = f", bootstrap={bootstrap_samples}" if bootstrap_samples else ""
    print(f"Computing gene-level metrics (mode={mode.value}{boot_str}, pvalue={pvalue_method}) ...")

    agg_results = compute_gene_metrics(
        config, raw_parquet, pairwise_parquet, mode, outdir,
        filters=filters, bootstrap_samples=bootstrap_samples,
        pvalue_method=pvalue_method,
    )

    manifest: dict[str, Any] = {
        "config": cfg_path,
        "mode": mode.value,
        "outdir": outdir,
        "aggregation_types": config.aggregation_types,
    }
    all_plot_outputs: list[str] = []

    for agg, family_dfs in agg_results.items():
        agg_dir = str(Path(outdir) / agg)
        Path(agg_dir).mkdir(parents=True, exist_ok=True)

        raw_df = family_dfs.get("raw")
        pw_df = family_dfs.get("pairwise")

        # Write metric TSVs
        if raw_df is not None and not raw_df.is_empty():
            raw_tsv = str(Path(agg_dir) / "metrics_raw.tsv")
            write_tsv(raw_df, raw_tsv)
            manifest[f"{agg}_raw_tsv"] = raw_tsv

        if pw_df is not None and not pw_df.is_empty():
            pw_tsv = str(Path(agg_dir) / "metrics_pairwise.tsv")
            write_tsv(pw_df, pw_tsv)
            manifest[f"{agg}_pairwise_tsv"] = pw_tsv

        # Panel table
        threshold = config.default_threshold
        raw_panel = build_gene_panel_table(raw_df, "raw", config, threshold)
        pw_panel = build_gene_panel_table(pw_df, "pairwise", config, threshold)
        panels = [p for p in [raw_panel, pw_panel] if not p.is_empty()]
        panel_df = pl.concat(panels, how="diagonal_relaxed") if panels else pl.DataFrame()

        panel_path = str(Path(agg_dir) / "panel_table.tsv")
        if not panel_df.is_empty():
            write_tsv(panel_df, panel_path)
        manifest[f"{agg}_panel_table"] = panel_path

        # QC
        qc_df, qc_warnings = build_gene_qc_summary(panel_df, config, mode)
        qc_path = str(Path(agg_dir) / "qc_summary.tsv")
        if not qc_df.is_empty():
            write_tsv(qc_df, qc_path)
        qc_report_path = str(Path(agg_dir) / "qc_report.md")
        write_gene_qc_report(qc_df, qc_warnings, qc_report_path, config_path=cfg_path, outdir=agg_dir, agg_type=agg)

        # Plots
        if render and not panel_df.is_empty():
            plot_outputs = render_gene_plots(panel_df, config, agg_dir, mode, agg)
            all_plot_outputs.extend(plot_outputs)

    manifest["plot_outputs"] = all_plot_outputs
    manifest["elapsed_seconds"] = time.perf_counter() - start
    manifest_path = str(Path(outdir) / "run_manifest.json")
    write_json(manifest, manifest_path)
    print(f"Pipeline outputs written to: {outdir} ({manifest['elapsed_seconds']:.1f}s)")
    return manifest


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Figure 2 pipeline (gene-level evaluation)")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--config", default="figure2_pipeline_config.json")
        sp.add_argument("--mode", choices=["raw", "pairwise", "both"], default="raw")
        sp.add_argument("--outdir", default=None)
        sp.add_argument("--filters", default=None, help="Comma-separated filter names (or omit for none).")
        sp.add_argument(
            "--bootstrap", type=int, default=None, metavar="N",
            help="Enable bootstrap std_error with N samples (e.g., --bootstrap 100).",
        )
        sp.add_argument(
            "--pvalue-method",
            choices=list(PVALUE_METHODS),
            default=DEFAULT_PVALUE_METHOD,
            help=f"P-value calculation method (default: {DEFAULT_PVALUE_METHOD})",
        )
        sp.add_argument("--overwrite", action="store_true")
        sp.add_argument("--dry-run", action="store_true")

    def add_parquet(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--raw-parquet", default=None, help="Gene-level parquet for raw metrics.")
        sp.add_argument("--pairwise-parquet", default=None, help="Gene-level parquet for pairwise metrics.")

    run_p = subparsers.add_parser("run", help="Compute + QC + plot.")
    add_common(run_p)
    add_parquet(run_p)

    comp_p = subparsers.add_parser("compute", help="Compute + QC only.")
    add_common(comp_p)
    add_parquet(comp_p)

    plot_p = subparsers.add_parser("plot", help="Plot from existing panel tables.")
    add_common(plot_p)
    plot_p.add_argument("--panel-table-mean", default=None)
    plot_p.add_argument("--panel-table-max", default=None)

    return parser


def _prepare_outdir(outdir: str, overwrite: bool) -> None:
    out_path = Path(outdir)
    if out_path.exists():
        has_existing = any(out_path.iterdir())
        if has_existing and not overwrite:
            raise ValueError(f"Output dir exists and is not empty: {outdir}. Use --overwrite.")
        if has_existing and overwrite:
            shutil.rmtree(out_path)
    out_path.mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = _build_parser()
    ns = parser.parse_args()

    config, cfg_path = load_gene_pipeline_config(ns.config)
    mode = PipelineMode(ns.mode)

    if ns.command in {"run", "compute"}:
        if mode.includes_raw() and not (ns.raw_parquet and ns.raw_parquet.strip()):
            parser.error("--raw-parquet is required when --mode is raw or both")
        if mode.includes_pairwise() and not (ns.pairwise_parquet and ns.pairwise_parquet.strip()):
            parser.error("--pairwise-parquet is required when --mode is pairwise or both")

    outdir = ns.outdir or str(Path.cwd() / "results" / f"figure2_{mode.value}")
    _prepare_outdir(outdir, ns.overwrite)

    if ns.dry_run:
        print("Dry run summary:")
        print(f"  command={ns.command}")
        print(f"  mode={mode.value}")
        print(f"  config={cfg_path}")
        print(f"  outdir={outdir}")
        print(f"  cohorts={[c.name for c in config.cohorts]}")
        print(f"  aggregation_types={config.aggregation_types}")
        print(f"  thresholds={config.thresholds}")
        if ns.command in {"run", "compute"}:
            print(f"  raw_parquet={getattr(ns, 'raw_parquet', None)}")
            print(f"  pairwise_parquet={getattr(ns, 'pairwise_parquet', None)}")
        return

    if ns.command == "plot":
        for agg in config.aggregation_types:
            panel_table_arg = getattr(ns, f"panel_table_{agg}", None)
            if panel_table_arg is None:
                print(f"  Skipping {agg}: --panel-table-{agg} not provided.")
                continue
            panel_df = pl.read_csv(panel_table_arg, separator="\t")
            agg_dir = str(Path(outdir) / agg)
            Path(agg_dir).mkdir(parents=True, exist_ok=True)
            render_gene_plots(panel_df, config, agg_dir, mode, agg)
        print(f"Plot outputs written to: {outdir}")
        return

    render = ns.command == "run"
    execute_gene_pipeline(
        config, cfg_path, mode, outdir,
        raw_parquet=ns.raw_parquet,
        pairwise_parquet=getattr(ns, "pairwise_parquet", None),
        filters=ns.filters,
        bootstrap_samples=ns.bootstrap,
        pvalue_method=getattr(ns, "pvalue_method", DEFAULT_PVALUE_METHOD),
        render=render,
    )


if __name__ == "__main__":
    main()
