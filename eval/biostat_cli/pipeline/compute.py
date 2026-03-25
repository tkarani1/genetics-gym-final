"""Pipeline compute logic: biostat_cli execution and resource materialization."""

from __future__ import annotations

import shutil
from dataclasses import asdict
from pathlib import Path
from typing import Any

import polars as pl

from biostat_cli import cli as biostat_cli
from biostat_cli.config import get_table_config, load_resources
from biostat_cli.io import write_json, write_tsv
from biostat_cli.pipeline.config import (
    FIGURE1_PAIRWISE_TABLE_KEY,
    FIGURE1_RAW_TABLE_KEY,
    PipelineArgs,
    PipelineConfig,
    TableRunConfig,
    parse_thresholds,
    validate_pipeline_config,
)
from biostat_cli.pipeline.panel import build_panel_table, build_qc_summary, write_qc_report
from biostat_cli.types import PanelLayoutConfig, PipelineMode, RateRatioDenominators


def prepare_outdir(outdir: str, overwrite: bool) -> None:
    """
    Prepare output directory, optionally removing existing contents.

    Args:
        outdir: Path to output directory.
        overwrite: If True, remove existing contents.

    Raises:
        ValueError: If directory exists with contents and overwrite is False.
    """
    out_path = Path(outdir)
    if out_path.exists():
        has_existing = any(out_path.iterdir())
        if has_existing and not overwrite:
            raise ValueError(
                f"Output directory already exists and is not empty: {outdir}. "
                "Use --overwrite to replace existing results."
            )
        if has_existing and overwrite:
            shutil.rmtree(out_path)
    out_path.mkdir(parents=True, exist_ok=True)


def parquet_path_for_resources(path: str) -> str:
    """Normalize parquet path for resources JSON."""
    if path.startswith("gs://"):
        return path
    return str(Path(path).expanduser().resolve())


def materialize_run_resources(
    outdir: str,
    mode: PipelineMode,
    *,
    raw_parquet: str | None,
    pairwise_parquet: str | None,
    raw_score_columns: list[str],
    pairwise_score_columns: list[str],
) -> tuple[str, TableRunConfig | None, TableRunConfig | None]:
    """
    Write ephemeral resources JSON for this pipeline run.

    Args:
        outdir: Output directory.
        mode: Pipeline mode.
        raw_parquet: Path to raw parquet file.
        pairwise_parquet: Path to pairwise parquet file.
        raw_score_columns: Score columns for raw table.
        pairwise_score_columns: Score columns for pairwise table.

    Returns:
        Tuple of (resources_path, raw_config, pairwise_config).
    """
    resources_path = str(Path(outdir) / "figure1_run_resources.json")
    table_info: dict[str, Any] = {}

    if mode.includes_raw():
        assert raw_parquet is not None
        table_info[FIGURE1_RAW_TABLE_KEY] = {
            "Path": parquet_path_for_resources(raw_parquet),
            "Level": "variant",
            "Score_cols": list(raw_score_columns),
        }

    if mode.includes_pairwise():
        assert pairwise_parquet is not None
        table_info[FIGURE1_PAIRWISE_TABLE_KEY] = {
            "Path": parquet_path_for_resources(pairwise_parquet),
            "Level": "variant",
            "Score_cols": list(pairwise_score_columns),
        }

    write_json({"Table_info": table_info}, resources_path)

    raw_cfg = (
        TableRunConfig(resources_json=resources_path, table_name=FIGURE1_RAW_TABLE_KEY)
        if mode.includes_raw()
        else None
    )
    pw_cfg = (
        TableRunConfig(resources_json=resources_path, table_name=FIGURE1_PAIRWISE_TABLE_KEY)
        if mode.includes_pairwise()
        else None
    )

    return resources_path, raw_cfg, pw_cfg


def format_totals_arg(denominators: RateRatioDenominators, key: str) -> str | None:
    """Format per-eval totals as CLI argument string."""
    pieces: list[str] = []
    totals_dict = denominators.case_totals if key == "case_total" else denominators.ctrl_totals
    for eval_name, value in totals_dict.items():
        pieces.append(f"{eval_name}:{value}")
    return ",".join(pieces) if pieces else None


def run_biostat(
    *,
    cfg: TableRunConfig,
    eval_set: list[str],
    thresholds: list[float],
    stat: str,
    out_prefix: str,
    denominators: RateRatioDenominators,
    bootstrap_samples: int | None,
) -> tuple[pl.DataFrame, list[dict[str, Any]], str]:
    """
    Execute biostat_cli.run() with the given configuration.

    Args:
        cfg: Table run configuration.
        eval_set: List of eval column names.
        thresholds: List of threshold values.
        stat: Comma-separated stat names.
        out_prefix: Output file prefix.
        denominators: Rate ratio denominators.
        bootstrap_samples: Number of bootstrap samples (or None).

    Returns:
        Tuple of (output DataFrame, timing info, output TSV path).
    """
    args = biostat_cli.RunArgs(
        resources_json=cfg.resources_json,
        table_name=cfg.table_name,
        eval_level="variant",
        stat=stat,
        eval_set=",".join(eval_set),
        filters="none",
        thresholds=",".join(f"{v:g}" for v in thresholds),
        case_total=None,
        ctrl_total=None,
        case_total_by_eval=format_totals_arg(denominators, "case_total"),
        ctrl_total_by_eval=format_totals_arg(denominators, "ctrl_total"),
        bootstrap_samples=bootstrap_samples,
        out_fname=out_prefix,
        write_missing="none",
    )

    out_df, eval_filter_timings, _ = biostat_cli.run(args)
    out_tsv = f"{out_prefix}.tsv"
    log_json = f"{out_prefix}_log.json"

    write_tsv(out_df, out_tsv)
    write_json(
        {
            "run_args": asdict(args),
            "table_path": get_table_config(load_resources(cfg.resources_json), cfg.table_name).path,
            "output_files": {"tsv": out_tsv, "log": log_json},
            "elapsed_seconds": sum(float(t["elapsed_seconds"]) for t in eval_filter_timings),
            "eval_filter_elapsed_seconds": eval_filter_timings,
        },
        log_json,
    )

    return out_df, eval_filter_timings, out_tsv


def execute_compute(
    args: PipelineArgs,
    config: PipelineConfig,
    outdir: str,
    cfg_path: str,
    eval_set: list[str],
    panel_layout: PanelLayoutConfig,
) -> dict[str, Any]:
    """
    Execute the compute phase of the pipeline.

    Args:
        args: Pipeline arguments.
        config: Pipeline configuration.
        outdir: Output directory.
        cfg_path: Path to config file.
        eval_set: List of eval column names.
        panel_layout: Panel layout configuration.

    Returns:
        Dictionary of output paths and metadata.

    Raises:
        ValueError: If preflight validation fails.
    """
    errors, warnings = validate_pipeline_config(
        config,
        args.mode,
        raw_parquet=args.raw_parquet,
        pairwise_parquet=args.pairwise_parquet,
        eval_set=eval_set,
        panel_layout=panel_layout,
        paper_strict=args.paper_strict,
    )
    if errors:
        raise ValueError("Preflight validation failed:\n- " + "\n- ".join(errors))

    resources_path, raw_run_cfg, pairwise_run_cfg = materialize_run_resources(
        outdir,
        args.mode,
        raw_parquet=args.raw_parquet,
        pairwise_parquet=args.pairwise_parquet,
        raw_score_columns=config.raw_score_columns,
        pairwise_score_columns=config.pairwise_score_columns,
    )

    threshold_for_panels = args.threshold if args.threshold is not None else config.default_threshold
    thresholds_for_compute = parse_thresholds(args.thresholds) or config.thresholds

    raw_df: pl.DataFrame | None = None
    pairwise_df: pl.DataFrame | None = None
    metric_outputs: dict[str, Any] = {}

    # Run combined compute
    if args.output_layout.value in {"combined", "both"}:
        raw_df, pairwise_df, combined_outputs = _run_combined_compute(
            args, config, outdir, eval_set, thresholds_for_compute, raw_run_cfg, pairwise_run_cfg
        )
        metric_outputs.update(combined_outputs)

    # Run per-eval compute
    if args.output_layout.value in {"per_eval", "both"}:
        per_eval_outputs = _run_per_eval_compute(
            args, config, outdir, eval_set, thresholds_for_compute, raw_run_cfg, pairwise_run_cfg
        )
        metric_outputs["per_eval"] = per_eval_outputs

        # If per_eval only, aggregate DataFrames for panel table
        if args.output_layout.value == "per_eval" and raw_df is None and pairwise_df is None:
            raw_df, pairwise_df = _aggregate_per_eval_outputs(args, outdir, eval_set)

    # Build panel table
    panel_df = build_panel_table(
        raw_df=raw_df,
        pairwise_df=pairwise_df,
        config=config,
        threshold=threshold_for_panels,
        eval_set=eval_set,
        panel_layout=panel_layout,
    )
    panel_path = str(Path(outdir) / "panel_table.tsv")
    write_tsv(panel_df, panel_path)

    # Build QC summary
    qc_summary, qc_warnings = build_qc_summary(
        panel_df, config, args.mode, panel_layout.panel_order, panel_layout.panel_metrics
    )
    qc_path = str(Path(outdir) / "qc_summary.tsv")
    write_tsv(qc_summary, qc_path)

    all_warnings = [*warnings, *qc_warnings]

    # Write QC report
    qc_report_path = str(Path(outdir) / "qc_report.md")
    write_qc_report(qc_summary, all_warnings, qc_report_path, config_path=cfg_path, outdir=outdir)

    return {
        "figure1_run_resources_json": resources_path,
        "raw_parquet": args.raw_parquet,
        "pairwise_parquet": args.pairwise_parquet,
        "panel_table_tsv": panel_path,
        "qc_summary_tsv": qc_path,
        "qc_report_md": qc_report_path,
        "warnings": all_warnings,
        "profile": args.profile,
        "output_layout": args.output_layout.value,
        "eval_set": eval_set,
        **metric_outputs,
    }


def _run_combined_compute(
    args: PipelineArgs,
    config: PipelineConfig,
    outdir: str,
    eval_set: list[str],
    thresholds: list[float],
    raw_run_cfg: TableRunConfig | None,
    pairwise_run_cfg: TableRunConfig | None,
) -> tuple[pl.DataFrame | None, pl.DataFrame | None, dict[str, str]]:
    """Run combined metrics computation."""
    raw_df: pl.DataFrame | None = None
    pairwise_df: pl.DataFrame | None = None
    outputs: dict[str, str] = {}

    if args.mode.includes_raw():
        assert raw_run_cfg is not None
        raw_prefix = str(Path(outdir) / "metrics_raw")
        raw_df, _, raw_path = run_biostat(
            cfg=raw_run_cfg,
            eval_set=eval_set,
            thresholds=thresholds,
            stat="enrichment,rate_ratio",
            out_prefix=raw_prefix,
            denominators=config.rate_ratio_denominators,
            bootstrap_samples=args.bootstrap,
        )
        outputs["raw_tsv"] = raw_path

    if args.mode.includes_pairwise():
        assert pairwise_run_cfg is not None
        pw_prefix = str(Path(outdir) / "metrics_pairwise")
        pairwise_df, _, pw_path = run_biostat(
            cfg=pairwise_run_cfg,
            eval_set=eval_set,
            thresholds=thresholds,
            stat="pairwise_enrichment,pairwise_rate_ratio",
            out_prefix=pw_prefix,
            denominators=config.rate_ratio_denominators,
            bootstrap_samples=args.bootstrap,
        )
        outputs["pairwise_tsv"] = pw_path

    return raw_df, pairwise_df, outputs


def _run_per_eval_compute(
    args: PipelineArgs,
    config: PipelineConfig,
    outdir: str,
    eval_set: list[str],
    thresholds: list[float],
    raw_run_cfg: TableRunConfig | None,
    pairwise_run_cfg: TableRunConfig | None,
) -> dict[str, list[str]]:
    """Run per-eval metrics computation."""
    outputs: dict[str, list[str]] = {"raw_tsvs": [], "pairwise_tsvs": []}
    per_eval_dir = Path(outdir) / "per_eval"
    per_eval_dir.mkdir(parents=True, exist_ok=True)

    for eval_name in eval_set:
        safe_eval = eval_name.replace("/", "_").replace(" ", "_")

        if args.mode.includes_raw():
            assert raw_run_cfg is not None
            raw_prefix = str(per_eval_dir / f"metrics_raw_{safe_eval}")
            _, _, raw_path = run_biostat(
                cfg=raw_run_cfg,
                eval_set=[eval_name],
                thresholds=thresholds,
                stat="enrichment,rate_ratio",
                out_prefix=raw_prefix,
                denominators=config.rate_ratio_denominators,
                bootstrap_samples=args.bootstrap,
            )
            outputs["raw_tsvs"].append(raw_path)

        if args.mode.includes_pairwise():
            assert pairwise_run_cfg is not None
            pw_prefix = str(per_eval_dir / f"metrics_pairwise_{safe_eval}")
            _, _, pw_path = run_biostat(
                cfg=pairwise_run_cfg,
                eval_set=[eval_name],
                thresholds=thresholds,
                stat="pairwise_enrichment,pairwise_rate_ratio",
                out_prefix=pw_prefix,
                denominators=config.rate_ratio_denominators,
                bootstrap_samples=args.bootstrap,
            )
            outputs["pairwise_tsvs"].append(pw_path)

    return outputs


def _aggregate_per_eval_outputs(
    args: PipelineArgs,
    outdir: str,
    eval_set: list[str],
) -> tuple[pl.DataFrame | None, pl.DataFrame | None]:
    """Aggregate per-eval outputs into combined DataFrames."""
    all_raw_dfs: list[pl.DataFrame] = []
    all_pw_dfs: list[pl.DataFrame] = []
    per_eval_dir = Path(outdir) / "per_eval"

    for eval_name in eval_set:
        safe_eval = eval_name.replace("/", "_").replace(" ", "_")

        if args.mode.includes_raw():
            raw_path = per_eval_dir / f"metrics_raw_{safe_eval}.tsv"
            if raw_path.exists():
                all_raw_dfs.append(pl.read_csv(str(raw_path), separator="\t"))

        if args.mode.includes_pairwise():
            pw_path = per_eval_dir / f"metrics_pairwise_{safe_eval}.tsv"
            if pw_path.exists():
                all_pw_dfs.append(pl.read_csv(str(pw_path), separator="\t"))

    raw_df = pl.concat(all_raw_dfs, how="diagonal_relaxed") if all_raw_dfs else None
    pairwise_df = pl.concat(all_pw_dfs, how="diagonal_relaxed") if all_pw_dfs else None

    return raw_df, pairwise_df


__all__ = [
    "prepare_outdir",
    "parquet_path_for_resources",
    "materialize_run_resources",
    "format_totals_arg",
    "run_biostat",
    "execute_compute",
]
