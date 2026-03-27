"""
Figure 1 Pipeline - One-command workflow for computing metrics and generating figures.

This module provides the main entry point for the figure1-pipeline command.

Usage:
    figure1-pipeline run --config config.json --mode both --raw-parquet ... --pairwise-parquet ...
    figure1-pipeline compute --config config.json --mode both --raw-parquet ... --pairwise-parquet ...
    figure1-pipeline plot --config config.json --panel-table panel_table.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

from biostat_cli.io import write_json
from biostat_cli.pipeline.compute import execute_compute, prepare_outdir
from biostat_cli.pipeline.config import (
    PipelineArgs,
    default_outdir,
    load_pipeline_config,
    resolve_eval_set,
    validate_pipeline_config,
)
from biostat_cli.pipeline.panel import resolve_panel_df_for_plot
from biostat_cli.pipeline.plot import render_plots
from biostat_cli.types import OutputLayout, PipelineMode


def _build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for figure1-pipeline."""
    parser = argparse.ArgumentParser(description="Figure 1 pipeline (compute + QC + plotting)")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add_common_args(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument("--config", default="figure1_pipeline_config.json")
        subparser.add_argument("--mode", choices=["raw", "pairwise", "both"], default="both")
        subparser.add_argument(
            "--profile",
            choices=["paper_figure1", "all_variant"],
            default="paper_figure1",
            help="paper_figure1: fixed 4-panel paper match; all_variant: auto-include all is_pos_* variant evals.",
        )
        subparser.add_argument(
            "--output-layout",
            choices=["combined", "per_eval", "both"],
            default="combined",
            help="combined: single aggregated outputs; per_eval: one-file-per-eval; both: write both.",
        )
        subparser.add_argument("--outdir", default=None)
        subparser.add_argument(
            "--threshold", type=float, default=None, help="Threshold used for panel table and plotting."
        )
        subparser.add_argument("--thresholds", default=None, help="Comma-separated compute thresholds override.")
        subparser.add_argument(
            "--bootstrap",
            type=int,
            default=None,
            help="Enable bootstrap std_error with N samples (e.g., 100).",
        )
        subparser.add_argument(
            "--eval-set",
            default=None,
            help="Optional override: comma-separated eval set to force a custom subset.",
        )
        subparser.add_argument(
            "--paper-strict",
            action="store_true",
            help="Fail if any required paper eval is missing (instead of fallback/warning).",
        )
        subparser.add_argument(
            "--overwrite",
            action="store_true",
            help="Overwrite existing output directory contents.",
        )
        subparser.add_argument("--dry-run", action="store_true")

    def add_parquet_args(subparser: argparse.ArgumentParser) -> None:
        subparser.add_argument(
            "--raw-parquet",
            default=None,
            help="Parquet for raw (all-VSM) metrics. Required when --mode is raw or both.",
        )
        subparser.add_argument(
            "--pairwise-parquet",
            default=None,
            help="Parquet for pairwise MPC-anchor metrics. Required when --mode is pairwise or both.",
        )

    # run subcommand
    run_parser = subparsers.add_parser("run", help="Compute metrics, build panel table/QC, and render plots.")
    add_common_args(run_parser)
    add_parquet_args(run_parser)

    # compute subcommand
    compute_parser = subparsers.add_parser("compute", help="Compute metrics and build panel table/QC only.")
    add_common_args(compute_parser)
    add_parquet_args(compute_parser)

    # plot subcommand
    plot_parser = subparsers.add_parser("plot", help="Render plots from an existing panel table TSV.")
    add_common_args(plot_parser)
    plot_parser.add_argument("--panel-table", default=None, help="Existing panel_table.tsv to plot from.")
    plot_parser.add_argument(
        "--raw-metrics-tsv",
        default=None,
        help="Existing raw biostat metrics TSV (from biostat_cli) to build panel table and plot.",
    )
    plot_parser.add_argument(
        "--pairwise-metrics-tsv",
        default=None,
        help="Existing pairwise biostat metrics TSV (from biostat_cli) to build panel table and plot.",
    )

    return parser


def _enforce_parquet_args(parser: argparse.ArgumentParser, args: PipelineArgs) -> None:
    """Validate that required parquet arguments are provided."""
    if args.command not in {"run", "compute"}:
        return
    if args.mode.includes_raw() and not (args.raw_parquet and args.raw_parquet.strip()):
        parser.error("--raw-parquet is required when --mode is raw or both")
    if args.mode.includes_pairwise() and not (args.pairwise_parquet and args.pairwise_parquet.strip()):
        parser.error("--pairwise-parquet is required when --mode is pairwise or both")


def _handle_dry_run(
    args: PipelineArgs,
    config: any,
    cfg_path: str,
    outdir: str,
    eval_set: list[str],
    panel_layout: any,
    resolve_warnings: list[str],
) -> None:
    """Print dry run summary."""
    print("Dry run summary:")
    print(f"  command={args.command}")
    print(f"  mode={args.mode.value}")
    print(f"  profile={args.profile}")
    print(f"  output_layout={args.output_layout.value}")
    print(f"  config={cfg_path}")
    print(f"  outdir={outdir}")
    print(f"  eval_set={','.join(eval_set)}")
    print(f"  bootstrap={args.bootstrap}")
    print(f"  paper_strict={args.paper_strict}")
    print(f"  overwrite={args.overwrite}")

    if resolve_warnings:
        print("  resolve warnings:")
        for w in resolve_warnings:
            print(f"    - {w}")

    if args.command in {"run", "compute"}:
        print(f"  raw_parquet={args.raw_parquet}")
        print(f"  pairwise_parquet={args.pairwise_parquet}")

        errs, warns = validate_pipeline_config(
            config,
            args.mode,
            raw_parquet=args.raw_parquet,
            pairwise_parquet=args.pairwise_parquet,
            eval_set=eval_set,
            panel_layout=panel_layout,
            paper_strict=args.paper_strict,
        )

        if warns:
            print("  validation warnings:")
            for w in warns:
                print(f"    - {w}")

        if errs:
            print("  validation errors:")
            for e in errs:
                print(f"    - {e}")
            raise SystemExit(2)


def _handle_plot_command(
    args: PipelineArgs,
    config: any,
    cfg_path: str,
    outdir: str,
    eval_set: list[str],
    panel_layout: any,
) -> None:
    """Handle the plot subcommand."""
    prepare_outdir(outdir, overwrite=args.overwrite)

    panel_df, panel_table = resolve_panel_df_for_plot(
        args,
        config=config,
        outdir=outdir,
        eval_set=eval_set,
        panel_layout=panel_layout,
    )

    plot_paths = render_plots(
        panel_df,
        config=config,
        outdir=outdir,
        mode=args.mode,
        output_layout=args.output_layout.value,
        panel_layout=panel_layout,
    )

    write_json(
        {
            "command": "plot",
            "mode": args.mode.value,
            "profile": args.profile,
            "output_layout": args.output_layout.value,
            "config": cfg_path,
            "panel_table": panel_table,
            "outdir": outdir,
            "plot_outputs": plot_paths,
        },
        str(Path(outdir) / "run_manifest.json"),
    )
    print(f"Plot outputs written to: {outdir}")


def _handle_compute_or_run(
    args: PipelineArgs,
    config: any,
    cfg_path: str,
    outdir: str,
    eval_set: list[str],
    panel_layout: any,
) -> None:
    """Handle the compute or run subcommands."""
    prepare_outdir(outdir, overwrite=args.overwrite)

    compute_outputs = execute_compute(
        args,
        config=config,
        outdir=outdir,
        cfg_path=cfg_path,
        eval_set=eval_set,
        panel_layout=panel_layout,
    )

    plot_outputs: list[str] = []
    if args.command == "run":
        panel_df = pl.read_csv(compute_outputs["panel_table_tsv"], separator="\t")
        plot_outputs = render_plots(
            panel_df,
            config=config,
            outdir=outdir,
            mode=args.mode,
            output_layout=args.output_layout.value,
            panel_layout=panel_layout,
        )

    manifest = {
        "command": args.command,
        "mode": args.mode.value,
        "profile": args.profile,
        "output_layout": args.output_layout.value,
        "config": cfg_path,
        "outdir": outdir,
        "outputs": compute_outputs,
        "plot_outputs": plot_outputs,
    }
    write_json(manifest, str(Path(outdir) / "run_manifest.json"))
    print(f"Pipeline outputs written to: {outdir}")


def main() -> None:
    """Main entry point for figure1-pipeline command."""
    parser = _build_parser()
    ns = parser.parse_args()

    args = PipelineArgs(
        command=ns.command,
        config=ns.config,
        mode=PipelineMode(ns.mode),
        profile=ns.profile,
        output_layout=OutputLayout(getattr(ns, "output_layout", "combined").replace("-", "_")),
        outdir=ns.outdir,
        threshold=ns.threshold,
        thresholds=ns.thresholds,
        bootstrap=ns.bootstrap,
        raw_parquet=getattr(ns, "raw_parquet", None),
        pairwise_parquet=getattr(ns, "pairwise_parquet", None),
        panel_table=getattr(ns, "panel_table", None),
        raw_metrics_tsv=getattr(ns, "raw_metrics_tsv", None),
        pairwise_metrics_tsv=getattr(ns, "pairwise_metrics_tsv", None),
        eval_set_override=getattr(ns, "eval_set", None),
        paper_strict=getattr(ns, "paper_strict", False),
        overwrite=ns.overwrite,
        dry_run=ns.dry_run,
    )

    _enforce_parquet_args(parser, args)

    config, cfg_path = load_pipeline_config(args.config)
    outdir = args.outdir or default_outdir(args, cfg_path)

    eval_set, panel_layout, resolve_warnings = resolve_eval_set(args, config)

    if args.command == "plot":
        _handle_plot_command(args, config, cfg_path, outdir, eval_set, panel_layout)
        return

    if args.dry_run:
        _handle_dry_run(args, config, cfg_path, outdir, eval_set, panel_layout, resolve_warnings)
        return

    _handle_compute_or_run(args, config, cfg_path, outdir, eval_set, panel_layout)


__all__ = ["main"]
