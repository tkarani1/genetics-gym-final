"""Pipeline configuration loading and validation."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import polars as pl

from biostat_cli.config import detect_pairwise_columns
from biostat_cli.types import (
    OutputLayout,
    PanelLayoutConfig,
    PipelineMode,
    ProfileType,
    RateRatioDenominators,
)

FIGURE1_RAW_TABLE_KEY = "FIGURE1_RAW"
FIGURE1_PAIRWISE_TABLE_KEY = "FIGURE1_PAIRWISE"

DEFAULT_THRESHOLD = 0.95
DEFAULT_THRESHOLDS = [0.90, 0.95, 0.98, 0.99]
DEFAULT_FILTER_NAME = "none"


@dataclass(frozen=True)
class PipelineConfig:
    """Configuration loaded from figure1_pipeline_config.json."""

    raw_score_columns: list[str]
    pairwise_score_columns: list[str]
    eval_set: list[str]
    panel_layout: PanelLayoutConfig
    method_display_names: dict[str, str]
    method_order: list[str]
    rate_ratio_denominators: RateRatioDenominators
    default_threshold: float
    thresholds: list[float]
    default_filter_name: str

    @property
    def panel_order(self) -> list[str]:
        return self.panel_layout.panel_order

    @property
    def panel_eval_map(self) -> dict[str, str]:
        return self.panel_layout.panel_eval_map

    @property
    def panel_titles(self) -> dict[str, str]:
        return self.panel_layout.panel_titles

    @property
    def panel_metrics(self) -> dict[str, dict[str, str]]:
        return self.panel_layout.panel_metrics


@dataclass(frozen=True)
class PipelineArgs:
    """Command-line arguments for the pipeline."""

    command: str
    config: str
    mode: PipelineMode
    profile: ProfileType
    output_layout: OutputLayout
    outdir: str | None
    threshold: float | None
    thresholds: str | None
    bootstrap: int | None
    raw_parquet: str | None
    pairwise_parquet: str | None
    panel_table: str | None
    raw_metrics_tsv: str | None
    pairwise_metrics_tsv: str | None
    eval_set_override: str | None
    paper_strict: bool
    overwrite: bool
    dry_run: bool


@dataclass(frozen=True)
class TableRunConfig:
    """Configuration for a single biostat_cli run."""

    resources_json: str
    table_name: str


def load_pipeline_config(config_path: str) -> tuple[PipelineConfig, str]:
    """
    Load pipeline configuration from JSON file.

    Args:
        config_path: Path to the figure1_pipeline_config.json file.

    Returns:
        Tuple of (PipelineConfig, resolved config path).

    Raises:
        ValueError: If required keys are missing from config.
    """
    cfg_path = Path(config_path).resolve()
    with cfg_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    if "raw_score_columns" not in payload or "pairwise_score_columns" not in payload:
        raise ValueError(
            f"Config {config_path} must define `raw_score_columns` and `pairwise_score_columns` "
            "(table paths are passed via --raw-parquet / --pairwise-parquet on run/compute)."
        )

    panel_layout = PanelLayoutConfig(
        panel_order=list(payload["panel_order"]),
        panel_eval_map=dict(payload["panel_eval_map"]),
        panel_titles=dict(payload["panel_titles"]),
        panel_metrics={
            str(panel_id): {"raw": str(v["raw"]), "pairwise": str(v["pairwise"])}
            for panel_id, v in dict(payload["panel_metrics"]).items()
        },
    )

    rate_ratio_data = {
        str(k): {
            kk: float(vv)
            for kk, vv in dict(v).items()
            if vv is not None and kk in {"case_total", "ctrl_total"}
        }
        for k, v in dict(payload["rate_ratio_denominators"]).items()
    }

    config = PipelineConfig(
        raw_score_columns=list(payload["raw_score_columns"]),
        pairwise_score_columns=list(payload["pairwise_score_columns"]),
        eval_set=list(payload["eval_set"]),
        panel_layout=panel_layout,
        method_display_names={str(k): str(v) for k, v in dict(payload["method_display_names"]).items()},
        method_order=list(payload["method_order"]),
        rate_ratio_denominators=RateRatioDenominators.from_dict(rate_ratio_data),
        default_threshold=float(payload.get("default_threshold", DEFAULT_THRESHOLD)),
        thresholds=[float(v) for v in payload.get("thresholds", DEFAULT_THRESHOLDS)],
        default_filter_name=str(payload.get("default_filter_name", DEFAULT_FILTER_NAME)),
    )
    return config, str(cfg_path)


def validate_pipeline_config(
    config: PipelineConfig,
    mode: PipelineMode,
    *,
    raw_parquet: str | None,
    pairwise_parquet: str | None,
    eval_set: list[str],
    panel_layout: PanelLayoutConfig,
    paper_strict: bool = False,
) -> tuple[list[str], list[str]]:
    """
    Validate pipeline configuration against parquet schemas.

    Args:
        config: Pipeline configuration.
        mode: Pipeline mode (raw/pairwise/both).
        raw_parquet: Path to raw parquet file.
        pairwise_parquet: Path to pairwise parquet file.
        eval_set: List of eval column names to check.
        panel_layout: Panel layout configuration.
        paper_strict: If True, fail on missing evals instead of warning.

    Returns:
        Tuple of (errors, warnings) lists.
    """
    errors: list[str] = []
    warnings: list[str] = []

    runs: list[tuple[str, str, list[str]]] = []
    if mode.includes_raw():
        assert raw_parquet is not None
        runs.append(("raw", raw_parquet, config.raw_score_columns))
    if mode.includes_pairwise():
        assert pairwise_parquet is not None
        runs.append(("pairwise", pairwise_parquet, config.pairwise_score_columns))

    for run_name, table_path, score_cols in runs:
        if table_path.startswith("gs://"):
            warnings.append(f"[{run_name}] table path is remote ({table_path}); local schema checks skipped.")
            continue
        parquet_path = Path(table_path)
        if not parquet_path.exists():
            errors.append(f"[{run_name}] parquet path does not exist: {table_path}")
            continue

        cols = pl.scan_parquet(table_path).collect_schema().names()
        missing_evals = [e for e in eval_set if e not in cols]
        if missing_evals:
            if paper_strict:
                errors.append(f"[{run_name}] missing eval columns (--paper-strict): {missing_evals}")
            else:
                warnings.append(f"[{run_name}] missing eval columns (will skip): {missing_evals}")
        missing_scores = [s for s in score_cols if s not in cols]
        if missing_scores:
            errors.append(f"[{run_name}] missing score columns: {missing_scores}")

        if run_name == "pairwise":
            pairwise_cols = detect_pairwise_columns(cols)
            if pairwise_cols is None:
                errors.append("[pairwise] pairwise score column structure not detected.")

    # Validate panel configuration
    panel_order = panel_layout.panel_order
    for panel_id in panel_order:
        if panel_id not in panel_layout.panel_eval_map:
            errors.append(f"panel `{panel_id}` missing from panel_eval_map.")
        if panel_id not in panel_layout.panel_titles:
            errors.append(f"panel `{panel_id}` missing from panel_titles.")
        if panel_id not in panel_layout.panel_metrics:
            errors.append(f"panel `{panel_id}` missing from panel_metrics.")

    # Check rate ratio denominators for panels using rate_ratio
    rate_ratio_required_evals: set[str] = set()
    for panel_id, eval_name in panel_layout.panel_eval_map.items():
        metrics = panel_layout.panel_metrics.get(panel_id, {})
        raw_metric = str(metrics.get("raw", ""))
        pairwise_metric = str(metrics.get("pairwise", ""))
        if "rate_ratio" in raw_metric or "rate_ratio" in pairwise_metric:
            rate_ratio_required_evals.add(eval_name)

    for eval_name in sorted(rate_ratio_required_evals):
        case_total, ctrl_total = config.rate_ratio_denominators.get_totals_for_eval(eval_name)
        if case_total is None or ctrl_total is None:
            warnings.append(
                f"rate-ratio denominators missing/incomplete for `{eval_name}`. "
                "Rate-ratio values for that eval will be NaN."
            )

    return errors, warnings


def discover_is_pos_evals(parquet_path: str) -> list[str]:
    """
    Discover all is_pos_* columns in a parquet file.

    Args:
        parquet_path: Path to parquet file (local only).

    Returns:
        Sorted list of is_pos_* column names.
    """
    if parquet_path.startswith("gs://"):
        return []
    cols = pl.scan_parquet(parquet_path).collect_schema().names()
    return sorted([c for c in cols if c.startswith("is_pos_")])


def resolve_eval_set(
    args: PipelineArgs,
    config: PipelineConfig,
) -> tuple[list[str], PanelLayoutConfig, list[str]]:
    """
    Resolve eval set based on profile and arguments.

    Args:
        args: Pipeline arguments.
        config: Pipeline configuration.

    Returns:
        Tuple of (eval_set, panel_layout, warnings).
    """
    warnings: list[str] = []

    if args.eval_set_override:
        evals = [e.strip() for e in args.eval_set_override.split(",") if e.strip()]
        panel_eval_map = {e: e for e in evals}
        panel_titles = {e: e.replace("is_pos_", "").replace("_", " ").title() for e in evals}
        default_raw_stat = config.panel_metrics.get(
            list(config.panel_metrics.keys())[0] if config.panel_metrics else "", {}
        ).get("raw", "enrichment")
        default_pairwise_stat = config.panel_metrics.get(
            list(config.panel_metrics.keys())[0] if config.panel_metrics else "", {}
        ).get("pairwise", "pairwise_enrichment")
        panel_metrics = {e: {"raw": default_raw_stat, "pairwise": default_pairwise_stat} for e in evals}
        panel_layout = PanelLayoutConfig(
            panel_order=evals,
            panel_eval_map=panel_eval_map,
            panel_titles=panel_titles,
            panel_metrics=panel_metrics,
        )
        return evals, panel_layout, warnings

    if args.profile == "paper_figure1":
        return config.eval_set, config.panel_layout, warnings

    # all_variant profile: discover evals from parquet
    discovered: list[str] = []
    if args.raw_parquet:
        discovered.extend(discover_is_pos_evals(args.raw_parquet))
    if args.pairwise_parquet:
        for e in discover_is_pos_evals(args.pairwise_parquet):
            if e not in discovered:
                discovered.append(e)
    discovered = sorted(set(discovered))

    if not discovered:
        warnings.append("all_variant profile: no is_pos_* columns found; falling back to config eval_set.")
        return config.eval_set, config.panel_layout, warnings

    default_raw_stat = "enrichment"
    default_pairwise_stat = "pairwise_enrichment"
    if config.panel_metrics:
        first_panel = list(config.panel_metrics.values())[0]
        default_raw_stat = first_panel.get("raw", default_raw_stat)
        default_pairwise_stat = first_panel.get("pairwise", default_pairwise_stat)

    panel_layout = PanelLayoutConfig(
        panel_order=discovered,
        panel_eval_map={e: e for e in discovered},
        panel_titles={e: e.replace("is_pos_", "").replace("_", " ").title() for e in discovered},
        panel_metrics={e: {"raw": default_raw_stat, "pairwise": default_pairwise_stat} for e in discovered},
    )

    return discovered, panel_layout, warnings


def parse_thresholds(raw: str | None) -> list[float] | None:
    """Parse comma-separated thresholds string."""
    if raw is None or not raw.strip():
        return None
    return [float(piece.strip()) for piece in raw.split(",") if piece.strip()]


def default_outdir(args: PipelineArgs, config_path: str) -> str:
    """Generate default output directory name based on arguments."""
    cfg = Path(config_path).stem
    mode = args.mode.value
    cmd = args.command
    threshold = args.threshold if args.threshold is not None else "default"
    thresholds = args.thresholds if args.thresholds else "cfg"
    bootstrap = args.bootstrap if args.bootstrap is not None else "none"
    name = f"figure1_{cfg}_{cmd}_{mode}_thr-{threshold}_thrs-{thresholds}_boot-{bootstrap}"
    safe_name = name.replace(",", "-").replace("/", "-").replace(" ", "").replace(":", "-")
    return str((Path.cwd() / "results" / safe_name).resolve())


__all__ = [
    "FIGURE1_RAW_TABLE_KEY",
    "FIGURE1_PAIRWISE_TABLE_KEY",
    "DEFAULT_THRESHOLD",
    "DEFAULT_THRESHOLDS",
    "DEFAULT_FILTER_NAME",
    "PipelineConfig",
    "PipelineArgs",
    "TableRunConfig",
    "load_pipeline_config",
    "validate_pipeline_config",
    "discover_is_pos_evals",
    "resolve_eval_set",
    "parse_thresholds",
    "default_outdir",
]
