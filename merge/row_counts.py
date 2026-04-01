"""Row-count diagnostics for the VSM table pipeline.

Tracks row counts at every stage where data can be gained or lost, and
writes a markdown report for post-hoc debugging of data coverage issues.
"""
from __future__ import annotations

import os
import posixpath
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime

import polars as pl


_KNOWN_EXTENSIONS = (".tsv.bgz", ".tsv.gz", ".tsv", ".parquet")


def _stem(uri: str) -> str:
    basename = posixpath.basename(uri.rstrip("/"))
    if not basename:
        basename = os.path.basename(uri.rstrip(os.sep))
    lower = basename.lower()
    for ext in _KNOWN_EXTENSIONS:
        if lower.endswith(ext):
            return basename[: len(basename) - len(ext)]
    return os.path.splitext(basename)[0]


def count_parquet_rows(path: str) -> int:
    """Row count from parquet file metadata (near-zero I/O cost)."""
    return pl.scan_parquet(path).select(pl.len()).collect().item()


def count_lazy(lf: pl.LazyFrame) -> int:
    """Materialize only the row count of a LazyFrame."""
    return lf.select(pl.len()).collect().item()


@dataclass
class RowCountsCollector:
    """Accumulates row counts at pipeline checkpoints."""

    config: dict[str, str] = field(default_factory=dict)
    inputs: dict[str, list[tuple[str, int]]] = field(
        default_factory=lambda: defaultdict(list),
    )
    stages: dict[str, dict[str, int]] = field(
        default_factory=lambda: defaultdict(dict),
    )

    def record_config(self, key: str, value: str) -> None:
        self.config[key] = value

    def record_input(self, category: str, source: str, count: int) -> None:
        self.inputs[category].append((source, count))

    def record(self, stage: str, metric: str, count: int) -> None:
        self.stages[stage][metric] = count


def write_report(collector: RowCountsCollector, path: str) -> None:
    """Render *collector* as a markdown row-count report to *path*."""
    lines: list[str] = []

    def _heading(level: int, text: str) -> None:
        lines.append("")
        lines.append(f"{'#' * level} {text}")
        lines.append("")

    def _table(headers: list[str], rows: list[list[str]]) -> None:
        lines.append("| " + " | ".join(headers) + " |")
        lines.append("|" + "|".join(" --- " for _ in headers) + "|")
        for row in rows:
            lines.append("| " + " | ".join(row) + " |")

    def _kv_table(pairs: list[tuple[str, str]]) -> None:
        _table(["Setting", "Value"], [[k, v] for k, v in pairs])

    def _fmt(n: int) -> str:
        return f"{n:,}"

    def _pct(num: int, denom: int) -> str:
        if denom == 0:
            return "N/A"
        return f"{num / denom * 100:.1f}%"

    # --- Title ---
    lines.append("# Row Count Report")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # --- Configuration ---
    _heading(2, "Configuration")
    cfg = collector.config
    config_rows: list[tuple[str, str]] = []
    for key in (
        "output_uri", "join_type", "percentile_order", "reference_score",
        "linker_table", "aggregate_genes", "collapse_genes", "smooth",
        "filter_count",
    ):
        if key in cfg:
            config_rows.append((key, cfg[key]))
    if config_rows:
        _kv_table(config_rows)

    # --- Input Tables ---
    _heading(2, "Input Tables")

    for category, label in [("pred", "Prediction Tables"), ("eval", "Evaluation Tables")]:
        entries = collector.inputs.get(category, [])
        if not entries:
            continue
        _heading(3, label)
        rows = []
        total = 0
        for uri, count in entries:
            rows.append([_stem(uri), _fmt(count)])
            total += count
        rows.append(["**Total**", f"**{_fmt(total)}**"])
        _table(["Source", "Rows"], rows)

    linker_entries = collector.inputs.get("linker", [])
    if linker_entries:
        _heading(3, "Linker Table")
        linker_rows: list[list[str]] = []
        for uri, count in linker_entries:
            linker_rows.append([_stem(uri), _fmt(count)])
        linker_stage = collector.stages.get("linker", {})
        if "unique_variants" in linker_stage:
            linker_rows.append(["Unique variants", _fmt(linker_stage["unique_variants"])])
        if "unique_ensgs" in linker_stage:
            linker_rows.append(["Unique ENSGs", _fmt(linker_stage["unique_ensgs"])])
        _table(["Metric", "Count"], linker_rows)

    # --- Pipeline Stages ---
    _heading(2, "Pipeline Stages")

    eval_stage = collector.stages.get("eval_merge", {})
    if eval_stage:
        _heading(3, "Eval merge (outer join)")
        _table(["Metric", "Count"], [
            ["Merged eval rows", _fmt(eval_stage["rows"])],
        ])

    pred_stage = collector.stages.get("pred_merge", {})
    if pred_stage:
        join_label = cfg.get("join_type", "?")
        _heading(3, f"Pred merge ({join_label} join)")
        pred_rows: list[list[str]] = []
        if "rows" in pred_stage:
            pred_rows.append(["Rows after merge", _fmt(pred_stage["rows"])])
        if "after_null_drop" in pred_stage:
            dropped = pred_stage["rows"] - pred_stage["after_null_drop"]
            pred_rows.append(["Rows dropped (null scores)", _fmt(dropped)])
            pred_rows.append(["Rows after null drop", _fmt(pred_stage["after_null_drop"])])
        _table(["Metric", "Count"], pred_rows)

    linker_stage = collector.stages.get("linker", {})
    if linker_stage:
        _heading(3, "Linker join")
        lk_rows: list[list[str]] = []
        if "pred_after_join" in linker_stage:
            lk_rows.append(["Pred rows after linker join", _fmt(linker_stage["pred_after_join"])])
        if "pred_no_ensg" in linker_stage:
            lk_rows.append(["Pred variants with no ENSG", _fmt(linker_stage["pred_no_ensg"])])
        if lk_rows:
            _table(["Metric", "Count"], lk_rows)

    gene_stage = collector.stages.get("gene_agg", {})
    if gene_stage:
        _heading(3, "Gene aggregation")
        _table(["Metric", "Count"], [
            ["Rows after aggregation", _fmt(gene_stage["rows"])],
        ])

    # --- Coverage Analysis ---
    coverage = collector.stages.get("coverage", {})
    if coverage:
        _heading(2, "Coverage Analysis")
        eval_total = eval_stage.get("rows", 0)
        cov_rows: list[list[str]] = []
        if "eval_in_pred" in coverage:
            cov_rows.append([
                "Eval variants in merged pred",
                _fmt(coverage["eval_in_pred"]),
                _pct(coverage["eval_in_pred"], eval_total),
            ])
        if "eval_not_in_pred" in coverage:
            cov_rows.append([
                "Eval variants NOT in merged pred",
                _fmt(coverage["eval_not_in_pred"]),
                _pct(coverage["eval_not_in_pred"], eval_total),
            ])
        if "eval_in_linker" in coverage:
            cov_rows.append([
                "Eval variants in linker",
                _fmt(coverage["eval_in_linker"]),
                _pct(coverage["eval_in_linker"], eval_total),
            ])
        if "eval_not_in_linker" in coverage:
            cov_rows.append([
                "Eval variants NOT in linker",
                _fmt(coverage["eval_not_in_linker"]),
                _pct(coverage["eval_not_in_linker"], eval_total),
            ])
        if cov_rows:
            _table(["Metric", "Count", "% of eval"], cov_rows)

    # --- Per-Score Null Counts ---
    nulls = collector.stages.get("per_score_nulls", {})
    if nulls:
        eval_total = eval_stage.get("rows", 0)
        _heading(2, "Per-Score Null Counts (after eval-pred join)")
        null_rows: list[list[str]] = []
        for col, n in nulls.items():
            null_rows.append([f"`{col}`", _fmt(n), _pct(n, eval_total)])
        _table(["Column", "Null rows", "% of eval"], null_rows)

    # --- Output ---
    output_stage = collector.stages.get("output", {})
    if output_stage:
        _heading(2, "Output")
        _table(["Metric", "Count"], [
            ["Rows written", _fmt(output_stage["rows"])],
        ])

    # --- Write file ---
    lines.append("")
    report_dir = os.path.dirname(path)
    if report_dir:
        os.makedirs(report_dir, exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(lines))
    print(f"  Row counts written to {path}", file=sys.stderr)
