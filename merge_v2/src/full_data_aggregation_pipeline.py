#!/usr/bin/env python3
"""Run the entire merge_v2 pipeline end-to-end and emit a full audit trail.

This is the single entry point for the ``merge_v2`` data build. Running it once
will (i) download every source table from GCS and (ii) produce every output
Parquet table, by invoking each script in ``merge_v2/src`` in **dependency
order** so that each script's inputs are already on disk by the time it runs.

Execution order (and why)
-------------------------
Each step is a normal CLI invocation of one of the sibling scripts, run as a
subprocess with the *same* Python interpreter that launched this pipeline (so
the children inherit this process's environment / DuckDB install). The order is
the topological sort of the data dependencies:

1. ``download_source_data.py``               -> data/raw_data/{scores,evals,linker,filters}
2. ``create_variant_scores_all_table.py``    -> scores/variant_scores_all_outer.parquet
3. ``create_variant_eval_all_table.py``      -> evals/variant_evals_all.parquet
4. ``create_ensg_eval_all_table.py``         -> evals/ensg_evals_all.parquet
5. ``create_percentile_score_tables.py``     -> scores/variant_scores_{outer_pre,inner_pre,
   (``--variant``)                              inner_post}_percentile.parquet + all_inner
6. ``create_variant_scores_gene_aggregation.py`` -> scores/gene_aggregated/*_ensg[_stats].parquet
   (``--gene-agg-stats``)                        (needs all_outer AND all_inner from step 5)
7. ``create_pairwise_score_tables.py``       -> scores/pairwise/pairwise_{raw,pre,post}/*
   (``--variant``)                              (needs all_outer + outer_pre_percentile)
8. ``create_pairwise_consolidated_score_tables.py`` -> scores/pairwise_consolidated/*
   (``--variant``, optional)                    (collapses each flavor's per-pair files
                                                 into one table; needed only for the
                                                 consolidated analysis route -- see below)
9. ``created_variant_scores_filtered_tables.py`` -> scores/filtered/*_filtered.parquet
   (needs outer_pre_percentile + linker + filters)
10. ``create_analysis_tables.py``            -> full_analysis_tables/** (joins scores x evals)

Pairwise analysis: two routes
------------------------------
There are two ways to turn the pairwise score tables into analysis tables, and
the pipeline supports choosing between them via ``config.json``:

(i)  *individual pairs* -- keep ``create_pairwise_consolidated_score_tables``
     disabled and run ``create_analysis_tables`` with ``--groups pairwise``. This
     writes one analysis table per anchor/non-anchor pair per flavor under
     ``full_analysis_tables/pairwise/<flavor>/``.

(ii) *consolidated* -- enable ``create_pairwise_consolidated_score_tables`` and
     run ``create_analysis_tables`` with ``--groups pairwise_consolidated``. The
     consolidation step first collapses each flavor's per-pair files into a
     single ``scores/pairwise_consolidated/pairwise_<flavor>_consolidated.parquet``
     table (outer-joining the pairs; coalescing the anchor for raw/pre and, for
     post, keeping every per-pair anchor percentile under a
     ``..._<anchor>_x_<non_anchor>`` name), then the analysis step writes one
     analysis table per flavor under
     ``full_analysis_tables/pairwise_consolidated/``.

Either route needs ``create_pairwise_score_tables`` to have produced the per-pair
files first.

Only the variant-level percentile / pairwise flavors are built because the
score manifest declares no ``gene_level`` score sources (so there is no
``ensg_scores_all.parquet`` to percentile or pair). The downstream gene group of
``create_analysis_tables.py`` is fed instead by the *gene-aggregated* variant
tables from step 6.

Configuration (``config.json``)
-------------------------------
Everything a user can tune lives in ``config.json`` next to this script, which
the pipeline reads by default. It has two sections: ``pipeline`` (the run-wide
options mirrored by the CLI flags below) and ``steps`` (the ordered list of
script invocations, each with its ``args``, ``passthrough`` options, and an
``enabled`` toggle). Point at a different file with ``--config PATH`` or ignore
it with ``--no-config``; any CLI flag overrides the matching config value, and
booleans accept ``--flag`` / ``--no-flag`` (e.g. ``--no-skip-download``).

The variant-level obs/exp sources were previously excluded from step 3
(``--exclude _obs_exp_mis``) because they were an upstream preprocessing bug
(they enumerated ~all SNVs across full transcripts, not possible-missense CDS
variants). They have since been corrected -- filtered to missense via an inner
join with the variant->gene linker -- so that guard has been removed and the
obs/exp fields now flow into ``variant_evals_all``.
See ``merge_v2/.agent_reports/obs_exp_mis_preprocessing_assessment.md``.

The audit trail
---------------
Every run writes ``aggregation_pipeline_{timestamp}.json`` (under
``--report-dir``, default ``merge_v2/pipeline_runs``). It is a provenance /
forensics record so that, if anything looks wrong analytically, you can both
(i) replay the exact command that produced any table and (ii) see where a row
count changed (i.e. where data may have been lost or fanned out). The report is
re-written after **every** step, so even a crashed or interrupted run leaves a
usable partial audit. For each step it records:

* the wall-clock ``started_at`` timestamp and how long the step took;
* the exact ``command`` (argv + a copy-pasteable string) used;
* every input file the step consumed -- with its creation time, size, and row
  count -- and every output table it wrote, with its row count.

A ``files_written`` index inverts this into a per-output-file provenance map
(which command/script/inputs produced each table) for direct lookup.

Row counts are read from Parquet footer metadata, so they are cheap and exact
for every Parquet table (a partitioned ``*.parquet`` directory is counted as one
table). Counting the raw bgzipped/TSV *text* sources requires a full
decompressed scan, so it is **off by default**; enable it with
``--count-text-inputs`` (or disable all counting with ``--no-row-counts``).

Run from the ``src/`` directory::

    python full_data_aggregation_pipeline.py                     # uses config.json
    python full_data_aggregation_pipeline.py --config my.json    # alternate config
    python full_data_aggregation_pipeline.py --no-skip-download  # override config bool
    python full_data_aggregation_pipeline.py --memory-limit 20GB --threads 8
    python full_data_aggregation_pipeline.py --from create_percentile_score_tables
    python full_data_aggregation_pipeline.py --only create_analysis_tables
    python full_data_aggregation_pipeline.py --list-steps
    python full_data_aggregation_pipeline.py --dry-run           # print the plan only
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shlex
import socket
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: merge_v2/src/...)
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # merge_v2/
DATA_DIR = PROJECT_DIR / "data"
DATA_CONFIG_DIR = PROJECT_DIR / "data_config"
SCORES_MANIFEST = DATA_CONFIG_DIR / "score_input_data.json"
EVALS_MANIFEST = DATA_CONFIG_DIR / "evals_input_data.json"
DEFAULT_REPORT_DIR = PROJECT_DIR / "pipeline_runs"
DEFAULT_CONFIG_PATH = SRC_DIR / "config.json"

# Extensions treated as a (single) tabular text source for row counting.
TEXT_SUFFIXES = (".tsv", ".tsv.gz", ".tsv.bgz", ".csv", ".csv.gz", ".txt")

# Passthrough flag names a step's underlying script understands.
FLAG_MEMORY = "memory_limit"
FLAG_THREADS = "threads"
FLAG_TEMP = "temp_dir"
FLAG_OVERWRITE = "overwrite"

# Conservative built-in fallback used only when config.json is missing or
# ``--no-config`` is given. config.json next to this file is the user-facing
# source of truth; keep these in sync with it.
DEFAULT_PIPELINE = {
    "skip_download": False,
    "overwrite": False,
    "memory_limit": None,
    "threads": None,
    "temp_dir": None,
    "report_dir": None,
    "count_text_inputs": False,
    "no_row_counts": False,
    "continue_on_error": False,
    "only": None,
    "from_step": None,
    "to_step": None,
    "skip": None,
}

_DUCK = [FLAG_MEMORY, FLAG_THREADS, FLAG_TEMP]
_DUCK_OVERWRITE = _DUCK + [FLAG_OVERWRITE]
DEFAULT_STEPS = [
    {"name": "download_source_data", "script": "download_source_data.py",
     "args": [], "passthrough": [], "enabled": True},
    {"name": "create_variant_scores_all_table", "script": "create_variant_scores_all_table.py",
     "args": [], "passthrough": _DUCK, "enabled": True},
    {"name": "create_variant_eval_all_table", "script": "create_variant_eval_all_table.py",
     "args": [], "passthrough": _DUCK, "enabled": True},
    {"name": "create_ensg_eval_all_table", "script": "create_ensg_eval_all_table.py",
     "args": [], "passthrough": _DUCK, "enabled": True},
    {"name": "create_percentile_score_tables", "script": "create_percentile_score_tables.py",
     "args": ["--variant"], "passthrough": _DUCK, "enabled": True},
    {"name": "create_variant_scores_gene_aggregation", "script": "create_variant_scores_gene_aggregation.py",
     "args": ["--gene-agg-stats"], "passthrough": _DUCK_OVERWRITE, "enabled": True},
    {"name": "create_pairwise_score_tables", "script": "create_pairwise_score_tables.py",
     "args": ["--variant"], "passthrough": _DUCK_OVERWRITE, "enabled": True},
    {"name": "create_pairwise_consolidated_score_tables", "script": "create_pairwise_consolidated_score_tables.py",
     "args": ["--variant"], "passthrough": _DUCK_OVERWRITE, "enabled": False},
    {"name": "created_variant_scores_filtered_tables", "script": "created_variant_scores_filtered_tables.py",
     "args": [], "passthrough": _DUCK_OVERWRITE, "enabled": True},
    {"name": "create_analysis_tables", "script": "create_analysis_tables.py",
     "args": [], "passthrough": _DUCK_OVERWRITE, "enabled": True},
]

# DuckDB is optional here -- only used to read (free) Parquet row counts. If it
# is not importable, the pipeline still runs and just records null row counts.
try:
    import duckdb  # type: ignore

    _DUCKDB_VERSION = duckdb.__version__
except Exception:  # pragma: no cover - environment dependent
    duckdb = None  # type: ignore
    _DUCKDB_VERSION = None


def sql_str(value: str) -> str:
    """Quote a SQL string literal for DuckDB."""
    return "'" + value.replace("'", "''") + "'"


def now_iso() -> str:
    """Current local time as an ISO-8601 string (with timezone offset)."""
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def ts_iso(epoch: float) -> str:
    """Convert an epoch timestamp to a local ISO-8601 string."""
    return datetime.fromtimestamp(epoch, timezone.utc).astimezone().isoformat(timespec="seconds")


def rel(path: Path) -> str:
    """Display a path relative to the project dir when possible, else absolute."""
    try:
        return str(path.resolve().relative_to(PROJECT_DIR))
    except ValueError:
        return str(path)


# ---------------------------------------------------------------------------
# Step definition.
# ---------------------------------------------------------------------------
@dataclass
class Step:
    """One script invocation in the pipeline.

    ``inputs`` / ``outputs`` are *path specs* resolved at run time. A spec is one
    of:

    * ``"tree:<relpath>"``  -- recursively enumerate every table under a dir
      (a ``*.parquet`` file or partitioned dir is one table and is not
      descended into; ``*.tsv*``/``*.csv*`` files are tables too).
    * ``"<relpath with * ? [ ]>"`` -- a glob; each match is a table.
    * ``"<relpath>"`` -- a single table (Parquet file/dir, text file, ...).
    """

    name: str
    script: str
    args: list[str] = field(default_factory=list)
    passthrough: set[str] = field(default_factory=set)
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)

    def script_path(self) -> Path:
        return SRC_DIR / self.script


def read_manifest_paths(manifest: Path, section: str) -> list[str]:
    """Relative ``file_path`` values declared under one manifest section.

    Used so the score/eval merge steps record their *actual* declared source
    files as inputs in the audit. Falls back to an empty list if the manifest is
    missing or malformed (the report then just lists the manifest itself).
    """
    try:
        with manifest.open() as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return []
    return [e["file_path"] for e in data.get(section, []) if isinstance(e, dict) and e.get("file_path")]


def load_config(config_path: Path | None) -> dict:
    """Load the pipeline config, layering it over the built-in defaults.

    Returns ``{"pipeline": {...}, "steps": [...], "source": str|None}``. When
    ``config_path`` is ``None`` (``--no-config``) or the file is missing/invalid,
    the conservative built-in defaults are used and ``source`` is ``None``.
    """
    pipeline = dict(DEFAULT_PIPELINE)
    steps = [dict(s) for s in DEFAULT_STEPS]
    source: str | None = None

    if config_path is not None:
        if not config_path.exists():
            print(f"WARNING: config not found at {rel(config_path)}; using built-in "
                  "defaults (run with --no-config to silence).", file=sys.stderr)
        else:
            try:
                with config_path.open() as fh:
                    data = json.load(fh)
            except (OSError, json.JSONDecodeError) as exc:
                sys.exit(f"ERROR: could not read config {config_path}: {exc}")
            if isinstance(data.get("pipeline"), dict):
                pipeline.update({k: v for k, v in data["pipeline"].items()
                                 if k in DEFAULT_PIPELINE})
            if isinstance(data.get("steps"), list) and data["steps"]:
                steps = data["steps"]
            source = str(config_path)

    return {"pipeline": pipeline, "steps": steps, "source": source}


def extract_exclude_tokens(args: list[str]) -> list[str]:
    """Pull the values passed to a step's ``--exclude`` (nargs='+') option.

    Used so the audit's recorded inputs match what the step actually consumes
    (e.g. the obs/exp sources excluded from the variant eval merge).
    """
    tokens: list[str] = []
    for i, a in enumerate(args):
        if a == "--exclude":
            for nxt in args[i + 1:]:
                if nxt.startswith("-"):
                    break
                tokens.append(nxt)
    return tokens


def step_io_map() -> dict[str, dict[str, list[str]]]:
    """Audit I/O path specs for each known step (manifest-derived, unfiltered).

    Inputs/outputs describe how the audit introspects each script and are keyed
    by step name; the config controls the command/args/order. Manifest sources
    are returned in full -- per-step ``--exclude`` filtering is applied later in
    :func:`build_steps`.
    """
    score_srcs = read_manifest_paths(SCORES_MANIFEST, "variant_level")
    var_eval_srcs = read_manifest_paths(EVALS_MANIFEST, "variant_level")
    gene_eval_srcs = read_manifest_paths(EVALS_MANIFEST, "gene_level")

    p_scores = "data/processed_data/scores"
    p_evals = "data/processed_data/evals"
    linker = "data/raw_data/linker/linker_all.parquet"

    return {
        "download_source_data": {
            "inputs": ["data_config/input_data_locations.json"],
            "outputs": [
                "tree:data/raw_data/scores",
                "tree:data/raw_data/evals",
                "tree:data/raw_data/linker",
                "tree:data/raw_data/filters",
            ],
        },
        "create_variant_scores_all_table": {
            "inputs": [str(rel(SCORES_MANIFEST))] + score_srcs,
            "outputs": [f"{p_scores}/variant_scores_all_outer.parquet"],
        },
        "create_variant_eval_all_table": {
            "inputs": [str(rel(EVALS_MANIFEST))] + var_eval_srcs,
            "outputs": [f"{p_evals}/variant_evals_all.parquet"],
        },
        "create_ensg_eval_all_table": {
            "inputs": [str(rel(EVALS_MANIFEST))] + gene_eval_srcs,
            "outputs": [f"{p_evals}/ensg_evals_all.parquet"],
        },
        "create_percentile_score_tables": {
            "inputs": [f"{p_scores}/variant_scores_all_outer.parquet"],
            "outputs": [
                f"{p_scores}/variant_scores_outer_pre_percentile.parquet",
                f"{p_scores}/variant_scores_inner_pre_percentile.parquet",
                f"{p_scores}/variant_scores_all_inner.parquet",
                f"{p_scores}/variant_scores_inner_post_percentile.parquet",
                f"{p_scores}/variant_scores_percentile_thresholds.tsv",
            ],
        },
        "create_variant_scores_gene_aggregation": {
            "inputs": [
                f"{p_scores}/variant_scores_all_outer.parquet",
                f"{p_scores}/variant_scores_all_inner.parquet",
                linker,
            ],
            "outputs": ["tree:data/processed_data/scores/gene_aggregated"],
        },
        "create_pairwise_score_tables": {
            "inputs": [
                f"{p_scores}/variant_scores_all_outer.parquet",
                f"{p_scores}/variant_scores_outer_pre_percentile.parquet",
            ],
            "outputs": ["tree:data/processed_data/scores/pairwise"],
        },
        "create_pairwise_consolidated_score_tables": {
            "inputs": ["tree:data/processed_data/scores/pairwise"],
            "outputs": ["tree:data/processed_data/scores/pairwise_consolidated"],
        },
        "created_variant_scores_filtered_tables": {
            "inputs": [
                f"{p_scores}/variant_scores_outer_pre_percentile.parquet",
                linker,
                "tree:data/raw_data/filters/variant",
                "tree:data/raw_data/filters/ensg",
            ],
            "outputs": ["tree:data/processed_data/scores/filtered"],
        },
        "create_analysis_tables": {
            "inputs": [
                f"{p_evals}/variant_evals_all.parquet",
                f"{p_evals}/ensg_evals_all.parquet",
                f"{p_scores}/variant_scores_*.parquet",
                "tree:data/processed_data/scores/pairwise",
                "tree:data/processed_data/scores/pairwise_consolidated",
                "data/processed_data/scores/gene_aggregated/*_ensg_stats.parquet",
            ],
            "outputs": ["tree:data/processed_data/full_analysis_tables"],
        },
    }


def build_steps(config_steps: list[dict]) -> list[Step]:
    """Assemble the ordered pipeline steps from the config's ``steps`` list.

    Each config entry supplies the command (``script``/``args``), which
    pipeline options to forward (``passthrough``), and an ``enabled`` toggle.
    The audit I/O specs are looked up by step name from :func:`step_io_map`, and
    any manifest-derived input matching a step's ``--exclude`` tokens is dropped
    so the recorded inputs reflect what the step actually consumes.
    """
    io = step_io_map()
    steps: list[Step] = []
    for entry in config_steps:
        name = entry.get("name")
        if not name:
            sys.exit(f"ERROR: config step missing 'name': {entry!r}")
        if not entry.get("enabled", True):
            continue
        args = list(entry.get("args", []))
        spec = io.get(name, {"inputs": [], "outputs": []})
        excl = extract_exclude_tokens(args)
        inputs = [i for i in spec["inputs"] if not any(tok in i for tok in excl)]
        steps.append(Step(
            name=name,
            script=entry.get("script", f"{name}.py"),
            args=args,
            passthrough=set(entry.get("passthrough", [])),
            inputs=inputs,
            outputs=list(spec["outputs"]),
        ))
    return steps


# ---------------------------------------------------------------------------
# Path-spec expansion + table enumeration.
# ---------------------------------------------------------------------------
def is_parquet_table(path: Path) -> bool:
    """A ``*.parquet`` path is one table whether it is a file or a directory."""
    return path.name.endswith(".parquet")


def is_text_table(path: Path) -> bool:
    return path.is_file() and any(path.name.endswith(s) for s in TEXT_SUFFIXES)


def find_tables(root: Path) -> list[Path]:
    """Recursively list the table files/dirs under ``root`` (sorted).

    A ``*.parquet`` entry (file or partitioned directory) is treated as a single
    table and not descended into; ``*.tsv*``/``*.csv*`` files are tables; any
    other directory is recursed into.
    """
    tables: list[Path] = []
    if not root.exists():
        return tables
    for child in sorted(root.iterdir()):
        if child.name.startswith("."):
            continue  # skip .duckdb_build / .ckpt_* scratch dirs
        if is_parquet_table(child):
            tables.append(child)
        elif is_text_table(child):
            tables.append(child)
        elif child.is_dir():
            tables.extend(find_tables(child))
    return tables


def expand_spec(spec: str) -> list[Path]:
    """Resolve a path spec to concrete table paths (see :class:`Step`)."""
    if spec.startswith("tree:"):
        return find_tables((PROJECT_DIR / spec[len("tree:"):]).resolve())
    if any(ch in spec for ch in "*?[]"):
        return sorted((PROJECT_DIR).glob(spec))
    return [(PROJECT_DIR / spec).resolve()]


# ---------------------------------------------------------------------------
# Row counting + per-table snapshots.
# ---------------------------------------------------------------------------
def count_rows(con, path: Path, count_text: bool) -> tuple[int | None, str]:
    """Return ``(row_count, status)`` for a table path.

    Parquet row counts come from the footer metadata (cheap + exact). Text
    sources require a full scan and are only counted when ``count_text`` is set.
    """
    if con is None:
        return None, "skipped_no_duckdb"
    if not path.exists():
        return None, "missing"
    if is_parquet_table(path):
        pattern = str(path / "**" / "*.parquet") if path.is_dir() else str(path)
        try:
            n = con.execute(
                f"SELECT count(*) FROM read_parquet({sql_str(pattern)})"
            ).fetchone()[0]
            return int(n), "counted"
        except Exception as exc:  # malformed / unreadable parquet
            return None, f"error: {exc}"
    if is_text_table(path):
        if not count_text:
            return None, "skipped_text"
        compression = "gzip" if path.name.endswith((".gz", ".bgz")) else "auto"
        try:
            n = con.execute(
                f"SELECT count(*) FROM read_csv({sql_str(str(path))}, "
                f"all_varchar=true, compression={sql_str(compression)})"
            ).fetchone()[0]
            return int(n), "counted"
        except Exception as exc:
            return None, f"error: {exc}"
    return None, "unsupported"


def dir_size_bytes(path: Path) -> int:
    """Total size in bytes of a file, or of every file under a directory."""
    if path.is_file():
        return path.stat().st_size
    total = 0
    for p in path.rglob("*"):
        if p.is_file():
            try:
                total += p.stat().st_size
            except OSError:
                pass
    return total


def snapshot_table(con, path: Path, count_text: bool) -> dict:
    """Capture path / timestamps / size / row count for one table."""
    info: dict = {"path": rel(path)}
    if not path.exists():
        info.update(exists=False, row_count=None, row_count_status="missing")
        return info
    st = path.stat()
    # st_birthtime is the true creation time on macOS/BSD; fall back to mtime.
    created = getattr(st, "st_birthtime", None)
    n, status = count_rows(con, path, count_text)
    info.update(
        exists=True,
        is_directory=path.is_dir(),
        created_at=ts_iso(created) if created else None,
        modified_at=ts_iso(st.st_mtime),
        size_bytes=dir_size_bytes(path),
        row_count=n,
        row_count_status=status,
    )
    return info


def snapshot_specs(con, specs: list[str], count_text: bool) -> list[dict]:
    """Snapshot every table matched by a list of specs (missing specs noted)."""
    snaps: list[dict] = []
    for spec in specs:
        matches = expand_spec(spec)
        if not matches:
            snaps.append({"path": spec, "exists": False, "row_count": None,
                          "row_count_status": "no_match"})
            continue
        for path in matches:
            snaps.append(snapshot_table(con, path, count_text))
    return snaps


# ---------------------------------------------------------------------------
# Command assembly + execution.
# ---------------------------------------------------------------------------
def build_command(step: Step, args: argparse.Namespace) -> list[str]:
    """Build the subprocess argv for a step, forwarding supported passthroughs."""
    cmd = [sys.executable, str(step.script_path()), *step.args]
    if FLAG_MEMORY in step.passthrough and args.memory_limit:
        cmd += ["--memory-limit", args.memory_limit]
    if FLAG_THREADS in step.passthrough and args.threads:
        cmd += ["--threads", str(args.threads)]
    if FLAG_TEMP in step.passthrough and args.temp_dir:
        cmd += ["--temp-dir", str(args.temp_dir)]
    if FLAG_OVERWRITE in step.passthrough and args.overwrite:
        cmd += ["--overwrite"]
    return cmd


def select_steps(all_steps: list[Step], args: argparse.Namespace) -> list[Step]:
    """Apply --only / --from / --to / --skip / --skip-download selection."""
    names = [s.name for s in all_steps]

    def require(name: str) -> None:
        if name not in names:
            sys.exit(f"ERROR: unknown step '{name}'. Known steps:\n  " + "\n  ".join(names))

    steps = list(all_steps)
    if args.only:
        for n in args.only:
            require(n)
        steps = [s for s in steps if s.name in set(args.only)]
    else:
        if args.from_step:
            require(args.from_step)
            start = names.index(args.from_step)
            steps = [s for s in steps if names.index(s.name) >= start]
        if args.to_step:
            require(args.to_step)
            end = names.index(args.to_step)
            steps = [s for s in steps if names.index(s.name) <= end]

    skip = set(args.skip or [])
    for n in skip:
        require(n)
    if args.skip_download:
        skip.add("download_source_data")
    steps = [s for s in steps if s.name not in skip]

    if not steps:
        sys.exit("ERROR: step selection left nothing to run.")
    return steps


def build_files_index(step_records: list[dict]) -> dict:
    """Invert the step records into a per-output-file provenance index."""
    index: dict = {}
    for rec in step_records:
        provenance = {
            "step": rec["name"],
            "script": rec["script"],
            "command": rec["command_str"],
            "status": rec["status"],
            "written_at": rec.get("finished_at"),
            "duration_seconds": rec.get("duration_seconds"),
            "inputs": [
                {
                    "path": i["path"],
                    "created_at": i.get("created_at"),
                    "modified_at": i.get("modified_at"),
                    "row_count": i.get("row_count"),
                }
                for i in rec.get("inputs", [])
            ],
        }
        for out in rec.get("outputs", []):
            if not out.get("exists"):
                continue
            entry = dict(provenance)
            entry["row_count"] = out.get("row_count")
            entry["size_bytes"] = out.get("size_bytes")
            index[out["path"]] = entry
    return index


def write_report(report_path: Path, report: dict) -> None:
    """Atomically (re)write the audit JSON so partial runs stay readable."""
    report_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = report_path.with_suffix(report_path.suffix + ".tmp")
    tmp.write_text(json.dumps(report, indent=2))
    tmp.replace(report_path)


def main() -> int:
    # Stage 1: figure out which config file to read (so its values can seed the
    # main parser's defaults). --config / --no-config are parsed here and again
    # below (so they show up in --help and are accepted by the full parser).
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    pre.add_argument("--no-config", dest="no_config", action="store_true")
    pre_args, _ = pre.parse_known_args()

    config_path = None if pre_args.no_config else Path(pre_args.config)
    config = load_config(config_path)
    pipe = config["pipeline"]

    # Stage 2: the full parser, with defaults seeded from the config. Any flag
    # the user passes overrides the config value; booleans accept --flag/--no-flag.
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config", default=str(DEFAULT_CONFIG_PATH),
        help=f"Pipeline config JSON to read defaults from (default: {rel(DEFAULT_CONFIG_PATH)}).",
    )
    parser.add_argument(
        "--no-config", dest="no_config", action="store_true",
        help="Ignore config.json and use the built-in defaults.",
    )
    parser.add_argument(
        "--report-dir", type=Path,
        default=pipe["report_dir"] or str(DEFAULT_REPORT_DIR),
        help=f"Directory for the aggregation_pipeline_*.json audit (default: {DEFAULT_REPORT_DIR}).",
    )
    parser.add_argument(
        "--skip-download", action=argparse.BooleanOptionalAction,
        default=bool(pipe["skip_download"]),
        help="Skip the download step (source data already present locally).",
    )
    parser.add_argument(
        "--only", nargs="+", metavar="STEP", default=pipe["only"],
        help="Run only these step(s) by name (mutually exclusive with --from/--to).",
    )
    parser.add_argument(
        "--from", dest="from_step", metavar="STEP", default=pipe["from_step"],
        help="Start at this step (inclusive). Useful to resume after a failure.",
    )
    parser.add_argument(
        "--to", dest="to_step", metavar="STEP", default=pipe["to_step"],
        help="Stop after this step (inclusive).",
    )
    parser.add_argument(
        "--skip", nargs="+", metavar="STEP", default=pipe["skip"],
        help="Skip these step(s) by name.",
    )
    parser.add_argument(
        "--overwrite", action=argparse.BooleanOptionalAction,
        default=bool(pipe["overwrite"]),
        help="Forward --overwrite to steps that skip existing outputs by default "
             "(gene aggregation, pairwise, filtered, analysis). The score/eval/"
             "percentile steps always rewrite their outputs regardless.",
    )
    parser.add_argument(
        "--memory-limit", default=pipe["memory_limit"],
        help="DuckDB memory limit forwarded to every DuckDB step, e.g. '20GB'.",
    )
    parser.add_argument(
        "--threads", type=int, default=pipe["threads"],
        help="DuckDB worker threads forwarded to every DuckDB step.",
    )
    parser.add_argument(
        "--temp-dir", type=Path, default=pipe["temp_dir"],
        help="Spill/build scratch dir forwarded to every DuckDB step.",
    )
    parser.add_argument(
        "--count-text-inputs", action=argparse.BooleanOptionalAction,
        default=bool(pipe["count_text_inputs"]),
        help="Also count rows of raw TSV/bgz text sources (a full scan; slow). "
             "Parquet tables are always counted from footer metadata (cheap).",
    )
    # NB: a single ``--no-row-counts`` BooleanOptionalAction is rejected on
    # Python 3.14 (its option name may not start with ``--no-``), so the two
    # directions are declared explicitly against the same ``no_row_counts`` dest.
    row_counts_group = parser.add_mutually_exclusive_group()
    row_counts_group.add_argument(
        "--no-row-counts", dest="no_row_counts", action="store_true",
        default=bool(pipe["no_row_counts"]),
        help="Do not compute any row counts (fastest auditing).",
    )
    row_counts_group.add_argument(
        "--row-counts", dest="no_row_counts", action="store_false",
        help="Compute row counts (overrides a config no_row_counts=true).",
    )
    parser.add_argument(
        "--continue-on-error", action=argparse.BooleanOptionalAction,
        default=bool(pipe["continue_on_error"]),
        help="Keep running subsequent steps even if one fails (default: stop). "
             "Downstream steps will likely error on missing inputs.",
    )
    parser.add_argument(
        "--list-steps", action="store_true",
        help="Print the ordered step names and exit.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the resolved plan (steps + commands) without running anything.",
    )
    args = parser.parse_args()

    if args.only and (args.from_step or args.to_step):
        sys.exit("ERROR: --only is mutually exclusive with --from/--to.")

    if config["source"]:
        print(f"Config: {rel(Path(config['source']))}")

    all_steps = build_steps(config["steps"])

    if args.list_steps:
        for i, s in enumerate(all_steps, start=1):
            print(f"{i}. {s.name}  ({s.script})")
        return 0

    steps = select_steps(all_steps, args)

    # A short-lived in-memory DuckDB connection just for (cheap) row counts.
    count_text = args.count_text_inputs and not args.no_row_counts
    con = None
    if duckdb is not None and not args.no_row_counts:
        con = duckdb.connect()
        try:
            args.report_dir.mkdir(parents=True, exist_ok=True)
            con.execute(f"SET temp_directory = {sql_str(str(args.report_dir))}")
        except Exception:
            pass

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = args.report_dir / f"aggregation_pipeline_{stamp}.json"

    pipeline_start = time.monotonic()
    report: dict = {
        "schema_version": 1,
        "pipeline": "full_data_aggregation_pipeline",
        "run_id": stamp,
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "duckdb_version": _DUCKDB_VERSION,
        "project_dir": str(PROJECT_DIR),
        "invocation": {
            "argv": sys.argv,
            "command": shlex.join(sys.argv),
            "cwd": os.getcwd(),
            "args": {
                k: (str(v) if isinstance(v, Path) else v)
                for k, v in vars(args).items()
            },
            "config_file": config["source"],
            "config": config,
        },
        "started_at": now_iso(),
        "finished_at": None,
        "duration_seconds": None,
        "status": "running",
        "row_counts_enabled": not args.no_row_counts,
        "text_input_counts_enabled": count_text,
        "planned_steps": [s.name for s in steps],
        "steps": [],
        "files_written": {},
    }

    print(f"merge_v2 full data aggregation pipeline")
    print(f"  python:  {sys.executable}")
    print(f"  steps:   {len(steps)}  ({', '.join(s.name for s in steps)})")
    print(f"  report:  {report_path}\n")

    if args.dry_run:
        print("DRY RUN -- the following commands would be executed in order:\n")
        for i, step in enumerate(steps, start=1):
            cmd = build_command(step, args)
            print(f"[{i}/{len(steps)}] {step.name}")
            print(f"    {shlex.join(cmd)}\n")
        if con is not None:
            con.close()
        return 0

    overall_ok = True
    aborted = False
    for i, step in enumerate(steps, start=1):
        cmd = build_command(step, args)
        record: dict = {
            "name": step.name,
            "script": rel(step.script_path()),
            "command": cmd,
            "command_str": shlex.join(cmd),
            "status": "running",
            "started_at": now_iso(),
            "finished_at": None,
            "duration_seconds": None,
            "return_code": None,
            "inputs": [],
            "outputs": [],
        }
        report["steps"].append(record)

        print(f"{'=' * 78}")
        print(f"[{i}/{len(steps)}] {step.name}")
        print(f"    {record['command_str']}")
        print(f"{'=' * 78}", flush=True)

        # Snapshot inputs *before* the step (captures pre-existing provenance).
        record["inputs"] = snapshot_specs(con, step.inputs, count_text)

        start = time.monotonic()
        try:
            proc = subprocess.run(cmd, cwd=str(SRC_DIR))
            rc = proc.returncode
        except KeyboardInterrupt:
            record["status"] = "interrupted"
            record["finished_at"] = now_iso()
            record["duration_seconds"] = round(time.monotonic() - start, 3)
            report["status"] = "interrupted"
            report["finished_at"] = now_iso()
            report["duration_seconds"] = round(time.monotonic() - pipeline_start, 3)
            write_report(report_path, report)
            print("\nInterrupted. Partial audit written to:", report_path)
            if con is not None:
                con.close()
            return 130

        duration = round(time.monotonic() - start, 3)
        record["return_code"] = rc
        record["finished_at"] = now_iso()
        record["duration_seconds"] = duration
        # Snapshot outputs *after* the step.
        record["outputs"] = snapshot_specs(con, step.outputs, count_text)
        record["status"] = "success" if rc == 0 else "failed"

        # Keep the files index current after each step.
        report["files_written"] = build_files_index(report["steps"])
        write_report(report_path, report)

        out_ct = sum(1 for o in record["outputs"] if o.get("exists"))
        print(f"\n--> {step.name}: {record['status']} in {duration:.1f}s "
              f"(rc={rc}, {out_ct} output table(s))", flush=True)

        if rc != 0:
            overall_ok = False
            if not args.continue_on_error:
                # Mark the remaining selected steps as skipped for a complete audit.
                for skipped in steps[i:]:
                    report["steps"].append({
                        "name": skipped.name,
                        "script": rel(skipped.script_path()),
                        "command": build_command(skipped, args),
                        "command_str": shlex.join(build_command(skipped, args)),
                        "status": "skipped_after_failure",
                        "started_at": None,
                        "finished_at": None,
                        "duration_seconds": None,
                        "return_code": None,
                        "inputs": [],
                        "outputs": [],
                    })
                aborted = True
                break

    report["status"] = (
        "failed" if not overall_ok else "success"
    )
    report["finished_at"] = now_iso()
    report["duration_seconds"] = round(time.monotonic() - pipeline_start, 3)
    report["files_written"] = build_files_index(report["steps"])
    write_report(report_path, report)

    if con is not None:
        con.close()

    print(f"\n{'=' * 78}")
    print(f"Pipeline {report['status'].upper()} in {report['duration_seconds']:.1f}s")
    if aborted:
        print("(stopped at first failure; re-run with --from <step> to resume, "
              "or --continue-on-error to push through)")
    print(f"Audit trail: {report_path}")
    print(f"{'=' * 78}")
    return 0 if overall_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
