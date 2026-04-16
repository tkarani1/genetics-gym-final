"""Preprocessing module: normalize raw input files into clean staged parquets.

Runs once when data changes (new score version, new eval set, etc.),
producing staged parquets that the DuckDB engine reads via trivial views.
Does not run during interactive UI sessions.

See preprocessing_design.md for the full design rationale.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import sys
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import polars as pl

from ui.registry import AssetRegistry, EvaluationEntry, PredictionEntry

logger = logging.getLogger(__name__)

VARIANT_KEYS = ["chrom", "pos", "ref", "alt"]
GENE_KEY = "ensg"

NULL_VALUES = ["N/A", "N/a", "n/a", "NA", "na", "Na"]


@dataclass
class ManifestEntry:
    source: str
    staged: str
    role: str
    level: str
    columns: list[str]
    source_mtime: str
    rows: int


@dataclass
class Manifest:
    generated: str
    config_hash: str
    files: list[ManifestEntry] = field(default_factory=list)


def _config_hash(config_path: str) -> str:
    with open(config_path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()[:16]


def _source_mtime(path: str) -> str:
    """Get mtime as ISO string. Returns empty string for non-existent/GCS files."""
    if path.startswith("gs://"):
        return ""
    try:
        mtime = os.path.getmtime(path)
        return datetime.fromtimestamp(mtime, tz=timezone.utc).isoformat()
    except OSError:
        return ""


def _detect_format(path: str) -> str:
    lower = path.lower()
    if lower.endswith(".parquet"):
        if os.path.isdir(path):
            return "parquet_dir"
        return "parquet"
    if lower.endswith(".tsv.bgz") or lower.endswith(".tsv.gz"):
        return "tsv_compressed"
    if lower.endswith(".tsv"):
        return "tsv"
    raise ValueError(f"Cannot determine format for: {path}")


def _scan_file(path: str) -> pl.LazyFrame:
    """Scan a file in any supported format, returning a LazyFrame."""
    fmt = _detect_format(path)
    if fmt == "parquet":
        return pl.scan_parquet(path)
    if fmt == "parquet_dir":
        return pl.scan_parquet(os.path.join(path, "**/*.parquet"), allow_missing_columns=True)
    return pl.scan_csv(
        path,
        separator="\t",
        has_header=True,
        infer_schema_length=100_000,
        ignore_errors=False,
        null_values=NULL_VALUES,
    )


def _normalize_keys(
    lf: pl.LazyFrame,
    level: str,
    key_mapping: dict[str, str] | None = None,
) -> pl.LazyFrame:
    """Rename key columns to canonical names and cast pos to Int64."""
    names = lf.collect_schema().names()
    renames: dict[str, str] = {}

    if key_mapping:
        for src, dst in key_mapping.items():
            if src in names:
                renames[src] = dst

    if level == "variant":
        if "chr" in names and "chrom" not in names and "chr" not in renames.values():
            renames["chr"] = "chrom"
    elif level == "gene":
        if "ENSG" in names and "ensg" not in names and "ENSG" not in renames.values():
            renames["ENSG"] = "ensg"

    if renames:
        lf = lf.rename(renames)

    if level == "variant":
        schema = lf.collect_schema()
        if "pos" in schema.names() and schema["pos"] != pl.Int64:
            lf = lf.with_columns(pl.col("pos").cast(pl.Int64))

    return lf


def _derive_staged_name(source_path: str) -> str:
    """Derive a clean filename for the staged parquet from the source path."""
    basename = os.path.basename(source_path.rstrip("/"))
    for ext in (".tsv.bgz", ".tsv.gz", ".tsv", ".parquet"):
        if basename.lower().endswith(ext):
            basename = basename[: len(basename) - len(ext)]
            break
    return f"{basename}.parquet"


def _derive_filter_name(path: str) -> str:
    """Derive a filter name from the path using parent_dir + stem."""
    stripped = path.rstrip("/").rstrip(os.sep)
    parent = os.path.basename(os.path.dirname(stripped))
    stem = _derive_staged_name(path).replace(".parquet", "")
    if parent:
        return f"{parent}_{stem}"
    return stem


class Preprocessor:
    """Normalize raw input files into clean staged parquets.

    Usage::

        prep = Preprocessor("ui/config.json")
        prep.run()
    """

    def __init__(self, config_path: str = "ui/config.json") -> None:
        self.config_path = config_path
        self.registry = AssetRegistry(config_path)
        self.staging_dir = os.path.join(
            os.path.dirname(config_path), "local_cache", "staging"
        )
        self._manifest_path = os.path.join(self.staging_dir, "manifest.json")
        self._existing_manifest = self._load_manifest()

    def _load_manifest(self) -> Manifest | None:
        if not os.path.exists(self._manifest_path):
            return None
        try:
            with open(self._manifest_path, "r") as f:
                raw = json.load(f)
            return Manifest(
                generated=raw["generated"],
                config_hash=raw["config_hash"],
                files=[ManifestEntry(**e) for e in raw["files"]],
            )
        except (json.JSONDecodeError, KeyError):
            logger.warning("Corrupt manifest.json, will re-process all files.")
            return None

    def _save_manifest(self, manifest: Manifest) -> None:
        os.makedirs(os.path.dirname(self._manifest_path), exist_ok=True)
        with open(self._manifest_path, "w") as f:
            json.dump(
                {
                    "generated": manifest.generated,
                    "config_hash": manifest.config_hash,
                    "files": [asdict(e) for e in manifest.files],
                },
                f,
                indent=2,
            )

    def _is_stale(self, source_path: str, role: str) -> bool:
        """Check if a source file needs re-processing."""
        if self._existing_manifest is None:
            return True
        current_mtime = _source_mtime(source_path)
        for entry in self._existing_manifest.files:
            if entry.source == source_path and entry.role == role:
                if entry.source_mtime == current_mtime and current_mtime:
                    return not os.path.exists(entry.staged)
                return True
        return True

    def run(self, force: bool = False) -> None:
        """Process all files declared in config.json.

        Skips files whose source mtime matches the manifest (incremental).
        ``force=True`` reprocesses everything.
        """
        start = time.perf_counter()
        manifest = Manifest(
            generated=datetime.now(timezone.utc).isoformat(),
            config_hash=_config_hash(self.config_path),
            files=[],
        )

        for subdir in ("predictions", "evaluations", "filters"):
            os.makedirs(os.path.join(self.staging_dir, subdir), exist_ok=True)

        for entry in self.registry.list_predictions():
            if not force and not self._is_stale(entry.path, "prediction"):
                existing = self._find_existing_entry(entry.path, "prediction")
                if existing:
                    manifest.files.append(existing)
                    continue
            try:
                me = self.stage_prediction(entry)
                manifest.files.append(me)
            except Exception:
                logger.exception("Failed to stage prediction: %s", entry.path)

        for entry in self.registry.list_evaluations():
            if not force and not self._is_stale(entry.path, "evaluation"):
                existing = self._find_existing_entry(entry.path, "evaluation")
                if existing:
                    manifest.files.append(existing)
                    continue
            try:
                me = self.stage_evaluation(entry)
                manifest.files.append(me)
            except Exception:
                logger.exception("Failed to stage evaluation: %s", entry.path)

        for path in self.registry.list_filters():
            if not force and not self._is_stale(path, "filter"):
                existing = self._find_existing_entry(path, "filter")
                if existing:
                    manifest.files.append(existing)
                    continue
            try:
                me = self.stage_filter(path)
                manifest.files.append(me)
            except Exception:
                logger.exception("Failed to stage filter: %s", path)

        linker = self.registry.linker_table
        if linker:
            if force or self._is_stale(linker, "linker"):
                try:
                    me = self.stage_linker(linker)
                    manifest.files.append(me)
                except Exception:
                    logger.exception("Failed to stage linker: %s", linker)
            else:
                existing = self._find_existing_entry(linker, "linker")
                if existing:
                    manifest.files.append(existing)

        self._save_manifest(manifest)
        elapsed = time.perf_counter() - start
        logger.info(
            "Preprocessing complete: %d files staged in %.1fs",
            len(manifest.files),
            elapsed,
        )

    def _find_existing_entry(
        self, source_path: str, role: str
    ) -> ManifestEntry | None:
        if self._existing_manifest is None:
            return None
        for entry in self._existing_manifest.files:
            if entry.source == source_path and entry.role == role:
                return entry
        return None

    def stage_prediction(self, entry: PredictionEntry) -> ManifestEntry:
        """Normalize keys, project declared score_columns only, write clean parquet."""
        logger.info("Staging prediction: %s", entry.path)
        lf = _scan_file(entry.path)
        lf = _normalize_keys(lf, entry.level, entry.key_mapping or None)

        schema = lf.collect_schema()
        available_cols = schema.names()

        if entry.level == "variant":
            keys = VARIANT_KEYS
        else:
            keys = [GENE_KEY]

        for k in keys:
            if k not in available_cols:
                raise ValueError(
                    f"Key column '{k}' missing in {entry.path}. "
                    f"Available: {available_cols}"
                )

        missing_scores = [c for c in entry.score_columns if c not in available_cols]
        if missing_scores:
            raise ValueError(
                f"Declared score columns {missing_scores} not found in {entry.path}. "
                f"Available: {available_cols}"
            )

        undeclared_floats = [
            c for c in available_cols
            if c not in keys
            and c not in entry.score_columns
            and schema[c] == pl.Float64
        ]
        if undeclared_floats:
            logger.warning(
                "Undeclared Float64 columns in %s (likely coalescing categories, "
                "will be dropped): %s",
                entry.path,
                undeclared_floats,
            )

        select_cols = keys + entry.score_columns
        lf = lf.select(select_cols)

        for sc in entry.score_columns:
            if schema[sc] != pl.Float64:
                lf = lf.with_columns(pl.col(sc).cast(pl.Float64))

        staged_name = _derive_staged_name(entry.path)
        staged_path = os.path.join(self.staging_dir, "predictions", staged_name)
        lf.sink_parquet(staged_path, compression="zstd")

        row_count = pl.scan_parquet(staged_path).select(pl.len()).collect().item()

        return ManifestEntry(
            source=entry.path,
            staged=staged_path,
            role="prediction",
            level=entry.level,
            columns=entry.score_columns,
            source_mtime=_source_mtime(entry.path),
            rows=row_count,
        )

    def stage_evaluation(self, entry: EvaluationEntry) -> ManifestEntry:
        """Normalize keys, project keys + label columns, write clean parquet."""
        logger.info("Staging evaluation: %s", entry.path)
        lf = _scan_file(entry.path)
        lf = _normalize_keys(lf, entry.level, entry.key_mapping or None)

        schema = lf.collect_schema()
        available_cols = schema.names()

        if entry.level == "variant":
            keys = VARIANT_KEYS
            for k in keys:
                if k not in available_cols:
                    raise ValueError(
                        f"Key column '{k}' missing in eval {entry.path}. "
                        f"Available: {available_cols}"
                    )

            if "is_pos" in available_cols:
                label_col = "is_pos"
            elif "is_case" in available_cols:
                label_col = "is_case"
            else:
                bool_cols = [
                    c for c in available_cols
                    if c not in keys and schema[c] == pl.Boolean
                ]
                if len(bool_cols) == 1:
                    label_col = bool_cols[0]
                else:
                    raise ValueError(
                        f"Cannot identify label column in {entry.path}. "
                        f"Expected 'is_pos' or 'is_case'. Available: {available_cols}"
                    )

            select_exprs = [pl.col(k) for k in keys]
            if label_col != "is_pos":
                select_exprs.append(pl.col(label_col).alias("is_pos"))
            else:
                select_exprs.append(pl.col("is_pos"))

            if schema[label_col] != pl.Boolean:
                select_exprs[-1] = (
                    pl.col(label_col).cast(pl.Boolean).alias("is_pos")
                )

            lf = lf.select(select_exprs)
            payload_cols = ["is_pos"]
        else:
            if GENE_KEY not in available_cols:
                raise ValueError(
                    f"Key column '{GENE_KEY}' missing in gene eval {entry.path}. "
                    f"Available: {available_cols}"
                )
            non_key = [c for c in available_cols if c != GENE_KEY and c != "gene_symbol"]
            lf = lf.select([GENE_KEY] + non_key)
            payload_cols = non_key

        staged_name = _derive_staged_name(entry.path)
        staged_path = os.path.join(self.staging_dir, "evaluations", staged_name)
        lf.sink_parquet(staged_path, compression="zstd")

        row_count = pl.scan_parquet(staged_path).select(pl.len()).collect().item()

        return ManifestEntry(
            source=entry.path,
            staged=staged_path,
            role="evaluation",
            level=entry.level,
            columns=payload_cols,
            source_mtime=_source_mtime(entry.path),
            rows=row_count,
        )

    def stage_filter(self, path: str) -> ManifestEntry:
        """Normalize keys (or handle headerless gene TSVs), deduplicate, write parquet."""
        logger.info("Staging filter: %s", path)
        lf = _scan_file(path)
        schema = lf.collect_schema()
        names = schema.names()

        if all(k in names for k in VARIANT_KEYS):
            lf = _normalize_keys(lf, "variant")
            lf = lf.select(VARIANT_KEYS).unique()
            level = "variant"
            keys = VARIANT_KEYS
        elif GENE_KEY in names:
            lf = lf.select([GENE_KEY]).unique()
            level = "gene"
            keys = [GENE_KEY]
        elif len(names) == 1 and names[0].startswith("ENSG"):
            header_val = names[0]
            lf = lf.rename({header_val: GENE_KEY})
            first_row = pl.LazyFrame({GENE_KEY: [header_val]})
            lf = pl.concat([first_row, lf]).unique()
            level = "gene"
            keys = [GENE_KEY]
        else:
            logger.warning(
                "Skipping filter %s: unrecognized key columns %s", path, names
            )
            return ManifestEntry(
                source=path,
                staged="",
                role="filter",
                level="unknown",
                columns=[],
                source_mtime=_source_mtime(path),
                rows=0,
            )

        filter_name = _derive_filter_name(path)
        staged_path = os.path.join(
            self.staging_dir, "filters", f"{filter_name}.parquet"
        )
        lf.sink_parquet(staged_path, compression="zstd")

        row_count = pl.scan_parquet(staged_path).select(pl.len()).collect().item()

        return ManifestEntry(
            source=path,
            staged=staged_path,
            role="filter",
            level=level,
            columns=keys,
            source_mtime=_source_mtime(path),
            rows=row_count,
        )

    def stage_linker(self, path: str) -> ManifestEntry:
        """Project keys + ensg, deduplicate, write clean parquet."""
        logger.info("Staging linker: %s", path)
        lf = _scan_file(path)
        lf = _normalize_keys(lf, "variant")

        schema = lf.collect_schema()
        available = schema.names()

        required = VARIANT_KEYS + [GENE_KEY]
        missing = [c for c in required if c not in available]
        if missing:
            raise ValueError(
                f"Linker missing columns {missing}. Available: {available}"
            )

        lf = lf.select(required).unique()
        staged_path = os.path.join(self.staging_dir, "linker.parquet")
        lf.sink_parquet(staged_path, compression="zstd")

        row_count = pl.scan_parquet(staged_path).select(pl.len()).collect().item()

        return ManifestEntry(
            source=path,
            staged=staged_path,
            role="linker",
            level="variant",
            columns=required,
            source_mtime=_source_mtime(path),
            rows=row_count,
        )

    def validate(self, entry: ManifestEntry) -> bool:
        """Confirm a staged file conforms to its expected schema."""
        if not entry.staged or not os.path.exists(entry.staged):
            logger.error("Staged file missing: %s", entry.staged)
            return False

        schema = pl.scan_parquet(entry.staged).collect_schema()
        actual_cols = schema.names()

        if entry.level == "variant" and entry.role != "filter":
            for k in VARIANT_KEYS:
                if k not in actual_cols:
                    logger.error(
                        "Staged file %s missing key column: %s", entry.staged, k
                    )
                    return False
        elif entry.level == "gene":
            if GENE_KEY not in actual_cols:
                logger.error(
                    "Staged file %s missing key column: %s", entry.staged, GENE_KEY
                )
                return False

        return True


def main() -> None:
    """CLI entry point: python -m ui.prepare [--config path] [--force]"""
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    parser = argparse.ArgumentParser(description="Preprocess raw files into clean staged parquets")
    parser.add_argument("--config", default="ui/config.json", help="Path to config.json")
    parser.add_argument("--force", action="store_true", help="Reprocess all files")
    args = parser.parse_args()

    prep = Preprocessor(args.config)
    prep.run(force=args.force)


if __name__ == "__main__":
    main()
