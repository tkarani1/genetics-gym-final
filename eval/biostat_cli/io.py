from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import polars as pl


def scan_table(path: str, storage_options: dict[str, Any] | None = None) -> pl.LazyFrame:
    if path.startswith("gs://"):
        return pl.scan_parquet(path, storage_options=storage_options or {})
    return pl.scan_parquet(path)


def collect_lazy(lf: pl.LazyFrame) -> pl.DataFrame:
    # Keep streaming enabled for memory-efficient terminal operations.
    return lf.collect(streaming=True)


def ensure_parent_dir(path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def write_tsv(df: pl.DataFrame, output_path: str) -> None:
    ensure_parent_dir(output_path)
    df.write_csv(output_path, separator="\t")


def write_json(payload: dict[str, Any], output_path: str) -> None:
    ensure_parent_dir(output_path)
    with Path(output_path).open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
