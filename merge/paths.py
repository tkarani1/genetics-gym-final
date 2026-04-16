"""
Shared path, URI, and schema utilities for the merge pipeline.

Centralises constants and helpers that were previously duplicated across
create_vsm_table.py, apply_filters.py, add_score.py, and row_counts.py.
"""
from __future__ import annotations

import os
import posixpath
import sys
import tempfile

import polars as pl

from .merge import JOIN_KEYS


KNOWN_EXTENSIONS = (".tsv.bgz", ".tsv.gz", ".tsv", ".parquet")

CACHE_DIR = os.path.join(tempfile.gettempdir(), "vsm_table_cache")


def expand_path(entry: str) -> list[str]:
    """If *entry* is a directory, return all files inside it with a known
    table extension.  Works for both local paths and GCS URIs (the latter
    require a trailing ``/`` to signal directory intent)."""
    if entry.startswith("gs://"):
        if not entry.endswith("/"):
            return [entry]
        import gcsfs

        fs = gcsfs.GCSFileSystem()
        blobs = fs.ls(entry, detail=False)
        found = sorted(
            f"gs://{b}"
            for b in blobs
            if any(b.lower().endswith(ext) for ext in KNOWN_EXTENSIONS)
        )
        if not found:
            print(
                f"  WARNING: GCS directory '{entry}' contains no files with "
                f"recognised extensions {KNOWN_EXTENSIONS}.",
                file=sys.stderr,
            )
        return found

    if os.path.isdir(entry):
        found = sorted(
            os.path.join(entry, name)
            for name in os.listdir(entry)
            if any(name.lower().endswith(ext) for ext in KNOWN_EXTENSIONS)
        )
        if not found:
            print(
                f"  WARNING: Directory '{entry}' contains no files with "
                f"recognised extensions {KNOWN_EXTENSIONS}.",
                file=sys.stderr,
            )
        return found

    return [entry]


def parse_uri_list(raw: str) -> list[str]:
    """Split a comma-separated string of URIs or directories and expand each."""
    entries = [s.strip() for s in raw.split(",") if s.strip()]
    expanded: list[str] = []
    for entry in entries:
        expanded.extend(expand_path(entry))
    return expanded


def derive_stem(uri: str) -> str:
    """Filename without directory or known table extensions."""
    basename = posixpath.basename(uri.rstrip("/"))
    if not basename:
        basename = os.path.basename(uri.rstrip(os.sep))
    lower = basename.lower()
    for ext in KNOWN_EXTENSIONS:
        if lower.endswith(ext):
            return basename[: len(basename) - len(ext)]
    return os.path.splitext(basename)[0]


def derive_filter_name(uri: str) -> str:
    """Extract a filter column stem from a URI or file path.

    Uses ``{parent_dir}_{filename_stem}`` to avoid collisions when
    different filter directories contain identically named files.
    """
    stripped = uri.rstrip("/").rstrip(os.sep)
    basename = posixpath.basename(stripped)
    if not basename:
        basename = os.path.basename(stripped)

    parent = posixpath.basename(posixpath.dirname(stripped))
    if not parent:
        parent = os.path.basename(os.path.dirname(stripped))

    stem = derive_stem(uri)

    if parent:
        return f"{parent}_{stem}"
    return stem


def score_columns(lf: pl.LazyFrame) -> list[str]:
    """Return Float64 non-key columns in *lf* (the actual numeric scores)."""
    schema = lf.collect_schema()
    return [
        c
        for c, dtype in schema.items()
        if c not in JOIN_KEYS and dtype == pl.Float64
    ]
