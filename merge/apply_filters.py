#!/usr/bin/env python3
"""
Annotate a table with boolean filter columns derived from filter tables.

Each filter table is a curated set of (chrom, pos, ref, alt) variant keys.
For each filter table a boolean column ``filter_{stem}`` is left-joined onto
the reference table: True where the variant key exists in the filter table,
False otherwise.
"""
import argparse
import os
import posixpath
import sys
import tempfile
import time

import polars as pl

from table_io import ensure_parquet, write_parquet
from merge import JOIN_KEYS


_KNOWN_EXTENSIONS = (".tsv.bgz", ".tsv.gz", ".tsv", ".parquet")


def derive_filter_name(uri: str) -> str:
    """
    Extract a filter column stem from a URI or file path.

    Strips the directory component, then removes known compound extensions
    to produce a bare stem suitable for use as ``filter_{stem}``.
    """
    basename = posixpath.basename(uri.rstrip("/"))
    if not basename:
        basename = os.path.basename(uri.rstrip(os.sep))
    lower = basename.lower()
    for ext in _KNOWN_EXTENSIONS:
        if lower.endswith(ext):
            return basename[: len(basename) - len(ext)]
    return os.path.splitext(basename)[0]


def apply_filters(
    lf: pl.LazyFrame,
    filter_uris: list[str],
    cache_dir: str,
) -> pl.LazyFrame:
    """
    For each URI in *filter_uris*, left-join a boolean ``filter_{stem}``
    column onto *lf*.
    """
    for uri in filter_uris:
        stem = derive_filter_name(uri)
        flag_col = f"filter_{stem}"

        pq_path = ensure_parquet(uri, cache_dir)
        storage_opts = (
            {"storage_options": {"token": "google_default"}}
            if pq_path.startswith("gs://")
            else {}
        )

        filter_keys = (
            pl.scan_parquet(pq_path, **storage_opts)
            .select(JOIN_KEYS)
            .unique()
            .with_columns(pl.lit(True).alias(flag_col))
        )

        print(f"  Applying filter '{flag_col}' from {uri}", file=sys.stderr)

        lf = (
            lf.join(filter_keys, on=JOIN_KEYS, how="left", coalesce=True)
            .with_columns(pl.col(flag_col).fill_null(False))
        )

    return lf


def _parse_uri_list(raw: str) -> list[str]:
    return [s.strip() for s in raw.split(",") if s.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate a reference table with boolean filter columns "
            "indicating variant membership in each filter table."
        ),
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="URI (GCS or local) of the reference parquet/tsv file.",
    )
    parser.add_argument(
        "--filter_tables",
        required=True,
        help="Comma-separated URIs of filter tables.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Destination URI for the annotated parquet file.",
    )

    args = parser.parse_args()
    start = time.perf_counter()

    cache_dir = os.path.join(tempfile.gettempdir(), "vsm_table_cache")

    pq_path = ensure_parquet(args.reference, cache_dir)
    storage_opts = (
        {"storage_options": {"token": "google_default"}}
        if pq_path.startswith("gs://")
        else {}
    )
    lf = pl.scan_parquet(pq_path, **storage_opts)

    lf = apply_filters(lf, _parse_uri_list(args.filter_tables), cache_dir)

    write_parquet(lf, args.output)

    elapsed = time.perf_counter() - start
    print(f"Done in {elapsed:.1f}s.", file=sys.stderr)


if __name__ == "__main__":
    main()
