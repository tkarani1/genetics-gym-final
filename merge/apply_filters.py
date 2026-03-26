#!/usr/bin/env python3
"""
Annotate a table with boolean filter columns derived from filter tables.

Filter tables can be keyed at two levels:

* **Variant-level** -- contain (chrom, pos, ref, alt) columns.
* **Gene-level** -- contain an ``ensg`` column, or are single-column
  headerless files of ENSG identifiers.

For each filter table a boolean column ``filter_{stem}`` is left-joined onto
the reference table: True where the key exists in the filter table, False
otherwise.
"""
import argparse
import os
import posixpath
import sys
import tempfile
import time

import polars as pl

from table_io import ensure_parquet, normalize_chrom_key, write_parquet
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


def _detect_filter_keys(
    filter_lf: pl.LazyFrame,
) -> tuple[pl.LazyFrame, list[str]]:
    """Determine whether *filter_lf* is variant-keyed or gene-keyed.

    Returns the (possibly renamed) LazyFrame and the list of join key
    column names to use.
    """
    schema = filter_lf.collect_schema()
    names = schema.names()

    if all(k in names for k in JOIN_KEYS):
        return filter_lf, JOIN_KEYS

    if "ensg" in names:
        return filter_lf, ["ensg"]

    if len(names) == 1 and names[0].startswith("ENSG"):
        header_val = names[0]
        filter_lf = filter_lf.rename({header_val: "ensg"})
        first_row = pl.LazyFrame({"ensg": [header_val]})
        filter_lf = pl.concat([first_row, filter_lf])
        return filter_lf, ["ensg"]

    raise ValueError(
        f"Cannot determine filter key type. Columns: {names}. "
        f"Expected either variant keys {JOIN_KEYS} or an 'ensg' column."
    )


def apply_filters(
    lf: pl.LazyFrame,
    filter_uris: list[str],
    cache_dir: str,
) -> pl.LazyFrame:
    """
    For each URI in *filter_uris*, left-join a boolean ``filter_{stem}``
    column onto *lf*.

    Automatically detects whether each filter table is variant-keyed
    (chrom, pos, ref, alt) or gene-keyed (ensg) and joins accordingly.
    """
    for uri in filter_uris:
        stem = derive_filter_name(uri)
        flag_col = f"filter_{stem}"

        pq_path = ensure_parquet(uri, cache_dir)
        raw_lf = normalize_chrom_key(pl.scan_parquet(pq_path))
        filter_lf, keys = _detect_filter_keys(raw_lf)

        filter_keys = (
            filter_lf
            .select(keys)
            .unique()
            .with_columns(pl.lit(True).alias(flag_col))
        )

        key_label = "ensg" if keys == ["ensg"] else "variant"
        print(
            f"  Applying filter '{flag_col}' from {uri} ({key_label}-keyed)",
            file=sys.stderr,
        )

        lf = (
            lf.join(filter_keys, on=keys, how="left", coalesce=True)
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
    lf = normalize_chrom_key(pl.scan_parquet(pq_path))

    lf = apply_filters(lf, _parse_uri_list(args.filter_tables), cache_dir)

    write_parquet(lf, args.output)

    elapsed = time.perf_counter() - start
    print(f"Done in {elapsed:.1f}s.", file=sys.stderr)


if __name__ == "__main__":
    main()
