import os
import sys
import hashlib
import tempfile
from typing import Literal

import polars as pl


NULL_VALUES = ["N/A", "N/a", "n/a", "NA", "na", "Na"]


def normalize_chrom_key(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Rename 'chr' -> 'chrom' if present, so all tables use a consistent key."""
    if "chr" in lf.collect_schema().names():
        return lf.rename({"chr": "chrom"})
    return lf


def _is_gcs_uri(uri: str) -> bool:
    return uri.startswith("gs://")


def detect_format(uri: str) -> Literal["parquet", "tsv", "tsv_bgz"]:
    lower = uri.lower()
    if lower.endswith(".parquet"):
        return "parquet"
    if lower.endswith(".tsv.bgz") or lower.endswith(".tsv.gz"):
        return "tsv_bgz"
    if lower.endswith(".tsv"):
        return "tsv"
    raise ValueError(
        f"Cannot determine format for '{uri}'. "
        "Expected extension: .parquet, .tsv, .tsv.bgz, or .tsv.gz"
    )


def _download_gcs_file(uri: str, dest: str) -> None:
    import gcsfs

    fs = gcsfs.GCSFileSystem()
    print(f"  Downloading {uri} ...", file=sys.stderr)
    fs.get(uri, dest)


def _upload_gcs_file(local_path: str, uri: str) -> None:
    import gcsfs

    fs = gcsfs.GCSFileSystem()
    print(f"  Uploading to {uri} ...", file=sys.stderr)
    fs.put(local_path, uri)


def scan_table(uri: str) -> pl.LazyFrame:
    """Return a Polars LazyFrame for any supported format and location."""
    fmt = detect_format(uri)

    if fmt == "parquet":
        return normalize_chrom_key(pl.scan_parquet(uri))

    # TSV or TSV.BGZ -- for GCS we need a local copy first
    path = uri
    if _is_gcs_uri(uri):
        suffix = ".tsv.bgz" if fmt == "tsv_bgz" else ".tsv"
        tmp = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
        tmp.close()
        _download_gcs_file(uri, tmp.name)
        path = tmp.name

    return normalize_chrom_key(pl.scan_csv(
        path,
        separator="\t",
        has_header=True,
        infer_schema_length=100_000,
        ignore_errors=False,
        low_memory=True,
        null_values=NULL_VALUES,
    ))


def _cache_key(uri: str) -> str:
    return hashlib.sha256(uri.encode()).hexdigest()[:16]


def ensure_parquet(
    uri: str,
    cache_dir: str,
    *,
    use_cache: bool = False,
    store_cache: bool = False,
) -> str:
    """
    If *uri* is already parquet, return it unchanged.
    Otherwise convert TSV/TSV.BGZ to parquet.

    When *use_cache* is True, look for an existing conversion in *cache_dir*
    first.  When *store_cache* is True, persist the converted file in
    *cache_dir* for future runs.  If neither flag is set the converted file
    is written to a disposable temp file.
    """
    fmt = detect_format(uri)
    if fmt == "parquet":
        return uri

    parquet_name = f"{_cache_key(uri)}.parquet"
    parquet_path = os.path.join(cache_dir, parquet_name)

    if use_cache and os.path.exists(parquet_path):
        print(f"  Using cached parquet for {uri}", file=sys.stderr)
        return parquet_path

    if store_cache:
        os.makedirs(cache_dir, exist_ok=True)
        dest = parquet_path
    else:
        dest = tempfile.NamedTemporaryFile(
            suffix=".parquet", delete=False,
        ).name

    print(f"  Converting {uri} -> {dest}", file=sys.stderr)
    lf = scan_table(uri)
    lf.sink_parquet(dest, compression="zstd", statistics=True)
    return dest


def write_parquet(lf: pl.LazyFrame, uri: str) -> None:
    """Sink a LazyFrame to a parquet file at *uri* (local or GCS)."""
    if _is_gcs_uri(uri):
        with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            print(f"  Sinking to local temp file ...", file=sys.stderr)
            lf.sink_parquet(tmp_path, compression="zstd", statistics=True)
            _upload_gcs_file(tmp_path, uri)
        finally:
            os.unlink(tmp_path)
    else:
        os.makedirs(os.path.dirname(uri) or ".", exist_ok=True)
        lf.sink_parquet(uri, compression="zstd", statistics=True)
