
## Reads everything_raw_by_chrom_*.tsv.bgz and writes two Parquet linker sets.
## Default: GCS (genetics-gym). Use --local for the original local Download paths.

import argparse
import os
from pathlib import Path

# Dataproc PySpark jobs often omit HOME/USER; Polars needs POLARS_TEMP_DIR (or HOME) for scratch space.
if not os.environ.get("POLARS_TEMP_DIR"):
    os.environ["POLARS_TEMP_DIR"] = os.environ.get("TMPDIR") or "/tmp"

import polars as pl

GCS_BUCKET = "genetics-gym"
INPUT_PREFIX = "linkers/everything_raw_by_chrom_tsv"
OUTPUT_SLIM_PREFIX = "linkers/linker_SNP_ensg_only_parquet"
OUTPUT_ALL_PREFIX = "linkers/linker_SNP_ensg_all_info_parquet"

# Local layout (original paths)
LOCAL_INPUT_DIR = Path("/Users/tk508/Downloads/everything_raw_by_chr/everything_raw_by_chrom_tsv")
LOCAL_OUTPUT_SLIM_DIR = Path("/Users/tk508/Downloads/linker_SNP_ensg_only")
LOCAL_OUTPUT_ALL_DIR = Path("/Users/tk508/Downloads/linker_SNP_ensg_all_info")

# Only read these columns — the TSV is wide; loading all VEP fields blows driver RAM on Dataproc.
SCAN_COLS = ["chrom", "pos", "ref", "alt", "transcript_consequences"]

NULL_VALUES = ["", "NA", "na", "NAN", "NaN", "nan", "NULL", "null"]


def _list_input_files_gcs() -> list[str]:
    import gcsfs

    fs = gcsfs.GCSFileSystem()
    pattern = f"{GCS_BUCKET}/{INPUT_PREFIX}/*.tsv.bgz"
    return sorted(f"gs://{p}" for p in fs.glob(pattern))


def _list_input_files_local() -> list[str]:
    return sorted(str(p) for p in LOCAL_INPUT_DIR.glob("*.tsv.bgz"))


def _join_output_dir(output_dir: str, filename: str) -> str:
    return f"{output_dir.rstrip('/')}/{filename}"


def run_slim(files: list[str], output_dir: str) -> None:
    minimal_dtype = pl.Struct([
        pl.Field("gene_id", pl.String)
    ])

    print(f"Found {len(files)} files. Generating slim version...")

    for file_path in files:
        base = file_path.rsplit("/", 1)[-1]
        output_filename = base.replace(".tsv.bgz", ".slim.parquet")
        output_path = _join_output_dir(output_dir, output_filename)

        print(f"Processing: {base}...")

        lf_slim = (
            pl.scan_csv(
                file_path,
                separator="\t",
                null_values=NULL_VALUES,
                has_header=True,
            )
            .select(SCAN_COLS)
            .with_columns(
                pl.col("transcript_consequences").str.json_decode(dtype=minimal_dtype)
            )
            .unnest("transcript_consequences")
            .with_columns(
                pl.when(pl.col("gene_id").str.starts_with("ENSG"))
                .then(pl.col("gene_id"))
                .otherwise(None)
                .alias("ensg")
            )
            .filter(pl.col("ensg").is_not_null())
            .select([
                "chrom",
                "pos",
                "ref",
                "alt",
                "ensg",
            ])
            .unique()
        )

        lf_slim.sink_parquet(output_path)

    print(f"\nDone! Slim files are in: {output_dir}")


def run_all_info(files: list[str], output_dir: str) -> None:
    consequence_dtype = pl.Struct([
        pl.Field("gene_id", pl.String),
        pl.Field("transcript_id", pl.String),
        pl.Field("mane_select", pl.String),
        pl.Field("canonical", pl.Int64),
    ])

    group_keys = ["chrom", "pos", "ref", "alt", "ensg"]

    print(f"Found {len(files)} files. Starting batch process...")

    for file_path in files:
        base = file_path.rsplit("/", 1)[-1]
        output_filename = base.replace(".tsv.bgz", ".collapsed.parquet")
        output_path = _join_output_dir(output_dir, output_filename)

        print(f"Processing: {base}...")

        lf = (
            pl.scan_csv(
                file_path,
                separator="\t",
                null_values=NULL_VALUES,
                has_header=True,
            )
            .select(SCAN_COLS)
            .with_columns(
                pl.col("transcript_consequences").str.json_decode(dtype=consequence_dtype)
            )
            .unnest("transcript_consequences")
            .with_columns(
                # Keep only Ensembl gene IDs as the grouping key.
                pl.when(pl.col("gene_id").str.starts_with("ENSG"))
                .then(pl.col("gene_id"))
                .otherwise(None)
                .alias("ensg")
            )
            .filter(pl.col("ensg").is_not_null())
            .group_by(group_keys)
            .agg([
                pl.col("gene_id").filter(pl.col("transcript_id").str.starts_with("NM")).unique().sort().alias("refseq_gene"),
                pl.col("gene_id").filter(pl.col("mane_select").is_not_null() & pl.col("gene_id").str.starts_with("ENSG")).unique().sort().alias("mane_ensg"),
                pl.col("gene_id").filter(pl.col("mane_select").is_not_null() & pl.col("transcript_id").str.starts_with("NM")).unique().sort().alias("mane_refseq_gene"),
                pl.col("gene_id").filter((pl.col("canonical") == 1) & pl.col("gene_id").str.starts_with("ENSG")).unique().sort().alias("canonical_ensg"),
                pl.col("gene_id").filter((pl.col("canonical") == 1) & pl.col("transcript_id").str.starts_with("NM")).unique().sort().alias("canonical_refseq_gene"),
                pl.col("transcript_id").filter(pl.col("transcript_id").str.starts_with("ENST")).unique().sort().alias("all_enst"),
                pl.col("transcript_id").filter(pl.col("transcript_id").str.starts_with("NM")).unique().sort().alias("all_NM"),
                pl.col("transcript_id").filter(pl.col("mane_select").is_not_null() & pl.col("transcript_id").str.starts_with("ENST")).unique().sort().alias("mane_enst"),
                pl.col("transcript_id").filter(pl.col("mane_select").is_not_null() & pl.col("transcript_id").str.starts_with("NM")).unique().sort().alias("mane_NM"),
                pl.col("transcript_id").filter((pl.col("canonical") == 1) & pl.col("transcript_id").str.starts_with("ENST")).unique().sort().alias("canonical_enst"),
                pl.col("transcript_id").filter((pl.col("canonical") == 1) & pl.col("transcript_id").str.starts_with("NM")).unique().sort().alias("canonical_NM"),
                pl.col("mane_select").is_not_null().any().alias("has_mane"),
                (pl.col("canonical") == 1).any().alias("has_canonical"),
            ])
        )

        lf.sink_parquet(output_path)

    print(f"\nSuccess! All files are now in: {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build SNP–gene linker Parquet tables (GCS or local).")
    parser.add_argument(
        "--local",
        action="store_true",
        help="Use local paths under ~/Downloads (everything_raw_by_chrom_tsv → linker outputs).",
    )
    args = parser.parse_args()

    if args.local:
        files = _list_input_files_local()
        out_slim = str(LOCAL_OUTPUT_SLIM_DIR)
        out_all = str(LOCAL_OUTPUT_ALL_DIR)
        print(f"Mode: local\n  input: {LOCAL_INPUT_DIR}\n  slim:  {out_slim}\n  all:   {out_all}\n")
    else:
        files = _list_input_files_gcs()
        out_slim = f"gs://{GCS_BUCKET}/{OUTPUT_SLIM_PREFIX}"
        out_all = f"gs://{GCS_BUCKET}/{OUTPUT_ALL_PREFIX}"
        print(f"Mode: GCS\n  input: gs://{GCS_BUCKET}/{INPUT_PREFIX}/\n  slim:  {out_slim}\n  all:   {out_all}\n")

    run_slim(files, out_slim)
  #  run_all_info(files, out_all)


if __name__ == "__main__":
    main()

# df = pl.scan_parquet("gs://genetics-gym/linkers/linker_SNP_ensg_all_info_parquet/*.parquet")
