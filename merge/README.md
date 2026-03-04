# VSM Table Builder

Command-line tool for consolidating Variant Scoring Method (VSM) prediction
and evaluation data into a single Parquet table.  The resulting table is
intended for downstream statistical analyses that compare prediction score
columns against evaluation label columns.

## Quick start

```bash
pip install -r requirements.txt

python create_vsm_table.py \
  --prediction_tables "gs://bucket/pred_esm.parquet,gs://bucket/pred_am.parquet" \
  --evaluation_tables "gs://bucket/eval_clinvar.parquet" \
  --output "gs://bucket/output/merged.parquet"
```

## Requirements

* Python 3.10+
* Dependencies (install via `pip install -r requirements.txt`):
  * **polars** -- columnar data processing (Lazy API for streaming large files)
  * **gcsfs** -- Google Cloud Storage filesystem access
  * **fsspec** -- filesystem abstraction used by gcsfs

## Input data

### Prediction tables

Each prediction table must contain the four genomic key columns plus one or
more float-valued score columns:

| chrom | pos | ref | alt | *score_name* |
|-------|-----|-----|-----|-------------|
| chr1  | 100 | A   | T   | 0.92        |

The score column name varies per table and corresponds to the VSM that
produced the score (e.g. `esm_score`, `am_pathogenicity`).

### Evaluation tables

Each evaluation table must contain the four genomic key columns plus one or
more boolean label columns:

| chrom | pos | ref | alt | *is_pos_name* |
|-------|-----|-----|-----|--------------|
| chr1  | 100 | A   | T   | true         |

The label column name varies per table and corresponds to the evaluation
method that assigned the label (e.g. `is_pos_clinvar`).

### Filter tables

Each filter table is a curated list of variants containing only the four
genomic key columns:

| chrom | pos | ref | alt |
|-------|-----|-----|-----|
| chr22 | 2001 | G  | A   |

For each filter table a boolean column `filter_{stem}` is appended to the
output, where `{stem}` is the filename without its extension (e.g.
`ordered_variants_filter.parquet` produces column
`filter_ordered_variants_filter`).  The column is `true` for rows whose
variant key appears in the filter table and `false` otherwise.

### Supported formats

| Format     | Extension(s)          | Notes |
|------------|-----------------------|-------|
| Parquet    | `.parquet`            | Preferred; read directly via Polars Lazy API |
| TSV        | `.tsv`                | Auto-converted to a cached Parquet file before processing |
| TSV (bgzipped) | `.tsv.bgz`, `.tsv.gz` | Same as TSV; bgzip decompression handled transparently |

Both local paths and `gs://` URIs are accepted.  For TSV/TSV.BGZ files
stored in GCS, the file is downloaded to a local temp directory before
conversion.

## CLI reference

```
python create_vsm_table.py [OPTIONS]
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `--prediction_tables` | Comma-separated URIs of prediction tables |
| `--evaluation_tables` | Comma-separated URIs of evaluation tables |
| `--output` | Destination URI for the merged Parquet file |

### Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--join_type` | `inner` | `inner` or `outer` join across all tables |
| `--percentile_order` | `post` | `pre` or `post` -- when to compute percentile ranks (see below) |
| `--filter_tables` | *(none)* | Comma-separated URIs of filter tables -- each adds a boolean column (see below) |

## How merging works

All tables are joined on the composite key `(chrom, pos, ref, alt)`.
Prediction tables and evaluation tables are merged together in a single
pass using a sequential pairwise join.

### Inner join (default)

Only rows whose key appears in **every** input table are retained.  This is
the default and the most conservative option.

### Outer join

All rows from all tables are retained.  Columns from tables that lack a
matching key for a given row will contain null values.

## Percentile rank calculation

For each prediction score column, a companion `{score}_percentile` column
is appended to the output.  The percentile rank uses the **average** method
for resolving ties and is computed as:

```
rank(method="average") / count_of_non_null_values
```

Rows with null/NA scores receive a null percentile value; the row itself is
never dropped.

### Percentile order

The `--percentile_order` flag controls whether percentile ranks are
calculated **before** or **after** the join:

| Mode | Flag value | Behavior |
|------|-----------|----------|
| Pre-merge | `pre` | Percentiles are computed on each prediction table independently, then all tables are inner-joined. Percentile values reflect each table's full distribution. |
| Post-merge | `post` | All tables are joined first, then percentiles are computed on the merged result. Percentile values reflect only the rows that survived the join. |
| Outer join | either | Because no rows are dropped, pre and post produce identical results. The tool always computes post-merge in this case. |

## Filter columns

When `--filter_tables` is provided, the tool left-joins each filter table
onto the merged result (after percentile computation).  Because this is a
left join, no rows are added or removed -- only new boolean columns are
appended.

Filter tables can also be applied standalone via `apply_filters.py`:

```
python apply_filters.py [OPTIONS]
```

| Argument | Description |
|----------|-------------|
| `--reference` | URI of the reference parquet/tsv file to annotate |
| `--filter_tables` | Comma-separated URIs of filter tables |
| `--output` | Destination URI for the annotated Parquet file |

## Examples

**Inner join, percentiles after merge (default):**

```bash
python create_vsm_table.py \
  --prediction_tables "gs://my-bucket/esm.parquet,gs://my-bucket/am.parquet" \
  --evaluation_tables "gs://my-bucket/clinvar.parquet" \
  --output "./merged_inner_post.parquet"
```

**Inner join, percentiles before merge:**

```bash
python create_vsm_table.py \
  --prediction_tables "gs://my-bucket/esm.parquet,gs://my-bucket/am.parquet" \
  --evaluation_tables "gs://my-bucket/clinvar.parquet" \
  --output "./merged_inner_pre.parquet" \
  --join_type inner \
  --percentile_order pre
```

**Outer join (all rows preserved):**

```bash
python create_vsm_table.py \
  --prediction_tables "gs://my-bucket/esm.parquet,gs://my-bucket/am.tsv.bgz" \
  --evaluation_tables "gs://my-bucket/clinvar.parquet,gs://my-bucket/omim.parquet" \
  --output "gs://my-bucket/output/merged_outer.parquet" \
  --join_type outer
```

**With filter tables:**

```bash
python create_vsm_table.py \
  --prediction_tables "gs://my-bucket/esm.parquet,gs://my-bucket/am.parquet" \
  --evaluation_tables "gs://my-bucket/clinvar.parquet" \
  --filter_tables "gs://my-bucket/ordered_variants_filter.parquet,gs://my-bucket/chr22_filter.parquet" \
  --output "./merged_with_filters.parquet"
```

**Standalone filter annotation (no merge):**

```bash
python apply_filters.py \
  --reference "./merged.parquet" \
  --filter_tables "./ordered_variants_filter.parquet,./chr22_filter.parquet" \
  --output "./merged_filtered.parquet"
```

**Mixed formats (TSV + Parquet):**

TSV inputs are automatically converted to cached Parquet files before
processing.  No extra flags are needed.

```bash
python create_vsm_table.py \
  --prediction_tables "./local_scores.tsv,gs://my-bucket/am.parquet" \
  --evaluation_tables "gs://my-bucket/clinvar.parquet" \
  --output "./merged.parquet"
```

## Project layout

```
merge/
├── create_vsm_table.py   # CLI entry point and pipeline orchestration
├── apply_filters.py       # Boolean filter column annotation (also standalone CLI)
├── merge.py               # Multi-table join on genomic keys
├── percentile.py          # Null-safe percentile rank calculation
├── table_io.py            # Format detection, reading, writing (local + GCS)
├── requirements.txt       # Python dependencies
└── README.md              # This file
```

## Performance notes

* All I/O uses the Polars **Lazy API** (`scan_parquet`, `scan_csv`).
  The full computation graph is built lazily and only materialized when
  the final `sink_parquet` writes the output.  This keeps memory usage
  manageable even for tables exceeding 30 GB.
* TSV-to-Parquet conversions are cached in a temp directory
  (`$TMPDIR/vsm_table_cache/`).  Repeated runs with the same TSV inputs
  skip the conversion step.
* Output Parquet files use **zstd** compression.
