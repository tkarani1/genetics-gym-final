# Merge Pipeline

Command-line tools for building consolidated variant-level (VSM) and
gene-level (GSM) tables from prediction scores and evaluation labels.
The resulting Parquet tables are intended for downstream statistical
analyses that compare prediction score columns against evaluation label
columns.

The pipeline supports multiple join strategies, percentile rank
computation, score negation, pairwise anchor comparisons, gene-level
aggregation, spatial smoothing, filter annotation, and a precomputed
fast path for iterating quickly on table combinations.

## Quick start

```bash
pip install -r requirements.txt

# Variant-level inner join (simplest case)
python create_vsm_table.py \
  --prediction_tables "path/to/scores/" \
  --evaluation_tables "path/to/eval1.parquet,path/to/eval2.parquet" \
  --reference_score none \
  --output "merged.parquet"

# Gene-level collapsed table
python create_vsm_table.py \
  --prediction_tables "path/to/scores/" \
  --evaluation_tables "path/to/gene_eval1.parquet,path/to/gene_eval2.parquet" \
  --linker_table "path/to/linker.parquet" \
  --aggregate_genes \
  --join_type inner \
  --reference_score none \
  --output "gsm_inner_collapsed.parquet"
```

## Requirements

* Python 3.9+
* Dependencies (install via `pip install -r requirements.txt`):
  * **polars** -- columnar data processing (Lazy API for streaming)
  * **gcsfs** -- Google Cloud Storage filesystem access
  * **fsspec** -- filesystem abstraction used by gcsfs
  * **pyarrow** -- Parquet metadata inspection and chunked writing

## Input data

### Prediction tables

Each prediction table must contain the four genomic key columns plus one
or more `Float64` score columns:

| chrom | pos | ref | alt | *score_name* |
|-------|-----|-----|-----|-------------|
| chr1  | 100 | A   | T   | 0.92        |

Prediction tables can be provided as individual files, comma-separated
lists, or directories (all files with recognised extensions inside the
directory are loaded).

### Evaluation tables

**Variant-level** evaluation tables contain the four genomic key columns
plus a boolean `is_pos` column:

| chrom | pos | ref | alt | is_pos |
|-------|-----|-----|-----|--------|
| chr1  | 100 | A   | T   | true   |

The `is_pos` column is automatically renamed to `is_pos_{filename_stem}`
in the output.

**Gene-level** evaluation tables are keyed on `ensg` and contain numeric
label columns (e.g. `n_case`, `n_ctrl`). When multiple gene-level eval
tables share column names, they are automatically suffixed with the
filename stem to avoid collisions.

### Linker table

Maps variant keys to gene identifiers. Required for gene-level
aggregation; optional for variant-level tables (adds an `ensg` column
for gene-filter support).

| chrom | pos | ref | alt | ensg |
|-------|-----|-----|-----|------|
| chr1  | 100 | A   | T   | ENSG00000000001 |

A single variant can map to multiple genes (1:many).

### Filter tables

**Variant-level** filter tables contain only the four genomic key columns.
**Gene-level** filter tables contain an `ensg` column (or are single-column
headerless files of ENSG identifiers). The key type is auto-detected.

For each filter table a boolean column `filter_{parent_dir}_{stem}` is
appended to the output (`true` where the key exists, `false` otherwise).
Filters can also be provided as directories, in which case all files
with recognised extensions inside are used.

### Supported formats

| Format         | Extension(s)              | Notes |
|----------------|---------------------------|-------|
| Parquet        | `.parquet`                | Preferred; read directly via Polars Lazy API |
| TSV            | `.tsv`                    | Auto-converted to a cached Parquet before processing |
| TSV (compressed) | `.tsv.bgz`, `.tsv.gz`   | Decompression handled transparently |

Both local paths and `gs://` URIs are accepted.

## Pipeline modes

### Full pipeline

The default mode merges prediction tables, merges evaluation tables,
joins them together, and writes a single output Parquet.

```bash
python create_vsm_table.py \
  --prediction_tables PRED_DIR_OR_URIS \
  --evaluation_tables EVAL_URIS \
  --output OUTPUT_PATH \
  [options]
```

### Subtable modes

Write only the prediction or evaluation side, skipping the final
eval-pred join.

```bash
# Evaluation subtable (outer join of eval tables only)
python create_vsm_table.py \
  --subtable eval \
  --evaluation_tables EVAL_URIS \
  --output eval_merged.parquet

# Prediction subtable (merge + percentiles, no eval)
python create_vsm_table.py \
  --subtable pred \
  --prediction_tables PRED_DIR \
  --join_type inner \
  --reference_score none \
  --percentile_order post \
  --output pred_inner_post.parquet
```

### Precomputed fast path

When both a precomputed prediction and evaluation subtable already
exist, load them directly, left-join, apply filters, and write output.
Skips all merging, percentile computation, negation, and aggregation.

```bash
python create_vsm_table.py \
  --precomputed_prediction pred_inner_post.parquet \
  --precomputed_evaluation eval_merged.parquet \
  --filter_tables FILTER_URIS \
  --output fast_inner.parquet
```

The fast path also supports a runtime linker join (to avoid
materializing large linked intermediates):

```bash
python create_vsm_table.py \
  --precomputed_prediction pred_inner_post.parquet \
  --precomputed_evaluation eval_merged.parquet \
  --linker_table linker_no_multimapped.parquet \
  --filter_tables GENE_FILTER_URIS \
  --output fast_inner_gene_filtered.parquet
```

### Incremental score addition

Add a new score to an existing merged prediction table without
re-running the full pipeline:

```bash
python add_score.py \
  --base pred_inner_post.parquet \
  --new_table new_score.parquet \
  --join_type inner \
  --percentile_order post \
  --output pred_inner_post_plus_new.parquet
```

For pairwise tables, use `--pairwise_anchor` to also create the
pairwise anchor-masked columns for the new score:

```bash
python add_score.py \
  --base pred_pairwise_post.parquet \
  --new_table new_score.parquet \
  --join_type outer \
  --pairwise_anchor polyphen_score \
  --percentile_order post \
  --output pred_pairwise_plus_new.parquet
```

## Join strategies

All tables are joined on the composite key `(chrom, pos, ref, alt)`.

| Strategy | `--join_type` | Behaviour |
|----------|---------------|-----------|
| Inner | `inner` | Only rows whose key appears in every prediction table are retained. Rows with any null score are dropped. |
| Outer | `outer` | All rows from all tables are retained. Missing scores are null. |
| Pairwise | `pairwise` | The anchor column (via `--anchor`) is inner-joined with each other score column individually, then all pairs are outer-joined. Preserves the full anchor column while pairing it with each other score. |

Evaluation tables are always outer-joined with each other. The final
eval-pred join is a left join (eval left, pred right).

## Percentile rank calculation

For each prediction score column, a companion `{score}_percentile`
column is computed as:

```
rank(method="average") / count_of_non_null_values
```

Rows with null scores receive a null percentile. The raw score columns
are dropped by default (use `--keep_raw_scores` to retain them).

| `--percentile_order` | Behaviour |
|----------------------|-----------|
| `pre` | Percentiles are computed on each prediction table independently before merging. Incompatible with negation and pairwise (falls back to `post`). |
| `post` (default) | Percentiles are computed on the merged prediction frame after joining. |
| `none` | Skip percentile computation; retain raw scores. |

When `--aggregate_genes` is set, `--percentile_order` is ignored;
percentiles are computed on the gene-level aggregates (mean/max)
instead.

### Pairwise percentile naming

In pairwise mode, percentile columns follow a specific naming
convention:

- `{anchor}_anchor_percentile` -- percentile of the anchor column
- `{score}_percentile_with_anchor` -- percentile of a non-anchor score
- `{anchor}_anchor_percentile_with_{score}` -- percentile of the
  pairwise product column

### Percentile thresholds

Use `--percentile_thresholds 0.5,0.9,0.95` to write a sidecar TSV
containing the raw score values at each specified quantile, saved
alongside the output as `{output}.percentile_thresholds.tsv`.

## Score negation

By default, scores are directionally aligned against a reference column
(`--reference_score AM`). Scores with negative Pearson correlation to
the reference are negated so all scores point in the same direction. Set
`--reference_score none` to disable.

## Gene-level aggregation

With `--aggregate_genes`, variant-level scores are aggregated to one row
per gene via mean and max. Requires an `ensg` column (provided by
`--linker_table`).

| Flag | Behaviour |
|------|-----------|
| `--aggregate_genes` | Collapse to one row per gene (mean + max per score) |
| `--aggregate_genes --no-collapse_genes` | Preserve original variant rows; mean/max values are repeated for all variants in a gene |

## Spatial smoothing

With `--smooth`, percentile-ranked scores are spatially smoothed using a
Gaussian kernel over 3D protein structure distances:

```bash
python create_vsm_table.py \
  --prediction_tables PRED_DIR \
  --evaluation_tables EVAL_URIS \
  --smooth \
  --smooth_reference_dir path/to/sir-reference-data/ \
  --smooth_sigma 10.0 \
  --output smoothed.parquet
```

## Filter annotation

When `--filter_tables` is provided, boolean membership columns are
left-joined onto the output. No rows are added or removed.

Filters can also be applied standalone:

```bash
python apply_filters.py \
  --reference merged.parquet \
  --filter_tables "variant_filters/,gene_filters/" \
  --output merged_filtered.parquet
```

## CLI reference: `create_vsm_table.py`

### Input arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--prediction_tables` | Yes* | Comma-separated URIs/directories of prediction tables |
| `--evaluation_tables` | Yes* | Comma-separated URIs/directories of evaluation tables |
| `--precomputed_prediction` | No | Path to precomputed prediction parquet (fast path) |
| `--precomputed_evaluation` | No | Path to precomputed evaluation parquet |
| `--linker_table` | No | Linker parquet mapping variant keys to ensg |
| `--filter_tables` | No | Comma-separated URIs/directories of filter tables |
| `--output` | Yes | Destination path for the output Parquet |

\* Not required when using `--subtable` or `--precomputed_*` alternatives.

### Mode arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--subtable` | *(none)* | `pred` or `eval` -- write only one side |
| `--join_type` | `inner` | `inner`, `outer`, or `pairwise` |
| `--anchor` | *(none)* | Anchor score column for pairwise joins |

### Processing arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--percentile_order` | `post` | `pre`, `post`, or `none` |
| `--reference_score` | `AM` | Reference column for negation (`none` to disable) |
| `--output_table_fields` | *(all)* | Comma-separated score columns to include |
| `--keep_raw_scores` | off | Retain raw score columns alongside percentiles |
| `--aggregate_genes` | off | Aggregate scores to gene level (mean + max) |
| `--collapse_genes` | on | Collapse to one row per gene (with `--aggregate_genes`) |
| `--smooth` | off | Apply spatial smoothing |
| `--smooth_reference_dir` | *(none)* | Path to sir-reference-data directory |
| `--smooth_sigma` | `10.0` | Gaussian kernel scale in angstroms |
| `--percentile_thresholds` | *(none)* | Comma-separated quantiles for sidecar TSV |

### Diagnostic arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--row_counts` | *(none)* | Write a markdown row-count report |
| `--use_cache` | off | Reuse cached TSV-to-Parquet conversions |
| `--store_cache` | off | Persist TSV-to-Parquet conversions for reuse |

## CLI reference: `add_score.py`

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--base` | Yes | | Path to existing merged prediction parquet |
| `--new_table` | Yes | | Path to the new score table |
| `--output` | Yes | | Destination path for the output |
| `--join_type` | Yes | | `inner` or `outer` |
| `--score_columns` | No | *(all Float64)* | Comma-separated score columns from new table |
| `--percentile_order` | No | `none` | `pre`, `post`, or `none` |
| `--pairwise_anchor` | No | *(none)* | Anchor name for pairwise column creation |
| `--use_cache` | No | off | Reuse cached TSV-to-Parquet conversions |
| `--store_cache` | No | off | Persist TSV-to-Parquet conversions |

## CLI reference: `apply_filters.py`

| Argument | Required | Description |
|----------|----------|-------------|
| `--reference` | Yes | URI of the reference parquet/tsv to annotate |
| `--filter_tables` | Yes | Comma-separated URIs/directories of filter tables |
| `--output` | Yes | Destination path for the annotated output |
| `--use_cache` | No | Reuse cached TSV-to-Parquet conversions |
| `--store_cache` | No | Persist TSV-to-Parquet conversions |

## Project layout

```
merge/
├── create_vsm_table.py    # Main pipeline CLI and orchestration
├── add_score.py            # Incremental score addition CLI
├── apply_filters.py        # Standalone filter annotation CLI
├── merge.py                # Multi-table join logic (inner/outer/pairwise)
├── percentile.py           # Null-safe percentile rank calculation
├── negate.py               # Correlation-based score negation
├── smooth.py               # Spatial smoothing over 3D protein structures
├── table_io.py             # Format detection, reading, writing (local + GCS)
├── row_counts.py           # Pipeline stage row-count diagnostics
├── requirements.txt        # Python dependencies
└── README.md               # This file
```

## Performance notes

* All I/O uses the Polars **Lazy API** (`scan_parquet`, `scan_csv`).
  The full computation graph is built lazily and only materialized when
  the final `sink_parquet` writes the output.
* Operations that force materialization include `rank()` (percentiles),
  `corr()` (negation), and `count().collect()` (row counts / diagnostics).
* TSV-to-Parquet conversions are not cached by default. Use
  `--store_cache` to persist conversions in `$TMPDIR/vsm_table_cache/`
  and `--use_cache` to reuse them on subsequent runs.
* The precomputed fast path avoids redundant score merging and percentile
  computation when generating multiple output tables from the same
  prediction/evaluation inputs.
* Output Parquet files use **zstd** compression.
* Gene-level non-collapsed outer joins and pairwise variant-level tables
  with many score columns are the most memory-intensive operations and
  may OOM on machines with limited RAM.
