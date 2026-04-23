# DuckDB Score Database

A persistent DuckDB-based pipeline for ingesting, deduplicating, inspecting, and merging variant/gene scores from Parquet files. Each step is a standalone Python script that can be run individually from the command line or composed into larger workflows.

All scripts accept a `--db <path>` argument to specify the database file. Override the default with this flag.

## Prerequisites

- Python 3.10+
- `duckdb` Python package (`pip install duckdb`)
- Input data as Parquet files with either variant key columns (`chrom`, `pos`, `ref`, `alt`) or a gene key column (`ensg`)

## Modules

### `initialize_db.py` — Create the database

Creates a new `.duckdb` file containing a single `metadata` table that tracks all registered score, eval, and merged tables. Refuses to overwrite an existing database.

The metadata table schema:

| Column          | Type    | Description                                              |
|-----------------|---------|----------------------------------------------------------|
| `source_column` | VARCHAR | Source column name from the Parquet file — the score column name for score tables, the eval column name for eval tables, or a comma-separated list of source names for merged tables |
| `source_path`   | VARCHAR | Absolute path to the source Parquet file for score/eval tables, or the set operation name (e.g. `intersection`) for merged tables |
| `table_name`    | VARCHAR | Internal DuckDB table name (primary key)                 |
| `table_type`    | VARCHAR | One of `'score'`, `'eval'`, `'merged_scores'`, `'merged_evals'`, `'merged_analysis'` |
| `analysis_level`    | VARCHAR | Key type: `'variant'` (chrom/pos/ref/alt) or `'gene'` (ensg) |
| `deduped`       | BOOLEAN | Whether duplicate keys have been resolved                |

```bash
python duckdb/initialize_db.py --db scores.duckdb
```

---

### `ingest_score.py` — Load a score column from Parquet

Reads a single numeric score column from a Parquet file and materializes it as a DuckDB table. The table is physically sorted by score and includes a precomputed dense rank.

Each row has the following columns:

| Column   | Type     | Description |
|----------|----------|-------------|
| `chrom`  | VARCHAR  | Chromosome (NULL for gene-keyed tables) |
| `pos`    | BIGINT   | Position (NULL for gene-keyed tables) |
| `ref`    | VARCHAR  | Reference allele (NULL for gene-keyed tables) |
| `alt`    | VARCHAR  | Alternate allele (NULL for gene-keyed tables) |
| `ensg`   | VARCHAR  | Ensembl gene ID (NULL for variant-keyed tables) |
| `key`    | UBIGINT  | Hash of key columns: `hash(chrom\|pos\|ref\|alt)` for variant, `hash(ensg)` for gene |
| `score`  | DOUBLE   | The ingested score value |
| `temp_1` | DOUBLE   | Reserved for percentile computation (initially NULL) |
| `temp_2` | DOUBLE   | Reserved for percentile computation (initially NULL) |
| `rank`   | INTEGER  | Dense rank of non-null scores (NULL for null scores) |

```bash
python duckdb/ingest_score.py \
  --db scores.duckdb \
  --score_name revel \
  --score_path data/gnomad_chr22.parquet \
  --table_name revel_chr22 \
  --analysis_level variant
```

After ingestion the script checks for duplicate keys and sets `metadata.deduped` accordingly. A warning is printed if duplicates are detected.

---

### `ingest_eval.py` — Load a boolean label column from Parquet

Similar to `ingest_score.py`, but stores a boolean `is_pos` column instead of score/temp/rank. Used for evaluation labels (e.g., pathogenicity truth sets).

```bash
python duckdb/ingest_eval.py \
  --db scores.duckdb \
  --eval_name is_pathogenic \
  --eval_path data/clinvar_labels.parquet \
  --table_name clinvar_eval \
  --analysis_level variant
```

---

### `eject_table.py` — Remove a table from the database

Drops the named table and deletes its row from the metadata table. Works for score, eval, and merged tables alike.

```bash
python duckdb/eject_table.py --db scores.duckdb --table_name revel_chr22
```

---

### `inspect_db.py` — View database contents and statistics

Lists all registered tables grouped by type and reports summary statistics:

- **Score tables**: row count, scored/null counts, unique keys, min/max/mean/median of score, dedup status, duplicate warnings
- **Eval tables**: row count, positive/negative/null counts, unique keys, dedup status
- **Merged tables**: row count, column count, data column list, source score names

Optionally print randomly sampled rows with `--sample N`.

```bash
python duckdb/inspect_db.py --db scores.duckdb
python duckdb/inspect_db.py --db scores.duckdb --sample 5
```

---

### `remove_duplicates.py` — Deduplicate a table by key

Removes rows with duplicate `key` values. Required before pairwise operations or merging.

Two strategies are available:

| Strategy        | Behavior |
|-----------------|----------|
| `keep_random`   | Keeps one random row per duplicate key (works on both score and eval tables) |
| `prefer_scored` | Keeps the row with a non-null score, breaking ties at random (score tables only) |

```bash
python duckdb/remove_duplicates.py \
  --db scores.duckdb \
  --table_name revel_chr22 \
  --strategy prefer_scored
```

---

### `pairwise_percentile.py` — In-place pairwise percentile computation

Computes intersection-based `CUME_DIST` between an anchor and a non-anchor score table, writing the results into the non-anchor table's `temp_1` and `temp_2` columns:

- `temp_1` = anchor score percentile within the key intersection
- `temp_2` = non-anchor score percentile within the key intersection
- Rows outside the intersection retain NULL
- NULL scores within the intersection are excluded from the ranking denominator

Both tables must be of type `'score'` and must be deduped.

```bash
python duckdb/pairwise_percentile.py \
  --db scores.duckdb \
  --anchor_table am_chr22 \
  --table_name revel_chr22
```

---

### `merge_scores.py` — Combine multiple score tables

Merges two or more deduped score tables into a single wide output table. The merge is controlled by two orthogonal flags:

**Set operations** (`--set_operation`):

| Operation      | Join type       | Description |
|----------------|-----------------|-------------|
| `intersection` | INNER JOIN      | Only rows present in every input table |
| `union`        | FULL OUTER JOIN | All rows from any input table; missing scores are NULL |
| `pairwise`     | Per-pair INNER, then FULL OUTER | Anchor-based pairwise percentiles across pair intersections |

**Percentile modes** (`--percentile`):

| Mode   | Description |
|--------|-------------|
| `pre`  | Percentile each table's score (full population) *before* joining |
| `post` | Join first, then percentile each score column across the joined result |
| `none` | Raw scores only, no percentile columns |

Constraints:
- All input tables must be deduped score tables
- `--anchor_table` is required for pairwise and forbidden otherwise
- `--percentile` must be `none` for pairwise (it computes its own intersection-based percentiles)
- `--output_table` must not already exist
- At least two input tables are required

The output table is registered in metadata as `table_type='merged_scores'`.

```bash
# Intersection with post-percentile
python duckdb/merge_scores.py \
  --db scores.duckdb \
  --tables revel_chr22 am_chr22 cadd_chr22 \
  --output_table merged_intersection \
  --set_operation intersection \
  --percentile post

# Union with pre-percentile
python duckdb/merge_scores.py \
  --db scores.duckdb \
  --tables revel_chr22 am_chr22 \
  --output_table merged_union \
  --set_operation union \
  --percentile pre

# Pairwise with AM as anchor
python duckdb/merge_scores.py \
  --db scores.duckdb \
  --tables am_chr22 revel_chr22 cadd_chr22 \
  --output_table merged_pairwise \
  --set_operation pairwise \
  --percentile none \
  --anchor_table am_chr22
```

**Output column structure by mode:**

- **intersection/union, `none`**: key columns + one raw score column per input (`{score_name}`)
- **intersection/union, `pre` or `post`**: key columns + raw score + percentile per input (`{score_name}`, `{score_name}_percentile`)
- **pairwise**: key columns + per non-anchor pair: anchor raw score, non-anchor raw score, anchor pairwise percentile (`{anchor}_pairwise_{C}`), non-anchor pairwise percentile (`{C}_pairwise_{anchor}`)

---

### `merge_evals.py` — Combine multiple eval tables

Merges two or more deduped eval tables into a single wide table via `FULL OUTER JOIN` on `key`. Each input eval table's `is_pos` column is renamed to the table's `source_column` name so they are distinct in the output.

The output table is registered in metadata as `table_type='merged_evals'`.

Constraints:
- All input tables must be deduped eval tables
- `--output_table` must not already exist
- At least two input tables are required

```bash
python duckdb/merge_evals.py \
  --db scores.duckdb \
  --tables clinvar_eval omim_eval \
  --output_table merged_labels
```

**Output column structure:** key columns (`chrom`, `pos`, `ref`, `alt`, `ensg`, `key`) + one boolean column per input eval, named after the eval's `source_column`.

---

### `create_analysis_table.py` — Join merged scores and merged evals

Combines a `merged_scores` table and a `merged_evals` table into a single `merged_analysis` table via a `FULL OUTER JOIN`. This is the final assembly step that brings score data and evaluation labels together into one wide table for downstream analysis.

**Same-key join** (default): When both input tables share the same `analysis_level` (both `variant` or both `gene`), the join is a straightforward `FULL OUTER JOIN` on `"key"`, coalescing key columns. The output `analysis_level` matches the shared level.

**Cross-key join** (`--linker_path`): When the two input tables have different `analysis_level` values (one `variant`, one `gene`), a linker Parquet file is required. The linker must contain `chrom`, `pos`, `ref`, `alt`, and `ensg` columns, establishing the many-to-one mapping from variant keys to gene keys. The output `analysis_level` is `variant`.

The output table is registered in metadata as `table_type='merged_analysis'`.

Constraints:
- `--scores_table` must be `table_type='merged_scores'` in metadata
- `--evals_table` must be `table_type='merged_evals'` in metadata
- `--output_table` must not already exist
- `--linker_path` is required when `analysis_level` differs between the two inputs, and forbidden when it matches

```bash
# Same-key join (both tables are variant-level)
python duckdb/create_analysis_table.py \
  --db scores.duckdb \
  --scores_table merged_scores_tbl \
  --evals_table merged_evals_tbl \
  --output_table my_analysis

# Cross-key join (variant scores + gene evals via linker)
python duckdb/create_analysis_table.py \
  --db scores.duckdb \
  --scores_table merged_variant_scores \
  --evals_table merged_gene_evals \
  --output_table cross_analysis \
  --linker_path data/variant_gene_linker.parquet
```

**Output column structure:** key columns (`chrom`, `pos`, `ref`, `alt`, `ensg`, `key`) + all data columns from the scores table (raw scores, percentiles, pairwise columns) + all data columns from the evals table (boolean label columns).

---

### `join_linker.py` — Enrich a merged table with key columns from a linker

Fills in the "other side" key columns on any merged table (`merged_scores`, `merged_evals`, or `merged_analysis`) using a linker Parquet that maps between variant keys and gene keys.

**Variant-level table + linker (add `ensg`):** LEFT JOIN on variant `key` to pull in `ensg` from the linker. This is many-to-one, so the row count does not change — it simply fills in the `ensg` column where a match exists. The `analysis_level` remains `variant`.

**Gene-level table + linker (add `chrom/pos/ref/alt`):** LEFT JOIN on `ensg` to pull in variant key columns from the linker. This is one-to-many, so the row count **expands** (each gene row becomes N variant rows). The `analysis_level` changes to `variant`.

The operation is **in-place** — it replaces the existing table rather than creating a new one.

Constraints:
- The table must be a merged type (`merged_scores`, `merged_evals`, or `merged_analysis`)
- The linker Parquet must contain `chrom`, `pos`, `ref`, `alt`, and `ensg` columns

```bash
# Add ensg to a variant-level merged table
python duckdb/join_linker.py \
  --db scores.duckdb \
  --table_name merged_scores_tbl \
  --linker_path data/linker.parquet

# Expand a gene-level merged table to variant granularity
python duckdb/join_linker.py \
  --db scores.duckdb \
  --table_name merged_gene_evals \
  --linker_path data/linker.parquet
```

---

## Typical Workflows

### Basic: Ingest, inspect, and merge scores

```bash
# 1. Create the database
python duckdb/initialize_db.py --db my_scores.duckdb

# 2. Ingest scores from Parquet files
python duckdb/ingest_score.py --db my_scores.duckdb \
  --score_name revel --score_path data/revel_chr22.parquet \
  --table_name revel_chr22 --analysis_level variant

python duckdb/ingest_score.py --db my_scores.duckdb \
  --score_name AM --score_path data/alphamissense_chr22.parquet \
  --table_name am_chr22 --analysis_level variant

python duckdb/ingest_score.py --db my_scores.duckdb \
  --score_name cadd_score --score_path data/cadd_chr22.parquet \
  --table_name cadd_chr22 --analysis_level variant

# 3. Deduplicate all score tables
python duckdb/remove_duplicates.py --db my_scores.duckdb \
  --table_name revel_chr22 --strategy prefer_scored
python duckdb/remove_duplicates.py --db my_scores.duckdb \
  --table_name am_chr22 --strategy prefer_scored
python duckdb/remove_duplicates.py --db my_scores.duckdb \
  --table_name cadd_chr22 --strategy prefer_scored

# 4. Inspect the database
python duckdb/inspect_db.py --db my_scores.duckdb --sample 3

# 5. Merge into a single wide table with post-join percentiles
python duckdb/merge_scores.py --db my_scores.duckdb \
  --tables revel_chr22 am_chr22 cadd_chr22 \
  --output_table merged_all \
  --set_operation intersection --percentile post

# 6. Inspect the result
python duckdb/inspect_db.py --db my_scores.duckdb
```

### Replacing a score and re-merging

If you realize a score table was ingested from the wrong source or needs to be updated:

```bash
# 1. Remove the old score and the existing merged table
python duckdb/eject_table.py --db my_scores.duckdb --table_name revel_chr22
python duckdb/eject_table.py --db my_scores.duckdb --table_name merged_all

# 2. Re-ingest from the corrected Parquet
python duckdb/ingest_score.py --db my_scores.duckdb \
  --score_name revel --score_path data/revel_chr22_v2.parquet \
  --table_name revel_chr22 --analysis_level variant

# 3. Deduplicate
python duckdb/remove_duplicates.py --db my_scores.duckdb \
  --table_name revel_chr22 --strategy prefer_scored

# 4. Re-merge
python duckdb/merge_scores.py --db my_scores.duckdb \
  --tables revel_chr22 am_chr22 cadd_chr22 \
  --output_table merged_all \
  --set_operation intersection --percentile post
```

### Pairwise merge with a specific anchor

When you want to compute intersection-based percentiles for each score relative to one anchor:

```bash
# Prerequisite: all tables ingested and deduped (see above)

# Merge with AM as the anchor score
python duckdb/merge_scores.py --db my_scores.duckdb \
  --tables am_chr22 revel_chr22 cadd_chr22 \
  --output_table pairwise_am_anchor \
  --set_operation pairwise --percentile none \
  --anchor_table am_chr22

# Inspect to see the pairwise columns
python duckdb/inspect_db.py --db my_scores.duckdb --sample 3
```

This produces columns like `AM_pairwise_revel`, `revel_pairwise_AM`, `AM_pairwise_cadd_score`, `cadd_score_pairwise_AM` — each representing the `CUME_DIST` percentile of the respective score within the key intersection of that anchor-nonanchor pair.

### Standalone pairwise percentile (in-place)

If you need to write pairwise percentiles directly into a score table's temp columns (e.g., for downstream analysis outside of the merge pipeline):

```bash
# Write anchor and non-anchor percentiles into revel_chr22's temp_1/temp_2
python duckdb/pairwise_percentile.py --db my_scores.duckdb \
  --anchor_table am_chr22 --table_name revel_chr22

# Inspect the temp columns
python duckdb/inspect_db.py --db my_scores.duckdb --sample 5
```

### Full analysis: scores + evals → analysis table

A complete workflow from raw Parquet files to a unified analysis table:

```bash
# 1. Create the database
python duckdb/initialize_db.py --db my_analysis.duckdb

# 2. Ingest scores
python duckdb/ingest_score.py --db my_analysis.duckdb \
  --score_name revel --score_path data/revel_chr22.parquet \
  --table_name revel_chr22 --analysis_level variant
python duckdb/ingest_score.py --db my_analysis.duckdb \
  --score_name AM --score_path data/am_chr22.parquet \
  --table_name am_chr22 --analysis_level variant

# 3. Ingest evals
python duckdb/ingest_eval.py --db my_analysis.duckdb \
  --eval_name is_pathogenic --eval_path data/clinvar.parquet \
  --table_name clinvar_eval --analysis_level variant
python duckdb/ingest_eval.py --db my_analysis.duckdb \
  --eval_name is_lof --eval_path data/lof_labels.parquet \
  --table_name lof_eval --analysis_level variant

# 4. Deduplicate everything
python duckdb/remove_duplicates.py --db my_analysis.duckdb \
  --table_name revel_chr22 --strategy prefer_scored
python duckdb/remove_duplicates.py --db my_analysis.duckdb \
  --table_name am_chr22 --strategy prefer_scored
python duckdb/remove_duplicates.py --db my_analysis.duckdb \
  --table_name clinvar_eval --strategy keep_random
python duckdb/remove_duplicates.py --db my_analysis.duckdb \
  --table_name lof_eval --strategy keep_random

# 5. Merge scores and evals separately
python duckdb/merge_scores.py --db my_analysis.duckdb \
  --tables revel_chr22 am_chr22 \
  --output_table merged_scores \
  --set_operation intersection --percentile post
python duckdb/merge_evals.py --db my_analysis.duckdb \
  --tables clinvar_eval lof_eval \
  --output_table merged_evals

# 6. Create the final analysis table
python duckdb/create_analysis_table.py --db my_analysis.duckdb \
  --scores_table merged_scores \
  --evals_table merged_evals \
  --output_table final_analysis

# 7. Inspect
python duckdb/inspect_db.py --db my_analysis.duckdb --sample 3
```

---

### Adding evaluation labels alongside scores

```bash
# Ingest a boolean label column
python duckdb/ingest_eval.py --db my_scores.duckdb \
  --eval_name is_pathogenic --eval_path data/clinvar_labels.parquet \
  --table_name clinvar_eval --analysis_level variant

# Deduplicate if needed
python duckdb/remove_duplicates.py --db my_scores.duckdb \
  --table_name clinvar_eval --strategy keep_random

# View all tables including eval
python duckdb/inspect_db.py --db my_scores.duckdb
```

## Design Notes

- **Key hashing**: Variant keys are hashed via `hash(chrom || '|' || pos || '|' || ref || '|' || alt)` using DuckDB's built-in `hash()` function. Gene-keyed tables are hashed via `hash(ensg)`. Both produce a `UBIGINT` key used for joins and deduplication.
- **Physical sort order**: Score tables are stored sorted by score at ingestion time. This benefits subsequent window-function operations (DENSE_RANK, CUME_DIST) by aligning with DuckDB's optimizer expectations.
- **NULL-safe percentiles**: All percentile computations use `CUME_DIST() OVER (PARTITION BY (score IS NOT NULL) ORDER BY score)` wrapped in `CASE WHEN ... IS NOT NULL` to exclude NULL scores from the ranking denominator and assign NULL percentiles to NULL scores. `CUME_DIST` is used rather than `PERCENT_RANK` to guarantee the maximum value in any percentile column is always 1.0.
- **Deduplication as explicit step**: Duplicate keys are detected at ingestion but not automatically removed. The `deduped` metadata flag gates access to pairwise operations and merging, ensuring data integrity without silent data loss.
- **Merged table metadata**: Merged tables store the comma-separated source column names in `source_column` and the set operation in `source_path` for traceability.
