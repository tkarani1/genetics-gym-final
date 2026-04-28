# DuckDB Score Database

A persistent DuckDB-based pipeline for ingesting, deduplicating, inspecting, merging, and exporting variant- and gene-level scores and evaluations from Parquet files. Each step is a standalone Python script that can be run individually from the command line or composed into larger workflows via shell scripts.

All scripts accept a `--db <path>` argument to specify the database file (default: `gg_data.duckdb`).

## Prerequisites

- Python 3.10+
- `duckdb` Python package (`pip install duckdb`)
- Input data as Parquet (or TSV/CSV for evals) with either variant key columns (`chrom`, `pos`, `ref`, `alt`) or a gene key column (`ensg`)

## Architecture Overview

The pipeline uses a **wide table model** for scores: all score columns for a given analysis level live in a single table (e.g. `variant_scores` or `gene_scores`) rather than one table per score. This dramatically reduces storage (key columns are stored once) and simplifies merges to single-table scans.

Eval tables remain individual (one table per evaluation source) and are merged into wide tables on demand via cascading outer joins.

All outputs from score merges and pairwise operations are exported directly to **Parquet files** rather than stored in the database, keeping the database compact and avoiding OOM issues during large joins.

Every mutating operation is recorded in an **audit log** for full provenance tracking.

## Modules

### `initialize_db.py` — Create the database

Creates a new `.duckdb` file containing two system tables:

**`metadata`** — tracks every ingested source:

| Column | Type | Description |
|---|---|---|
| `source_column` | VARCHAR | Score column name (scores), eval label (evals), or comma-separated list (merged) |
| `source_path` | VARCHAR | Absolute path to source file, or operation name for merged tables |
| `table_name` | VARCHAR | Internal DuckDB table name |
| `table_type` | VARCHAR | `'score'`, `'eval'`, `'merged_scores'`, `'merged_evals'`, or `'merged_analysis'` |
| `analysis_level` | VARCHAR | `'variant'` (chrom/pos/ref/alt) or `'gene'` (ensg) |
| `deduped` | BOOLEAN | Whether duplicate keys have been resolved |
| `eval_column` | VARCHAR | Source column read as `is_pos` (eval tables only, nullable) |
| `case_column` | VARCHAR | Source column read as `n_case` (eval tables only, nullable) |
| `ctrl_column` | VARCHAR | Source column read as `n_ctrl` (eval tables only, nullable) |

Unique constraint: `(table_name, source_column)` — supports wide tables with multiple metadata rows per table.

**`audit_log`** — append-only event journal:

| Column | Type | Description |
|---|---|---|
| `ts` | TIMESTAMP | Event timestamp (auto-populated) |
| `module` | VARCHAR | Python module that produced the event |
| `action` | VARCHAR | Action type (e.g. `create_table`, `add_column`, `export_parquet`) |
| `table_name` | VARCHAR | Affected table (nullable) |
| `details` | VARCHAR | Free-text details (nullable) |

```bash
python duckdb/initialize_db.py --db my_data.duckdb
```

---

### `ingest_score.py` — Add a score column to a wide table

Ingests a single numeric score column from a Parquet file into a **wide score table**. Defaults to `variant_scores` or `gene_scores` based on `--analysis_level`, overridable with `--table_name`.

- If the wide table doesn't exist yet, creates it with canonical key columns plus the score column.
- If it exists, adds the column via `ALTER TABLE ADD COLUMN`, then updates existing keys and inserts new keys (with NULLs for all other score columns).
- Deduplicates on ingest: keeps the row with a non-null score when duplicate keys exist.

Wide table schema:

| Column | Type | Description |
|---|---|---|
| `chrom` | VARCHAR | Chromosome (NULL for gene-keyed tables) |
| `pos` | BIGINT | Position (NULL for gene-keyed tables) |
| `ref` | VARCHAR | Reference allele (NULL for gene-keyed tables) |
| `alt` | VARCHAR | Alternate allele (NULL for gene-keyed tables) |
| `ensg` | VARCHAR | Ensembl gene ID (NULL for variant-keyed tables) |
| `key` | UBIGINT | Hash of key columns |
| *score_1* | DOUBLE | First ingested score |
| *score_2* | DOUBLE | Second ingested score |
| ... | ... | One column per ingested score |

```bash
python duckdb/ingest_score.py \
  --db my_data.duckdb \
  --score_name AM_score \
  --score_path data/alphamissense.parquet \
  --analysis_level variant

python duckdb/ingest_score.py \
  --db my_data.duckdb \
  --score_name loeuf_v2_score \
  --score_path data/loeuf.parquet \
  --table_name gene_scores \
  --analysis_level gene
```

---

### `ingest_eval.py` — Load evaluation data

Ingests an evaluation table from a Parquet, TSV, or CSV file. Each eval table has a fixed schema of three data columns — which are populated depends on the arguments:

| Column | Type | Populated when |
|---|---|---|
| `is_pos` | BOOLEAN | `--eval_column` is given (source column cast to boolean) |
| `n_case` | INTEGER | `--case_column` is given (must pair with `--ctrl_column`) |
| `n_ctrl` | INTEGER | `--ctrl_column` is given (must pair with `--case_column`) |

This allows flexible ingestion regardless of the source column names — e.g. `--eval_column is_case` reads a column named `is_case` as `is_pos`, and `--case_column n_case_subset` reads `n_case_subset` as `n_case`. The original source column names are recorded in metadata for provenance.

```bash
# Boolean eval (variant-level)
python duckdb/ingest_eval.py \
  --db my_data.duckdb \
  --eval_name clinvar \
  --eval_path data/clinvar.parquet \
  --table_name clinvar_eval \
  --analysis_level variant \
  --eval_column is_pos

# Count eval (gene-level)
python duckdb/ingest_eval.py \
  --db my_data.duckdb \
  --eval_name asc \
  --eval_path data/asc_ensg_eval.tsv \
  --table_name asc_gene_eval \
  --analysis_level gene \
  --case_column n_case --ctrl_column n_ctrl

# Both boolean and count columns
python duckdb/ingest_eval.py \
  --db my_data.duckdb \
  --eval_name genebass \
  --eval_path data/genebass.tsv \
  --table_name genebass_gene_eval \
  --analysis_level gene \
  --eval_column is_pos \
  --case_column n_case --ctrl_column n_ctrl
```

---

### `remove_duplicates.py` — Deduplicate a table by key

Removes rows with duplicate `key` values. Required for eval tables before merging (score tables are deduplicated at ingest).

| Strategy | Behavior |
|---|---|
| `keep_random` | Keeps one random row per duplicate key (works on both score and eval tables) |
| `prefer_scored` | Keeps the row with a non-null score, breaking ties at random (score tables only) |

```bash
python duckdb/remove_duplicates.py \
  --db my_data.duckdb \
  --table_name clinvar_eval \
  --strategy keep_random
```

---

### `eject_table.py` — Remove a table or score column

Two modes of operation:

- **Drop entire table** (default): removes the table and all its metadata rows.
- **Drop a score column** (`--score_column`): drops only that column from a wide score table and its metadata entry. If it's the last score column, the entire table is removed.

```bash
# Drop an entire eval table
python duckdb/eject_table.py --db my_data.duckdb --table_name clinvar_eval

# Drop a single score column from the wide table
python duckdb/eject_table.py --db my_data.duckdb \
  --table_name variant_scores --score_column AM_score
```

---

### `inspect_db.py` — View database contents and statistics

Lists all registered tables grouped by type:

- **Score tables**: per-column scored/null counts, min/max/mean/median, source path.
- **Eval tables**: positive/negative/null counts (boolean), sum case/ctrl (count), source path, column mapping.
- **Merged tables**: row count, column listing, source provenance.
- **Audit log** (`--audit`): event history with timestamps, modules, actions, and details.

```bash
python duckdb/inspect_db.py --db my_data.duckdb
python duckdb/inspect_db.py --db my_data.duckdb --sample 5
python duckdb/inspect_db.py --db my_data.duckdb --audit
python duckdb/inspect_db.py --db my_data.duckdb --audit-last 20
```

---

### `export_table.py` — Export a table to Parquet

Writes any database table to a Parquet file.

```bash
python duckdb/export_table.py --db my_data.duckdb \
  --table_name variant_scores \
  --output_path output/variant_scores.parquet
```

---

### `merge_evals.py` — Combine eval tables into a wide table

Merges N eval tables into a single wide database table via **cascading 2-way FULL OUTER JOINs** (bounded memory usage). Boolean columns are renamed to the source label; count columns are prefixed with the source label.

Output column structure: key columns + per eval source:
- `{label}` (BOOLEAN) — if the source had a boolean column
- `{label}_n_case` (INTEGER) — if the source had count columns
- `{label}_n_ctrl` (INTEGER) — if the source had count columns

```bash
python duckdb/merge_evals.py --db my_data.duckdb \
  --tables clinvar_eval asd_eval dd_eval \
  --output_table variant_evals_merged \
  --memory_limit 8GB
```

---

### `merge_scores.py` — Export score merges and pairwise analysis

The primary output engine. Reads from a wide score table and writes to Parquet files.

**Set operations** (`--set_operation`):

| Operation | Description |
|---|---|
| `intersection` | Only rows where all selected score columns are non-null → single Parquet |
| `union` | All rows (NULLs preserved) → single Parquet |
| `pairwise` | Anchor-based pairwise percentiles → one Parquet file per anchor-target pair |

**Percentile modes** (`--percentile`, for intersection/union only):

| Mode | Description |
|---|---|
| `pre`/`post` | Add `CUME_DIST` percentile columns for each score |
| `none` | Raw scores only |

**Pairwise mode** computes its own percentiles and supports additional features:

| Flag | Description |
|---|---|
| `--evals_table` | Left-join a merged-evals table, embedding all eval columns in each output file |
| `--gene_average` | Join a linker Parquet and add `AVG(percentile) OVER (PARTITION BY ensg)` columns |
| `--linker_path` | Path to linker Parquet (required with `--gene_average`) |
| `--memory_limit` | DuckDB memory limit (e.g. `8GB`) |

Each pairwise output file contains:
- Key columns: `chrom`, `pos`, `ref`, `alt`, `ensg`, `key`
- 2 raw scores (anchor + target)
- 2 pairwise `CUME_DIST` percentile columns
- (with `--gene_average`) 2 gene-averaged percentile columns (`_gene_avg` suffix)
- (with `--evals_table`) all eval columns

```bash
# Intersection with percentiles
python duckdb/merge_scores.py --db my_data.duckdb \
  --wide_table variant_scores \
  --columns AM_score cadd_score polyphen_score \
  --output_path output/intersection.parquet \
  --set_operation intersection --percentile post

# Pairwise with eval join
python duckdb/merge_scores.py --db my_data.duckdb \
  --wide_table variant_scores \
  --columns polyphen_score AM_score cadd_score \
  --output_dir output/pairwise/ \
  --set_operation pairwise --percentile none \
  --anchor_column polyphen_score \
  --evals_table variant_evals_merged \
  --memory_limit 8GB

# Pairwise with gene-averaged percentiles
python duckdb/merge_scores.py --db my_data.duckdb \
  --wide_table variant_scores \
  --columns polyphen_score AM_score cadd_score \
  --output_dir output/pairwise_gene_avg/ \
  --set_operation pairwise --percentile none \
  --anchor_column polyphen_score \
  --evals_table variant_evals_merged \
  --gene_average \
  --linker_path data/linker_all.parquet \
  --memory_limit 8GB
```

---

### `pairwise_percentile.py` — Single-pair percentile export

Computes pairwise `CUME_DIST` percentiles between two specific score columns in a wide table and writes the result to a single Parquet file. Useful for ad-hoc single-pair analysis; for batch operations use pairwise mode in `merge_scores.py`.

```bash
python duckdb/pairwise_percentile.py --db my_data.duckdb \
  --wide_table variant_scores \
  --anchor_column polyphen_score \
  --target_column AM_score \
  --output_path output/polyphen_x_am.parquet
```

---

### `create_analysis_table.py` — Join scores and evals in-database

Combines a scores table and a merged-evals table into a `merged_analysis` database table via FULL OUTER JOIN.

- **Same-level join**: Both tables share `analysis_level` → direct key join.
- **Cross-level join** (`--linker_path`): Different levels (e.g. variant scores + gene evals) → three-way join via linker.

Accepts both `score` and `merged_scores` types for the scores input.

```bash
# Same-level join
python duckdb/create_analysis_table.py --db my_data.duckdb \
  --scores_table variant_scores \
  --evals_table variant_evals_merged \
  --output_table variant_analysis

# Cross-level join
python duckdb/create_analysis_table.py --db my_data.duckdb \
  --scores_table variant_scores \
  --evals_table gene_evals_merged \
  --output_table cross_analysis \
  --linker_path data/linker_all.parquet
```

---

### `join_linker.py` — Enrich a merged table with key columns

Fills in the "other side" key columns on any merged table using a linker Parquet. The operation is **in-place**.

- **Variant table**: fills `ensg` (many-to-one, row count unchanged).
- **Gene table**: expands to variant granularity (one-to-many, row count increases, `analysis_level` changes to `variant`).

```bash
python duckdb/join_linker.py --db my_data.duckdb \
  --table_name variant_evals_merged \
  --linker_path data/linker_all.parquet
```

---

## Typical Workflows

### Wide table ingestion

```bash
# 1. Create the database
python duckdb/initialize_db.py --db my_data.duckdb

# 2. Ingest multiple scores into a single wide table
python duckdb/ingest_score.py --db my_data.duckdb \
  --score_name polyphen_score --score_path data/all.parquet --analysis_level variant
python duckdb/ingest_score.py --db my_data.duckdb \
  --score_name AM_score --score_path data/julia.parquet --analysis_level variant
python duckdb/ingest_score.py --db my_data.duckdb \
  --score_name cadd_score --score_path data/all.parquet --analysis_level variant

# 3. Inspect
python duckdb/inspect_db.py --db my_data.duckdb
```

### Full pairwise analysis with evals and gene averaging

```bash
# 1. Ingest scores (as above)

# 2. Ingest evals
python duckdb/ingest_eval.py --db my_data.duckdb \
  --eval_name clinvar --eval_path data/clinvar.parquet \
  --table_name clinvar_eval --analysis_level variant --eval_column is_pos
python duckdb/ingest_eval.py --db my_data.duckdb \
  --eval_name asd --eval_path data/asd.parquet \
  --table_name asd_eval --analysis_level variant --eval_column is_pos

# 3. Deduplicate evals
python duckdb/remove_duplicates.py --db my_data.duckdb \
  --table_name clinvar_eval --strategy keep_random
python duckdb/remove_duplicates.py --db my_data.duckdb \
  --table_name asd_eval --strategy keep_random

# 4. Merge evals into a wide table
python duckdb/merge_evals.py --db my_data.duckdb \
  --tables clinvar_eval asd_eval \
  --output_table merged_evals

# 5. Export pairwise with evals and gene averaging
python duckdb/merge_scores.py --db my_data.duckdb \
  --wide_table variant_scores \
  --columns polyphen_score AM_score cadd_score \
  --output_dir output/pairwise/ \
  --set_operation pairwise --percentile none \
  --anchor_column polyphen_score \
  --evals_table merged_evals \
  --gene_average \
  --linker_path data/linker_all.parquet \
  --memory_limit 8GB
```

### Replacing a score column

```bash
# 1. Remove the old column
python duckdb/eject_table.py --db my_data.duckdb \
  --table_name variant_scores --score_column AM_score

# 2. Re-ingest from corrected source
python duckdb/ingest_score.py --db my_data.duckdb \
  --score_name AM_score --score_path data/am_v2.parquet --analysis_level variant
```

---

## Design Notes

- **Wide table model**: All score columns for a given analysis level coexist in a single table. Key columns are stored once, reducing database size by ~10x compared to one-table-per-score. New scores are added via `ALTER TABLE ADD COLUMN` + `UPDATE`/`INSERT`.
- **Key hashing**: Variant keys are hashed via `hash(chrom || '|' || pos || '|' || ref || '|' || alt)` using DuckDB's built-in `hash()` function. Gene keys use `hash(ensg)`. Both produce a `UBIGINT` used for joins and deduplication.
- **NULL-safe percentiles**: All percentile computations use `CUME_DIST() OVER (PARTITION BY (score IS NOT NULL) ORDER BY score)` wrapped in `CASE WHEN ... IS NOT NULL` to exclude NULL scores from the ranking denominator. `CUME_DIST` is used rather than `PERCENT_RANK` to guarantee the maximum percentile is always 1.0.
- **Per-pair Parquet output**: Pairwise mode writes one file per anchor-target pair. This bounds memory usage to a single pair's worth of computation and provides natural resumability (existing files are skipped).
- **Gene averaging**: When `--gene_average` is set, pairwise percentiles are computed first, then a linker is joined to obtain `ensg`, and `AVG(percentile) OVER (PARTITION BY ensg)` is computed. Both per-variant and gene-averaged percentiles are included in the output.
- **Cascading joins**: Eval merges use cascading 2-way FULL OUTER JOINs rather than a single N-way join, bounding memory usage to one intermediate table at a time.
- **Explicit eval column tracking**: Metadata records the original source column names (`eval_column`, `case_column`, `ctrl_column`) for each eval table, allowing ingestion of non-standard column names (e.g. `is_case`, `is_pos_fine_mapped`, `n_case_subset`) while normalizing to `is_pos`/`n_case`/`n_ctrl` internally.
- **Audit log**: Every mutating operation appends to the `audit_log` table, recording the module, action, affected table, and details. This enables lineage tracking and dependency analysis.
- **Deduplication**: Score tables are deduplicated at ingest (prefer non-null score). Eval tables require explicit deduplication via `remove_duplicates.py` before merging.
