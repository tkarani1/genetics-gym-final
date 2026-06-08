# `src/`

Scripts for building the `write_tables_pipeline` dataset. All scripts are intended to be
run from this `src/` directory and read/write to the parallel `../data/` tree.
Paths are resolved relative to each script's own location, so they also work
when invoked from elsewhere.

The easiest way to build everything is the orchestrator
[`full_data_aggregation_pipeline.py`](#full_data_aggregation_pipelinepy), which
runs every script below in dependency order from a single
[`config.json`](#configjson) and writes an audit trail of the run. The
remaining sections document each underlying script, which can also be run on its
own.

```
write_tables_pipeline/
├── data/
│   ├── raw_data/
│   │   ├── input_data_locations.json       # GCS source manifest (downloads)
│   │   ├── scores/                          # downloaded score_tables
│   │   ├── evals/                           # downloaded eval_tables
│   │   ├── linker/                          # downloaded variant->gene linker
│   │   └── filters/                         # downloaded filter_tables
│   │       ├── variant/                     #   variant-level filters (.tsv.gz)
│   │       └── ensg/                        #   gene/ensg-level filter (.tsv)
│   └── processed_data/
│       ├── scores/
│       │   ├── score_input_data.json        # which score fields to merge
│       │   ├── variant_scores_all_outer.parquet  # wide score table (outer)
│       │   ├── variant_scores_all_inner.parquet  # drop-NA (inner) score table
│       │   ├── variant_scores_*_percentile.parquet  # percentile-transformed tables
│       │   ├── variant_scores_percentile_thresholds.tsv
│       │   ├── gene_aggregated/             # gene (ensg)-annotated + per-gene stats
│       │   ├── pairwise/                     # anchor-based pairwise tables
│       │   └── filtered/                     # score table + linker + filters joined
│       ├── evals/
│       │   ├── evals_input_data.json         # which eval fields to merge
│       │   ├── variant_evals_all.parquet     # wide variant eval table
│       │   └── ensg_evals_all.parquet        # wide gene (ensg) eval table
│       └── full_analysis_tables/             # scores LEFT-JOINed onto evals (final)
├── pipeline_runs/                            # aggregation_pipeline_*.json audit trails
└── src/
    ├── config.json                          # orchestrator configuration
    ├── full_data_aggregation_pipeline.py    # runs the whole pipeline + audit
    ├── download_source_data.py
    ├── create_variant_scores_all_table.py
    ├── create_variant_eval_all_table.py
    ├── create_ensg_eval_all_table.py
    ├── create_percentile_score_tables.py
    ├── create_variant_scores_gene_aggregation.py
    ├── create_pairwise_score_tables.py
    ├── created_variant_scores_filtered_tables.py
    └── create_analysis_tables.py
```

## `full_data_aggregation_pipeline.py`

The single entry point for the whole build. It runs every other script in
`src/` **in dependency order** (as subprocesses, with the same Python
interpreter that launched it) so each script's inputs are already on disk by the
time it runs, then writes a full audit trail of the run. Each step is just a
normal CLI invocation of one of the scripts documented below.

Execution order (each step's key output in parentheses):

1. `download_source_data.py` (`raw_data/{scores,evals,linker,filters}`)
2. `create_variant_scores_all_table.py` (`scores/variant_scores_all_outer.parquet`)
3. `create_variant_eval_all_table.py` (`evals/variant_evals_all.parquet`)
4. `create_ensg_eval_all_table.py` (`evals/ensg_evals_all.parquet`)
5. `create_percentile_score_tables.py --variant` (the percentile + `all_inner` tables)
6. `create_variant_scores_gene_aggregation.py --gene-agg-stats` (`scores/gene_aggregated/`)
7. `create_pairwise_score_tables.py --variant` (`scores/pairwise/`)
8. `created_variant_scores_filtered_tables.py` (`scores/filtered/`)
9. `create_analysis_tables.py` (`full_analysis_tables/`)

Only the variant-level percentile / pairwise flavors are built, because the
score manifest declares no `gene_level` score sources (so there is no
`ensg_scores_all.parquet` to percentile or pair); the gene group of
`create_analysis_tables.py` is fed instead by the gene-aggregated tables from
step 6.

### Configuration (`config.json`)

Everything a user can tune lives in [`config.json`](#configjson) next to the
script, which the pipeline reads by default. Any CLI flag overrides the matching
config value; booleans accept both spellings (e.g. `--skip-download` /
`--no-skip-download`). Point at a different file with `--config PATH` or ignore
it entirely (built-in defaults) with `--no-config`.

### The audit trail

Every run writes `aggregation_pipeline_{timestamp}.json` under `--report-dir`
(default `../pipeline_runs`). It is a provenance/forensics record so that, if
anything looks wrong analytically, you can both (i) replay the exact command
that produced any table and (ii) see where a row count changed (i.e. where data
may have been lost or fanned out). The report is re-written after **every** step,
so even a crashed or interrupted run leaves a usable partial audit. It records,
per step: the `started_at` timestamp and duration; the exact `command` (argv +
copy-pasteable string); every input consumed (with creation time, size, and row
count) and every output written (with row count). A `files_written` index
inverts this into a per-output-file provenance map. The resolved `config` (and
its file path) is also stored under `invocation` for reproducibility.

Row counts come from Parquet footer metadata, so they are cheap and exact for
every Parquet table (a partitioned `*.parquet` directory counts as one table).
Counting the raw bgzipped/TSV *text* sources requires a full decompressed scan,
so it is **off by default**; enable it with `--count-text-inputs` (or disable
all counting with `--no-row-counts`).

### Usage

```bash
# From write_tables_pipeline/src
python full_data_aggregation_pipeline.py                     # uses config.json
python full_data_aggregation_pipeline.py --config my.json    # alternate config
python full_data_aggregation_pipeline.py --no-skip-download  # override a config bool
python full_data_aggregation_pipeline.py --memory-limit 20GB --threads 8
python full_data_aggregation_pipeline.py --from create_percentile_score_tables
python full_data_aggregation_pipeline.py --only create_analysis_tables
python full_data_aggregation_pipeline.py --list-steps
python full_data_aggregation_pipeline.py --dry-run           # print the plan only
```

### Behavior notes

- `--from` / `--to` run an inclusive range of steps (by name); `--only` runs a
  specific subset; `--skip` drops named steps. `--only` is mutually exclusive
  with `--from` / `--to`. `--list-steps` prints the resolved step order.
- `--skip-download` (default `true` in the shipped config) drops the download
  step when the source data is already local.
- `--overwrite` is forwarded only to steps that skip existing outputs by default
  (gene aggregation, pairwise, filtered, analysis); the score/eval/percentile
  steps always rewrite regardless.
- `--memory-limit` / `--threads` / `--temp-dir` are forwarded to every DuckDB
  step. By default the run stops at the first failing step; pass
  `--continue-on-error` to push on (downstream steps will likely error on
  missing inputs).

## `config.json`

The all-encompassing configuration for `full_data_aggregation_pipeline.py` — the
single place to control how the pipeline runs. It has two sections:

- **`pipeline`** — the run-wide options, mirroring the CLI flags above
  (`skip_download`, `overwrite`, `memory_limit`, `threads`, `temp_dir`,
  `report_dir`, `count_text_inputs`, `no_row_counts`, `continue_on_error`,
  `only`, `from_step`, `to_step`, `skip`). These seed the parser's defaults; a
  CLI flag overrides the matching value. Inline `_pipeline_docs` documents each.
- **`steps`** — the ordered list of script invocations. Each entry has a `name`
  (audit label + selection key), a `script` (file in this `src/` dir), `args`
  (extra CLI args passed to that script), `passthrough` (which pipeline options
  to forward: any of `memory_limit`, `threads`, `temp_dir`, `overwrite`), and an
  `enabled` toggle (set `false` to skip a step entirely).

The obs/exp exclusion (see `create_variant_eval_all_table.py` below) is dictated
here, by the `create_variant_eval_all_table` step's
`args: ["--exclude", "_obs_exp_mis"]`. That is its single source of truth: the
pipeline forwards it to the script **and** parses it to drop those sources from
the audit's recorded inputs. If `config.json` is missing or `--no-config` is
given, the pipeline falls back to conservative built-in defaults.

## `download_source_data.py`

Downloads the GCS objects listed in
`../data/raw_data/input_data_locations.json` into the local `raw_data` tree:

| Manifest key    | Destination                                                       |
| --------------- | ----------------------------------------------------------------- |
| `score_tables`  | `../data/raw_data/scores`                                         |
| `eval_tables`   | `../data/raw_data/evals`                                          |
| `linker_table`  | `../data/raw_data/linker`                                         |
| `filter_tables` | `../data/raw_data/filters/variant` and `../data/raw_data/filters/ensg` |

The `score_tables`/`eval_tables` groups are organized into subcategories
(`variant_level`, `gene_level`); `linker_table` is a single standalone URI (the
variant→gene linker). `filter_tables` is also subcategorized, but its two
subcategories route to *different* destinations: `variant_level` →
`../data/raw_data/filters/variant` and `gene_level` →
`../data/raw_data/filters/ensg`. Both single files (`.tsv`, `.tsv.bgz`,
`.tsv.gz`) and partitioned directory prefixes (`.parquet/`) are supported;
downloads run via `gcloud storage cp -r`.

### Prerequisites

- The [Google Cloud SDK](https://cloud.google.com/sdk) (`gcloud`) on your `PATH`.
- Authentication with access to the source buckets:

```bash
gcloud auth login
```

### Usage

```bash
# From write_tables_pipeline/src
python download_source_data.py                              # download everything (default)
python download_source_data.py --score-tables variant_level # only score variant_level
python download_source_data.py --eval-tables variant_level,gene_level
python download_source_data.py --score-tables all --eval-tables none
python download_source_data.py --linker-table all          # only the linker table
python download_source_data.py --filter-tables all         # only the filter tables
python download_source_data.py --dry-run                    # preview gcloud commands only
```

### Selection semantics

- With **no** selection flags, everything is downloaded (scores, evals, filters,
  linker).
- If **any** of `--score-tables` / `--eval-tables` / `--filter-tables` /
  `--linker-table` is supplied, only the explicitly requested categories are
  downloaded; the unspecified groups default to `none`.
- `--score-tables` / `--eval-tables` / `--filter-tables` accept a comma-separated
  list of subcategory keys (e.g. `variant_level,gene_level`) or the keywords
  `all` / `none`. Subcategory keys are read dynamically from the manifest.
- `--linker-table` accepts `all` / `none` (it is a single URI, not subcategorized).
- Both dashed and underscored spellings are accepted (e.g. `--score-tables`
  and `--score_tables`).

### Behavior notes

- Destination directories are created automatically if missing.
- Downloads continue past individual failures; failed URIs are listed in a
  summary and the script exits non-zero if any failed.
- An unknown subcategory key produces an error listing the valid options.

## `create_variant_scores_all_table.py`

Merges every variant-level score field into a single wide
`../data/processed_data/scores/variant_scores_all_outer.parquet`.

For each entry under `variant_level` in
`../data/processed_data/scores/score_input_data.json`, the script projects the
requested score column(s) and **full-outer-joins** all sources on the variant
key `(chrom, pos, ref, alt)` so that no variant from any score file is dropped.

In each `score_fields` mapping, the **key** is the source column name in the
input file and the **value** is the name the column should have in the output
(e.g. `{"AM_score": "AM"}` reads `AM_score` and writes it as `AM`).

### Engine / memory strategy

Each score table is ~78M rows and the full set is ~22 GB, so the merge is built
**one source at a time into an on-disk DuckDB table** rather than in one shot:

- Each source is first deduped to one row per variant key (`max()` per score,
  ignoring NULLs), then `FULL JOIN`-ed onto the growing wide table, which is
  re-materialized each step. Only one 2-table join runs at a time and the row
  count stays ~one-per-variant, so peak memory/disk stays bounded.
- This deliberately avoids two approaches that exhausted the disk (~74 GB
  spill) at this scale: a single 16-way join (pipelines every hash build at
  once) and `UNION ALL` + `GROUP BY` (inflates the input ~16x to ~1.25B rows).
- `.tsv.bgz` sources are read **directly** via streaming gzip decompression —
  no multi-GB decompressed copies are written to disk.
- Rows with an incomplete key (any NULL `chrom/pos/ref/alt`) are dropped.
- The intermediates live in an on-disk build DB under `--temp-dir`, removed on
  completion; the final table is streamed to Parquet via `COPY ... TO`.

For reference, this builds ~79M variants x 17 scores (~3.7 GB zstd parquet) in
roughly 5 min with `--memory-limit 20GB` on a 36 GB / 83 GB-free machine.
Runtime is intentionally traded for a bounded memory/disk footprint.

### Source handling

- Partitioned Parquet directories (`part-*.parquet` + `_SUCCESS`) are read via a
  `**/*.parquet` glob (the `_SUCCESS` marker is ignored).
- Bgzipped TSVs are read with quoting disabled, since columns such as
  `alleles`/`values` contain unescaped JSON; only the key + score columns are
  projected.
- The chromosome key is normalized across formats (`chr` → `chrom`).

### Usage

```bash
# From write_tables_pipeline/src (requires the score sources to be downloaded first)
python create_variant_scores_all_table.py
python create_variant_scores_all_table.py --memory-limit 8GB --threads 4
python create_variant_scores_all_table.py --output /tmp/scores.parquet
python create_variant_scores_all_table.py --dry-run   # print the SQL and exit

# Subset the manifest (mutually exclusive; match by file_path, basename, or substring)
python create_variant_scores_all_table.py --include AM_AM_score.parquet sift.parquet
python create_variant_scores_all_table.py --exclude coalesced_snp_misfit.tsv.bgz
```

### Behavior notes

- `--include` / `--exclude` select which manifest sources to merge. They are
  mutually exclusive; each token matches an entry by full `file_path`, basename,
  or substring, and a token that matches nothing is a hard error.
- The on-disk build DB + spill default to `<output dir>/.duckdb_build` and are
  removed on completion; point `--temp-dir` at a volume with enough free space.
- Duplicate output score names across the *selected* sources are rejected up front.
- A missing source, missing join key, or missing score column produces a clear
  error (listing the columns that are available).

## `create_variant_eval_all_table.py`

The eval counterpart of `create_variant_scores_all_table.py`: merges every
variant-level eval field from
`../data/processed_data/evals/evals_input_data.json` into a single wide
`../data/processed_data/evals/variant_evals_all.parquet` (full-outer on
`(chrom, pos, ref, alt)`; no variant dropped). Same engine and memory strategy
(sequential materialized full-outer joins, per-source dedup, on-disk build DB).

The one difference is **per-field typing**, since eval fields are heterogeneous:

- `is_pos`-style labels → normalized to `BOOLEAN` (handles `true`/`false`
  strings, native booleans, and `0`/`1`; deduped with `bool_or`). A field is
  treated as boolean if its output name starts with `is_pos` or its source
  column is BOOLEAN.
- `observed_*` / `expected_*` and other numerics → single-precision `FLOAT`
  (deduped with `max`). `FLOAT` rather than `DOUBLE` roughly halves the on-disk
  footprint of the numeric columns.

The per-source dedup is essential for the observed/expected tables: those are
keyed per *variant* but stored per *transcript*, so a raw outer join would fan
out up to ~23x per key. Collapsing each source to one row per key first (safe
since there are no within-key conflicts) prevents the blow-up.

**Key derivation.** Most sources carry `chrom/pos/ref/alt` directly. The
`hail_evaluation_tables` (e.g. `schema_no_dup`, `epi25_*`, `asc_*`,
`genebass_af_matched_controls`) instead expose Hail's `locus` (`chr1:925942`)
and `alleles` (`["A","G"]`). When the direct quartet is absent, the key is
derived from those: `chrom`/`pos` by splitting `locus` on `:`, and `ref`/`alt`
from the cleaned `alleles` array. The `chr` prefix is kept so derived keys match
sources that carry `chrom` directly.

### Usage

```bash
# From write_tables_pipeline/src (requires the eval sources to be downloaded first:
#   python download_source_data.py --eval-tables variant_level)
python create_variant_eval_all_table.py
python create_variant_eval_all_table.py --memory-limit 20GB
python create_variant_eval_all_table.py --dry-run   # print the SQL and exit

# Subset the manifest (mutually exclusive; match by file_path, basename, or substring)
python create_variant_eval_all_table.py --exclude _obs_exp_mis      # drop the 3 obs/exp sources
python create_variant_eval_all_table.py --include clinvar.tsv.bgz gnomad_independent_mis.tsv.bgz
```

> The current `variant_evals_all.parquet` is built **excluding** the three
> `*_obs_exp_mis.parquet` sources (`--exclude _obs_exp_mis`), since those inputs
> are wrong upstream — see
> `../.agent_reports/obs_exp_mis_preprocessing_assessment.md`. When building via
> the orchestrator this exclusion is set in
> [`config.json`](#configjson) (the `create_variant_eval_all_table` step's
> `args`), which is its single source of truth.

### Behavior notes

- `--include` / `--exclude` select which manifest sources to merge (same matching
  and mutual-exclusivity as the scores script above).

## `create_ensg_eval_all_table.py`

The gene-level counterpart of `create_variant_eval_all_table.py`: merges every
`gene_level` eval field from
`../data/processed_data/evals/evals_input_data.json` into a single wide
`../data/processed_data/evals/ensg_evals_all.parquet`, full-outer-joined on the
gene key `ensg` so no gene from any source is dropped. Same engine and memory
strategy (sequential materialized full-outer joins, per-source dedup to one row
per gene, on-disk build DB) — though gene-level sources are small (a handful of
MB, one row per gene), so in practice it is fast; the staged design is kept
purely for parity with the variant script.

Like the variant eval merge, fields are **typed per field**:

- `is_pos`-style labels → `BOOLEAN` (handles `true`/`false` strings, native
  booleans, and `0`/`1`); a field is boolean if its output name starts with
  `is_pos` or its source column is BOOLEAN. (The clinvar / genebass burden gene
  labels.)
- `n_case_*` / `n_ctrl_*` and any other numeric → single-precision `FLOAT`
  (parallel with the variant merge; `FLOAT` rather than `DOUBLE` roughly halves
  the on-disk footprint).

All current gene-level sources are plain TSVs (read with quoting disabled and
`all_varchar`, then cast explicitly); partitioned Parquet directories are
supported via a `**/*.parquet` glob.

### Usage

```bash
# From write_tables_pipeline/src (requires the gene-level eval sources to be downloaded first:
#   python download_source_data.py --eval-tables gene_level)
python create_ensg_eval_all_table.py
python create_ensg_eval_all_table.py --memory-limit 20GB
python create_ensg_eval_all_table.py --dry-run   # print the SQL and exit

# Subset the manifest (mutually exclusive; match by file_path, basename, or substring)
python create_ensg_eval_all_table.py --include dd-lof-vep asd-lof-vep
python create_ensg_eval_all_table.py --exclude schema_ensg_eval
```

### Behavior notes

- `--include` / `--exclude` select which manifest sources to merge (same matching
  and mutual-exclusivity as the scores/variant-eval scripts above).
- All current gene-level sources expose the Ensembl gene id directly as `ensg`;
  extra key aliases are accepted defensively for future sources.

## `create_percentile_score_tables.py`

Derives percentile-transformed tables from a wide `<base>_all_outer.parquet` score
table (built by `create_variant_scores_all_table.py`), writing them into the
same directory. For `variant_scores_all_outer.parquet` the base is `variant_scores`:

| Output | What it is |
| ------ | ---------- |
| `variant_scores_outer_pre_percentile.parquet`  | The original (outer) table, each score column replaced by its per-column percentile. NULL stays NULL; percentiles are computed over each column's non-NULL values independently. |
| `variant_scores_inner_pre_percentile.parquet`  | `..._outer_pre_percentile` with every row that has any NULL score dropped. Since percentiles preserve NULLs, this is the *outer* percentiles restricted to the fully-populated rows — equivalent to an INNER JOIN across the score columns. |
| `variant_scores_all_inner.parquet`             | The original table with every NA-score row dropped (an INNER JOIN / drop-NA of the input). A standalone deliverable that also serves as the intermediate for the `inner_post` table below. |
| `variant_scores_inner_post_percentile.parquet` | Built **from** `variant_scores_all_inner.parquet`: percentiles computed on the inner rows (so the ranking/denominator is the inner set, not the outer set). |
| `variant_scores_percentile_thresholds.tsv`     | Optional byproduct (on by default): for each score, the raw value at each requested percentile threshold. |

The key columns (`(chrom, pos, ref, alt)` for variant scores) are carried
through unchanged; every other column is treated as a score and percentiled.

In the percentile parquet outputs, each score column is renamed to make the
transform explicit: the pre-percentile tables (`outer_pre`, `inner_pre`) use
`<score>_pre_percentile` and the post table (`inner_post`) uses
`<score>_post_percentile` (e.g. `AM` -> `AM_pre_percentile` /
`AM_post_percentile`). `variant_scores_all_inner.parquet` and the threshold TSV
keep the original score names.

### Datasets (variant / ensg)

By default the script processes **both** `variant_scores_all_outer.parquet` and
`ensg_scores_all.parquet`. `--variant` / `--ensg` restrict it to one;
`--variant-score-path` / `--ensg-score-path` (also accepted as
`--variant_score_path` / `--ensg_score_path`) override the input locations. Each
requested dataset's input must exist or the script errors, so while the
gene-level `ensg_scores_all.parquet` is not yet built, pass `--variant`.

### Threshold TSV byproduct

For each score, the TSV lists the raw score value at each percentile threshold
(default `[0.85, 0.9, 0.95, 0.98, 0.99, 0.995]`) — the smallest raw value whose
percentile is ≥ the threshold, over that score's **full non-NULL distribution**
(the same basis as `outer_pre_percentile`). Override with `--tsv-thresholds`
(values in `(0, 1]`) or disable entirely with `--no-tsv`.

### Percentile definition (max tie-break)

For a column with `N` non-NULL values, a value `v` maps to
`pct(v) = (count of non-NULL values <= v) / N` — i.e. `cume_dist` with **max**
tie-breaking, so every member of a tie group gets the rank of the *last*
element. This guarantees the maximum value maps to exactly `1.0` even when it
is heavily duplicated (an "average" tie-break would pull a 4%-mass maximum down
to ~0.98).

### Engine / memory strategy

Same DuckDB out-of-core philosophy as the table it consumes. The percentile
transform is applied **one column at a time into an on-disk DuckDB table**: for
each score column a small distinct-value → percentile map is built (`GROUP BY`
to distinct values, then one ordered running `SUM` for the max-tie cumulative
fraction) and `LEFT JOIN`-ed back, re-materializing the wide table each step
and dropping the previous one. Only one column's ranking is in flight at a
time, so peak memory/disk is bounded; runtime is traded for that footprint.

### Usage

```bash
# From write_tables_pipeline/src (requires variant_scores_all_outer.parquet to exist first)
python create_percentile_score_tables.py                       # both datasets (ensg skipped)
python create_percentile_score_tables.py --variant             # variant only
python create_percentile_score_tables.py --memory-limit 20GB --threads 4
python create_percentile_score_tables.py --variant-score-path /tmp/variant_scores_all_outer.parquet
python create_percentile_score_tables.py --tsv-thresholds 0.9 0.99 0.999
python create_percentile_score_tables.py --no-tsv --dry-run    # print the plan and exit
```

### Behavior notes

- Score columns are discovered dynamically as all input columns except the keys.
- A score is treated as missing if it is NULL **or** NaN.
- The on-disk build DB + spill default to `<output dir>/.duckdb_build` and are
  removed on completion; point `--temp-dir` at a volume with enough free space.
- The run prints the row counts of all four parquet outputs and warns if
  `inner_pre`, `all_inner`, and `inner_post` disagree (they must all be the
  INNER set).

## `create_pairwise_score_tables.py`

Builds **anchor-based pairwise** score tables from the wide score parquet files,
writing one small file per pair into a `pairwise/` directory alongside
`variant_scores_all_outer.parquet`.

A *pairwise* operation pairs one **anchor** score with every other
(*non-anchor*) score and restricts to the rows where **both** scores are present
— the per-pair intersection. For anchor `polyphen` and 16 other scores this
yields 16 files (e.g. `AM_polyphen_anchor.parquet`), each with the key columns
plus the anchor and non-anchor columns. Because every pair's intersection is a
different subset, each file has its own row set.

Output files are named `{non_anchor}_{anchor}_anchor.parquet` (raw) and
`{non_anchor}_{anchor}_anchor_{pre,post}_percentile.parquet` (percentile flavors).

### Three flavors (subdirectories of `pairwise/`)

| Subdir | What each pair file holds | Input |
| ------ | ------------------------- | ----- |
| `pairwise_raw/`  | The **raw** anchor/non-anchor values on the pair's intersection — no percentile transform. | `<base>_all.parquet` |
| `pairwise_pre/`  | Percentiles computed **before** pairing: each column is its *outer* (whole-distribution) percentile, then the pair is intersected. No percentile is recomputed — the precomputed outer percentiles are just restricted to the intersection. | `<base>_outer_pre_percentile.parquet` |
| `pairwise_post/` | Percentiles computed **after** pairing: the pair is intersected first, then `CUME_DIST` is computed over that intersection. Recomputed from scratch because each pair's intersection differs from any precomputed table. | `<base>_all.parquet` |

The `pairwise_raw` files keep the original score names (`polyphen`, `AM`, …).
The percentile flavors rename each score column to record the transform:
`{score}_{anchor}_anchor_{pre,post}_percentile` (e.g.
`AM_polyphen_anchor_pre_percentile`), matching the column-renaming convention in
`create_percentile_score_tables.py`. The `pre` values are mathematically
identical to the corresponding non-pairwise outer percentiles; the rename is
purely for bookkeeping and downstream naming consistency.

For a given pair the three flavors share an identical row set (the pair's
intersection), so `raw`, `pre`, and `post` are three views of the same rows.

### Anchor

`--anchor` selects the anchor score; it defaults to `polyphen` when that column
is present in the dataset, otherwise it must be given explicitly. A file is
produced for every *other* (non-anchor) score.

### Datasets (variant / ensg)

By default the script processes **both** `variant_scores` and `ensg_scores`.
`--variant` / `--ensg` restrict it to one. Each requested dataset's inputs must
exist or the script errors, so while the gene-level tables are not yet built,
pass `--variant`. (The two datasets share the same `pairwise_*` subdirs; pick
distinct anchors per dataset to avoid filename collisions, or run them into
separate `--output-dir`s.)

### Percentile definition (max tie-break)

`pairwise_post` uses `CUME_DIST` over the pair's intersection — for `N` rows a
value `v` maps to `(count of values <= v) / N` (max tie-break, so the maximum is
exactly `1.0`). `pairwise_pre` inherits the same convention from the precomputed
outer-percentile table.

### Engine / memory strategy

Same DuckDB out-of-core philosophy as the tables it consumes. Each pair is a
single streaming `COPY (...) TO` over a `read_parquet` scan; `pairwise_post`
additionally sorts the pair's intersection for `CUME_DIST` (DuckDB spills to
`--temp-dir` under `--memory-limit`). Only one pair is in flight at a time, and
existing output files are **skipped** (resumable) unless `--overwrite` is set.

### Usage

```bash
# From write_tables_pipeline/src (requires variant_scores_all_outer.parquet and, for pairwise_pre,
# variant_scores_outer_pre_percentile.parquet to exist first)
python create_pairwise_score_tables.py --variant                    # raw + pre + post
python create_pairwise_score_tables.py --variant --anchor cadd      # different anchor
python create_pairwise_score_tables.py --variant --modes raw post   # subset of flavors
python create_pairwise_score_tables.py --variant --memory-limit 20GB --threads 4
python create_pairwise_score_tables.py --variant --dry-run          # print the plan and exit
```

### Behavior notes

- Score columns are discovered dynamically as all input columns except the keys.
  `--anchor` is always an **original** score name (e.g. `polyphen`), even for
  `pre`, whose input columns carry the `_pre_percentile` suffix from
  `create_percentile_score_tables.py` — the script maps names accordingly.
- For `raw`/`post`, a score is treated as present only if it is not NULL **and**
  not NaN; for `pre`, the precomputed percentile is NULL exactly where the raw
  value was, so the same intersection results.
- Each requested mode validates its own input exists, with a hint pointing at
  the script that builds it.
- Outputs default to `../data/processed_data/scores/pairwise/{pairwise_raw,
  pairwise_pre,pairwise_post}/`; override the base with `--output-dir`.
- The DuckDB spill dir defaults to `<output dir>/.duckdb_spill` and is removed on
  completion; point `--temp-dir` at a volume with enough free space.

## `create_variant_scores_gene_aggregation.py`

First step of the **gene-aggregation** feature: attaches a gene (`ensg`) label
to each non-percentile variant score table by INNER-JOINing it with the
variant→gene linker (`../data/raw_data/linker/linker_all.parquet`, downloaded by
`download_source_data.py`). Outputs land in
`../data/processed_data/scores/gene_aggregated/` as `{name}_ensg.parquet`.

The linker carries `(chrom, pos, ref, alt, ensg)`. Because a variant can map to
more than one gene, the INNER JOIN may **duplicate** score rows — one copy per
gene a variant belongs to. This is intended: downstream steps aggregate per
`ensg`, and a variant contributing to two genes is fine. Variants absent from
the linker are dropped (INNER JOIN), as required for a per-gene aggregation.

Each output schema is the variant key, then `ensg`, then every score column from
the input, in that order:

```
chrom, pos, ref, alt, ensg, <score_1>, <score_2>, ...
```

### Inputs

By default both non-percentile variant score tables are processed:

| Input | Built by |
| ----- | -------- |
| `variant_scores_all_outer.parquet`  | `create_variant_scores_all_table.py` |
| `variant_scores_all_inner.parquet`  | `create_percentile_score_tables.py`  |

The percentile tables (`*_pre_percentile` / `*_post_percentile`) are
intentionally excluded. Pass `--input` to gene-annotate a specific table (or set
of tables) instead; a bare filename is resolved against the scores directory.

### Engine / memory strategy

Same DuckDB out-of-core philosophy as the tables it consumes. Each output is a
single streaming `COPY (...) TO` over a 2-table `read_parquet` INNER JOIN
(DuckDB builds one hash table on the small 5-column linker and probes it with the
score table, spilling to `--temp-dir` under `--memory-limit`). Only one join is
in flight at a time.

### Usage

```bash
# From write_tables_pipeline/src (requires the score tables and the linker to exist first:
#   python download_source_data.py --linker-table all)
python create_variant_scores_gene_aggregation.py
python create_variant_scores_gene_aggregation.py --input variant_scores_all_outer.parquet
python create_variant_scores_gene_aggregation.py --gene-agg-stats               # also write *_ensg_stats.parquet
python create_variant_scores_gene_aggregation.py --gene-agg-stats --memory-limit 20GB --threads 4
python create_variant_scores_gene_aggregation.py --dry-run   # print the plan + SQL and exit
```

### Per-gene statistics (`--gene-agg-stats`)

With `--gene-agg-stats`, after each `{name}_ensg.parquet` is written the script
also computes per-`ensg` statistics for **every** score column and writes a
second file `{name}_ensg_stats.parquet` — the original columns plus the statistic
columns appended (named `{score}_{statistic}`). The table keeps one row per
`(variant, gene)`.

For each score `c` and gene `g` (only the non-NULL/non-NaN values of `c` within
`g` contribute):

| Column | Meaning |
| ------ | ------- |
| `{c}_mean`   | Mean of the gene's values (broadcast to every row of the gene). |
| `{c}_median` | Median (continuous/interpolated; broadcast). |
| `{c}_max`    | Maximum (broadcast). |
| `{c}_p90` / `{c}_p95` | 90th / 95th percentile score, using the **same discrete max-tie convention** as `create_percentile_score_tables.py` (smallest value `v` with `(count ≤ v)/N ≥ p`; broadcast). |
| `{c}_zscore` | **Per-row** standard score `(x − mean_g) / std_g`, with `std_g` the **population** standard deviation. NULL when the gene has zero variance or this row's score is missing. |
| `{c}_modified_zscore` | **Per-row** modified z-score `0.6745·(x − median_g) / MAD_g`, where `MAD_g = median(\|x_i − median_g\|)` (Iglewicz–Hoaglin / statology). NULL when `MAD_g = 0` or this row's score is missing. |

`std` and `MAD` are intermediates and are not emitted. mean/median/max/p90/p95
are gene-level properties, so they are populated for every row of a gene
(including rows whose own score is NULL).

All appended statistic columns are stored as single-precision `FLOAT` (computed
in DOUBLE, narrowed only for storage); these ~7-per-score columns are otherwise
high-entropy and compress poorly, so this roughly halves the output size. The
raw score columns are themselves already `FLOAT` (carried through from the wide
score table). For reference, the 17-score variant tables produce 141-column
outputs of ~14 GB (outer) and ~9.5 GB (inner) — sizes that predate the
DOUBLE → FLOAT narrowing and will be smaller once the tables are rebuilt.

#### Stats engine / memory strategy

The statistics step favors memory/disk safety. Per-gene statistics are computed
**one score column at a time** into a tiny per-gene dimension table (~one row per
gene), reading only `ensg` + that one score column from the columnar Parquet each
pass (holistic median/quantile aggregates spill to `--temp-dir`). The final
`{name}_ensg_stats.parquet` is then a single streaming `COPY` of the `_ensg`
table LEFT-JOINed to that tiny dimension table, computing the per-row z-scores on
the fly — no wide ~80M-row intermediate is ever materialized.

### Behavior notes

- Score columns are discovered dynamically as all input columns except the keys
  (`chrom, pos, ref, alt`); a missing key column or a pre-existing `ensg` column
  is a hard error.
- Existing outputs are **skipped** (resumable) unless `--overwrite` is set; this
  applies to both `{name}_ensg.parquet` and `{name}_ensg_stats.parquet`.
- The on-disk build DB + DuckDB spill default to `<output dir>/.duckdb_build` and
  are removed on completion; point `--temp-dir` at a volume with enough free space.
- The run prints input vs. output row counts (output ≥ input due to multi-gene
  duplication, minus any variants absent from the linker) and, for stats, the
  final row × column counts.

## `created_variant_scores_filtered_tables.py`

Attaches the variant→gene linker **and every filter table** to a wide variant
score table via a long chain of **FULL OUTER JOIN**s, writing
`../data/processed_data/scores/filtered/variant_scores_outer_pre_percentile_filtered.parquet`.

By default the input is `variant_scores_outer_pre_percentile.parquet`. The joins
run in this order:

1. **Linker** (`../data/raw_data/linker/linker_all.parquet`) on
   `(chrom, pos, ref, alt)`, in the *same manner* as
   `create_variant_scores_gene_aggregation.py` — it adds `ensg` and fans the
   table out to one row per `(variant, gene)` (full-outer here, so no variant is
   dropped).
2. **Every variant-level filter** under `../data/raw_data/filters/variant`. Each
   `variant_filters_*.tsv.gz` file *is* a set of variants — membership in the
   file is the filter — so each contributes a single **BOOLEAN** membership
   column named after the file (`variant_filters_<name>.tsv.gz` → `<name>`):
   `TRUE` when the variant is in that filter, `FALSE` otherwise. The filter is
   reduced to its distinct variant keys first, so it cannot fan the join out.
3. **The gene/ensg-level filter** (`../data/raw_data/filters/ensg/ensg_filters.tsv`)
   on `ensg`. Its non-key columns are carried through with an `ensg_` prefix
   (e.g. `CATH_class_1` → `ensg_CATH_class_1`) to avoid colliding
   (case-insensitively) with the variant-level filter names; `uniprot_id` stays
   a string, every other column is cast to `BOOLEAN`.

The result has the variant key, `ensg`, every original score, one boolean per
variant-level filter (~65), and the prefixed gene-level filter columns. For the
current inputs this is ~80.5M `(variant, gene)` rows × 155 columns (~5.4 GB
zstd).

### Engine / memory / disk strategy

65+ sequential joins over an ~80M-row (post fan-out) table is extremely taxing,
but every right-hand input here is **small**: each variant filter is reduced to
its distinct keys, the gene filter is ~one row per gene, and only the linker is a
large (single) hash build. So a batch of these joins can be **pipelined and
streamed in one pass** without ever materializing the wide table:

- Each batch of `--checkpoint-every` joins is nested into a single
  `COPY (...) TO` that streams straight from the previous checkpoint Parquet to
  the next. The wide table is **never stored in a DuckDB table** — that is what
  exhausts disk, since DuckDB's single build file does not release dropped-stage
  space mid-connection. Small build sides stay in memory; large spills go to
  `--temp-dir` under `--memory-limit`.
- The script **checkpoints to Parquet** every `--checkpoint-every` joins (fresh
  connection per batch) and keeps **at most two checkpoints on disk** (the new
  one is written and row-count-verified before the previous is deleted), so disk
  use is bounded to roughly two copies of the compressed wide table.
- A crash only loses the joins since the last checkpoint: on restart the script
  finds the latest checkpoint and **resumes** from the next join. A `plan.json`
  records the exact join order; a resume whose plan no longer matches the inputs
  is refused (re-run with `--overwrite` to start clean).

### Usage

```bash
# From write_tables_pipeline/src (requires the input score table, the linker, and the
# filters to exist locally first):
python download_source_data.py --linker-table all --filter-tables all
python create_percentile_score_tables.py --variant   # builds the default input

python created_variant_scores_filtered_tables.py
python created_variant_scores_filtered_tables.py --checkpoint-every 5 --memory-limit 24GB
python created_variant_scores_filtered_tables.py --join-type left   # keep only the score table's variants
python created_variant_scores_filtered_tables.py --input variant_scores_all_outer.parquet --output /tmp/out.parquet
python created_variant_scores_filtered_tables.py --dry-run          # print the join plan + sample SQL
```

### Behavior notes

- `--join-type` chooses the outer-join flavor: `full` (default) keeps variants
  present only in a filter or the linker (their score/`ensg` columns are NULL);
  `left` keeps only the input score table's variants.
- Variant filter membership columns are `FALSE` for non-members and `TRUE` for
  members; they are discovered from the filenames and a (case-insensitive) name
  collision between two filters is a hard error.
- If the output already exists the run is a no-op unless `--overwrite` is set
  (which also clears any checkpoints). Use `--keep-intermediates` to retain the
  checkpoint directory after a successful run.
- The DuckDB spill dir defaults to `<output dir>/.duckdb_spill` and the
  checkpoints to `<output dir>/.ckpt_<output stem>/`; both are removed on a
  successful, non-`--keep-intermediates` run. Point `--temp-dir` at a volume
  with enough free space for the linker step's spill.

## `create_analysis_tables.py`

The final assembly step: LEFT-JOINs most of the score tables under
`../data/processed_data/scores` onto the merged eval tables under
`../data/processed_data/evals`, writing the results into
`../data/processed_data/full_analysis_tables`. Every output is named
`{score file stem}_eval.parquet`.

Three groups of score tables are processed (select a subset with `--groups`):

| Group | Inputs | Joined onto | Output location |
| ----- | ------ | ----------- | --------------- |
| `variant`  | each top-level `scores/variant_scores_*.parquet` (the `_all_*` and `_*_percentile` wide tables) | `evals/variant_evals_all.parquet` on `(chrom, pos, ref, alt)` | flat in `full_analysis_tables/` |
| `pairwise` | every file under `scores/pairwise/{pairwise_raw,pairwise_pre,pairwise_post}` | `variant_evals_all.parquet` on the variant key | `full_analysis_tables/pairwise/<flavor>/` |
| `gene`     | the two `scores/gene_aggregated/*_ensg_stats.parquet` tables | `variant_evals_all.parquet` (variant key) **and** `evals/ensg_evals_all.parquet` (on `ensg`) | `full_analysis_tables/gene_aggregated/` |

The `filtered` score tables are intentionally **not** processed.

### Why LEFT JOIN is safe

Both eval tables are built with exactly **one row per key** — one per
`(chrom, pos, ref, alt)` for the variant evals and one per `ensg` for the gene
evals. Joining a unique-keyed table onto a score table therefore cannot
duplicate ("fan out") any score row, and the LEFT direction keeps every score
row (eval columns are simply NULL where there is no match). As a guard, each
output's row count is compared to its input's; any increase (which would mean an
eval table was not unique on its key) is reported loudly.

### Gene-aggregated dedup output

A gene-aggregated table has one row per `(variant, gene)`, so a multi-gene
variant appears in several rows sharing the same `(chrom, pos, ref, alt)`. For
each gene table two outputs are written: `{stem}_eval.parquet` (the join as-is)
and `{stem}_eval_deduped.parquet`, which keeps only rows whose
`(chrom, pos, ref, alt)` is **unique** in the table (any key occurring more than
once has *all* its rows dropped — an unambiguous resolution since there is no
principled single gene-copy to keep). The dedup is computed from the just-written
`_eval` file via a per-key count + semi-join (scanning only the four key columns),
so the join is not recomputed.

### Engine / memory strategy

Same DuckDB out-of-core philosophy as the upstream scripts. Each `_eval` output
is a single streaming `COPY (...) TO` over a `read_parquet` LEFT JOIN (the small
eval side is the hash build, the score table probes it, spilling to `--temp-dir`
under `--memory-limit`). The dedup step materializes the set of unique keys into
a TEMP table, then writes the deduped file with a streaming semi-join. Existing
outputs are skipped (resumable) unless `--overwrite` is set.

### Usage

```bash
# From write_tables_pipeline/src (requires the score tables and both eval tables to exist first)
python create_analysis_tables.py                          # all three groups
python create_analysis_tables.py --groups variant pairwise
python create_analysis_tables.py --groups gene --memory-limit 20GB
python create_analysis_tables.py --overwrite
python create_analysis_tables.py --dry-run                # print the plan and exit
```

### Behavior notes

- `--groups` chooses any subset of `variant` / `pairwise` / `gene` (default: all).
- `--scores-dir` / `--output-dir` override the input/output roots;
  `--variant-evals` / `--ensg-evals` override the eval table locations.
- Existing outputs are skipped unless `--overwrite` is set.
- `--compression` (default `zstd`) and `--row-group-size` (default 512000) tune
  the Parquet writer; `--memory-limit` / `--threads` / `--temp-dir` tune DuckDB.
