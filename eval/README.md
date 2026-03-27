# BioStat-CLI

Memory-efficient CLI for genomic statistics on merged Parquet tables using Polars lazy execution.

## Features

- Computes `auc`, `auprc`, `enrichment`, `rate_ratio`, `pairwise_enrichment`, and `pairwise_rate_ratio`
- Optional nonparametric bootstrap stderr (`std_error`) for all supported stats
- Supports `variant` and `gene` eval levels (strategy-based evaluators)
- Reads local paths and `gs://` parquet inputs via Polars/fsspec/gcsfs
- Writes:
  - main metrics TSV
  - run log JSON (args, resolved table path, total runtime, per eval/filter runtime)
- Optional missing-variant TSV to explain `rows_used` vs `total_eval_rows`
- Includes a parallel runner (`biostat_cli.cli_parallel`) for eval/filter concurrency

## Install

```bash
cd /Users/tk508/Work/new/gg-script-codex
pip install -e .
```

## Resources JSON format

`--resources-json` defaults to `resources.json`.

```json
{
  "Table_info": {
    "VSM_all_inner_per_v1": {
      "Path": "gs://bucket/path/to/data.parquet",
      "Level": "variant",
      "Score_cols": ["AM_percentile", "score_PAI3D_percentile"],
      "Filters": {"ordered": "filter_ordered"},
      "evals": ["is_pos__schema", "is_pos_dd"],
      "Case_totals": {
        "is_pos__schema": 1000,
        "is_pos_dd": 1200
      },
      "Ctrl_totals": {
        "is_pos__schema": 5000,
        "is_pos_dd": 5200
      }
    }
  }
}
```

## CLI arguments

- `--resources-json` path to resources file (default: `resources.json`)
- `--table-name` table key under `Table_info` (**required**)
- `--eval-level` `variant` or `gene` (**required**)
- `--stat` `all` or csv subset (`auc,auprc,enrichment,rate_ratio,pairwise_enrichment,pairwise_rate_ratio`)
- `--eval-set` optional csv eval override (defaults to all `evals` from resources)
- `--filters` optional csv logical filter names (from `Filters` keys); `none` is always included
- `--thresholds` optional csv thresholds
- `--case-total`, `--ctrl-total` optional denominators for rate ratio
- `--case-total-by-eval` optional per-eval case totals (`eval_name:value,eval2:value2`)
- `--ctrl-total-by-eval` optional per-eval control totals (`eval_name:value,eval2:value2`)
- `--bootstrap [N]` enable nonparametric row bootstrap stderr calculation; optional `N` sets sample count (e.g., `--bootstrap 50`, default `100` when `N` omitted)
- `--out-fname` output naming schema/prefix (**required**)
- `--write-missing` controls missing-entity report: `none`, `all`, or `any` (default: `none`)

Rate-ratio denominator resolution priority (high to low):

1. per-eval CLI overrides (`--case-total-by-eval`, `--ctrl-total-by-eval`)
2. per-eval table metadata (`Case_totals`, `Ctrl_totals` in resources JSON)
3. global CLI totals (`--case-total`, `--ctrl-total`)

Bootstrap behavior:

- Point `value` and `p_value` are computed on the original dataset.
- `std_error` is estimated from bootstrap resamples (rows sampled with replacement within each eval/filter subset).
- `--bootstrap N` must use `N >= 2` when bootstrap is enabled.

Output paths are derived from `--out-fname`:

- main TSV: `<schema>.tsv`
- log JSON: `<schema>_log.json`
- missing TSV: `<schema>_missing.tsv` (when `--write-missing` is `all` or `any`)

## Threshold behavior

- Thresholds are percentile-based **fractions in `[0,1]`**
- Default thresholds: `0.90,0.95,0.98,0.99`
- Passing any threshold `> 1.0` exits with error code `22`

Example:

```bash
--thresholds 0.90,0.95,0.98,0.99
```

## Main output TSV schema

Columns:

- `eval_name`
- `filter_name`
- `score_name`
- `threshold`
- `stat`
- `value`
- `std_error`
- `p_value`
- `tp`, `fp`, `tn`, `fn`
- `rows_used`
- `total_eval_rows`

For pairwise stats (`pairwise_enrichment`, `pairwise_rate_ratio`), additional columns:

- `anchor_value` - baseline value from anchor VSM on full set
- `adjustment_ratio` - ratio of VSM performance to anchor performance on pairwise intersection

## Pairwise statistics

Pairwise statistics compute adjusted enrichment/rate_ratio that maximize variant coverage per VSM comparison.

### Formula

```
enr(VSM_i) = enr(VSM*, S* ∩ S_e) × [enr(VSM_i, S_i ∩ S* ∩ S_e) / enr(VSM*, S_i ∩ S* ∩ S_e)]
```

Where:
- `VSM*` is the anchor VSM (the one with maximum variant coverage)
- `S*` is the set of variants defined by the anchor
- `S_i` is the set of variants defined by VSM_i
- `S_e` is the evaluation set (after filtering)

### Required input table format

Pairwise statistics require pre-computed percentile columns with specific naming:

| Column Pattern | Description |
|----------------|-------------|
| `{anchor}_anchor_percentile` | Anchor percentile on full set S* |
| `{vsm}_percentile_with_anchor` | VSM_i percentile on S_i ∩ S* |
| `{anchor}_anchor_percentile_with_{vsm_short}` | Anchor percentile on S_i ∩ S* |

Example columns for anchor `mpc_score` with VSMs `esm1b_score` and `MisFit_S_score`:

```
mpc_score_anchor_percentile
esm1b_score_percentile_with_anchor
mpc_score_anchor_percentile_with_esm1b
MisFit_S_score_percentile_with_anchor
mpc_score_anchor_percentile_with_MisFit_S
```

The column structure is auto-detected; no additional JSON configuration required.

### Example usage

```bash
python -m biostat_cli.cli \
  --resources-json ../files/vsm_all.json \
  --table-name VSM_pairwise_table \
  --eval-level variant \
  --stat "pairwise_enrichment,pairwise_rate_ratio" \
  --thresholds 0.90,0.95,0.98,0.99 \
  --case-total 1000 \
  --ctrl-total 5000 \
  --out-fname ../results/VSM_pairwise
```

## Missing-output TSV (optional)

Use to inspect score-null entities after eval/filter masking.

- One row per entity (per eval/filter)
- `all` mode: only variants missing score values in all methods
- `any` mode: variants missing in one or more methods, with category for all vs partial
- Missing report sort order:
  - first by `eval_name`
  - then by `filter_name` (when present)
  - then by `missing_category` (`all_methods` before `partial_methods`)
  - then variant-level chromosome/position (`chr1`..`chr22`, `chrX`, `chrY`) or gene identifier/name

Columns:

- `eval_name`
- `filter_name`
- id columns are auto-detected for both levels:
  - variant: `chrom,pos,ref,alt` or `locus,alleles`
  - gene: `GENE_ID`, `gene_id`, `ensg`, or `gene_symbol`
- `missing_category` (`all_methods` or `partial_methods`)
- `missing_score_count`
- `missing_score_names`

## Runtime logging

The generated log JSON (`<schema>_log.json`) includes:

- `run_args`
- `table_path`
- `output_files`
- `elapsed_seconds`
- `eval_filter_elapsed_seconds` (per eval/filter runtime)

## Examples

### Standard runner

```bash
python -m biostat_cli.cli \
  --resources-json ../files/vsm_all.json \
  --table-name VSM_all_inner_per_v1 \
  --eval-level variant \
  --stat "enrichment,auc" \
  --thresholds 0.90,0.95,0.98,0.99 \
  --out-fname ../results/VSM_v1 \
  --write-missing any
```

### Parallel runner

```bash
python -m biostat_cli.cli_parallel \
  --resources-json ../files/vsm_all.json \
  --table-name VSM_all_inner_per_v1 \
  --eval-level variant \
  --stat "enrichment,auc" \
  --out-fname ../results/VSM_v1_parallel
```

## Figure 1 Pipeline (local-data-first)

This repository now includes a one-command Figure 1-style pipeline that:

- computes raw metrics (`enrichment`, `rate_ratio`)
- computes pairwise-adjusted metrics (`pairwise_enrichment`, `pairwise_rate_ratio`)
- builds a panel-ready table (default threshold `0.95`)
- writes QC outputs
- renders `raw`, `pairwise`, and `combined` plots

Default config file:

- `figure1_pipeline_config.json` (panels, evals, score column names — **not** parquet paths)

Parquet inputs are **required on the CLI** for `run` and `compute` so you can swap tables without editing the preset.

### Quick start

```bash
figure1-pipeline run \
  --config figure1_pipeline_config.json \
  --mode both \
  --raw-parquet /path/to/vsm_all_figure1_ready.parquet \
  --pairwise-parquet /path/to/vsm_pairwise_mpc_anchor_figure1_ready.parquet
```

### Other modes

```bash
# Compute TSV/QC only (no plotting)
figure1-pipeline compute \
  --config figure1_pipeline_config.json \
  --mode both \
  --raw-parquet /path/to/vsm_all_figure1_ready.parquet \
  --pairwise-parquet /path/to/vsm_pairwise_mpc_anchor_figure1_ready.parquet

# Plot from an existing panel table
figure1-pipeline plot --config figure1_pipeline_config.json --mode both --panel-table /path/to/panel_table.tsv
```

### Helpful options

- `--raw-parquet` / `--pairwise-parquet` — required for `run`/`compute` when the mode includes that family (`raw`, `pairwise`, or `both`)
- `--outdir /path/to/output_dir` to choose output location
- `--threshold 0.95` to select panel-plot threshold
- `--thresholds 0.90,0.95,0.98,0.99` to override compute thresholds
- `--dry-run` to preview what will run (for `run`/`compute`, also runs preflight validation on the parquets)

### Expected outputs

Each run directory contains:

- `figure1_run_resources.json` (ephemeral resources file written for this run; points at the CLI parquets)
- `metrics_raw.tsv` (when mode includes `raw`)
- `metrics_pairwise.tsv` (when mode includes `pairwise`)
- `panel_table.tsv`
- `qc_summary.tsv`
- `qc_report.md`
- `run_manifest.json`
- `figure1_raw.png/pdf` (when mode includes `raw`)
- `figure1_pairwise.png/pdf` (when mode includes `pairwise`)
- `figure1_combined.png/pdf` (when mode is `both`)
