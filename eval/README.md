# BioStat-CLI

Memory-efficient CLI for genomic statistics on merged Parquet tables using Polars lazy execution.

## Features

- Computes `auc`, `auprc`, `enrichment`, `rate_ratio`, `pairwise_enrichment`, and `pairwise_rate_ratio`
- **Gene-averaged statistics**: `gene_avg_enrichment`, `gene_avg_rate_ratio`, `gene_avg_auc`, `gene_avg_auprc` — computes per-gene stat values and averages them, giving each gene equal weight
- Optional nonparametric bootstrap stderr (`std_error`) for all supported stats
- Supports `variant` and `gene` eval levels (strategy-based evaluators)
- Reads local paths and `gs://` parquet inputs via Polars/fsspec/gcsfs
- Writes:
  - main metrics TSV
  - run log JSON (args, resolved table path, total runtime, per eval/filter runtime)
- Optional missing-variant TSV to explain `rows_used` vs `total_eval_rows`
- Optional per-gene variant coverage TSV (`--write-gene-variant-coverage`)
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
- `--stat` `all` or csv subset (`auc,auprc,enrichment,rate_ratio,pairwise_enrichment,pairwise_rate_ratio,gene_avg_enrichment,gene_avg_rate_ratio,gene_avg_auc,gene_avg_auprc`)
- `--eval-set` optional csv eval override (defaults to all `evals` from resources)
- `--filters` optional csv logical filter names (from `Filters` keys); `none` is always included
- `--thresholds` optional csv thresholds
- `--case-total-by-eval` optional per-eval case totals (`eval_name:value,eval2:value2`); required for rate-ratio stats unless totals appear under `Case_totals` / `Ctrl_totals` in resources JSON
- `--ctrl-total-by-eval` optional per-eval control totals (`eval_name:value,eval2:value2`); same resolution as case totals
- `--bootstrap [N]` enable nonparametric row bootstrap stderr calculation; optional `N` sets sample count (e.g., `--bootstrap 50`, default `100` when `N` omitted)
- `--pvalue-method` p-value calculation method: `fisher` (default) or `poisson`. Fisher's exact test is recommended for 2×2 contingency tables; Poisson is the legacy approximation
- `--vsm-comparison-method` method for `vsm_comparison` pairwise VSM table: `fisher` (default) or `poisson`
- `--gene-col` gene identifier column for gene-averaged stats (default: `ensg`)
- `--write-gene-variant-coverage` write per-gene variant coverage report as a separate TSV (default: off)
- `--out-fname` output naming schema/prefix (**required**)
- `--write-missing` controls missing-entity report: `none`, `all`, or `any` (default: `none`)

Rate-ratio denominator resolution priority (high to low):

1. per-eval CLI (`--case-total-by-eval`, `--ctrl-total-by-eval`)
2. per-eval table metadata (`Case_totals`, `Ctrl_totals` in resources JSON)

Bootstrap behavior:

- Point `value` and `p_value` are computed on the original dataset.
- `std_error` precedence:
  - if `--bootstrap` is enabled: bootstrap SD from resampled `value`s (rows sampled with replacement within each eval/filter subset)
  - else if `stat=enrichment`: analytic $\mathrm{SE}(\ln \mathrm{LR}^{+})$ (see below); `p_value` still follows `--pvalue-method` (`fisher` or `poisson`)
  - else if `stat=rate_ratio` and `--pvalue-method poisson`: analytic Poisson stderr on RR scale
  - else: `NaN`
- **Enrichment** (`value` is unchanged: case rate / control rate). Analytic uncertainty uses the positive likelihood ratio $\mathrm{LR}^{+}=\frac{\mathrm{TP}/(\mathrm{TP}+\mathrm{FN})}{\mathrm{FP}/(\mathrm{FP}+\mathrm{TN})}$ (same ratio as `value` when defined).
  - $\mathrm{SE}(\ln \mathrm{LR}^{+})=\sqrt{\left(\frac{1}{\mathrm{TP}}-\frac{1}{\mathrm{TP}+\mathrm{FN}}\right)+\left(\frac{1}{\mathrm{FP}}-\frac{1}{\mathrm{FP}+\mathrm{TN}}\right)}$ on **raw** counts when this is defined.
  - If that expression is undefined on raw counts (e.g. $\mathrm{TP}=0$ or $\mathrm{FP}=0$), apply **+0.5 to all four cells** once and recompute $\mathrm{SE}$ and CI from the corrected table (stderr/CI only; `value` remains from raw counts).
  - **95% CI on the `value` (LR+) scale:** $\exp\left(\ln(\mathrm{LR}^{+})\pm 1.96\cdot \mathrm{SE}(\ln \mathrm{LR}^{+})\right)$ using the same cell table as for $\mathrm{SE}$.
  - Main TSV columns `enrichment_ci_lower` / `enrichment_ci_upper` store these bounds; they are set to `NaN` when `--bootstrap` is used (analytic CI not reported alongside bootstrap `std_error`).
- **Rate ratio** `value` is $\mathrm{TP}/\mathrm{case\_total}$ divided by $\mathrm{FP}/\mathrm{ctrl\_total}$ (when totals are set).
  - Analytic Poisson `std_error` (RR scale) when `--pvalue-method poisson`: $\mathrm{SE}[\log(\mathrm{RR})]=\sqrt{1/\mathrm{TP}+1/\mathrm{FP}}$, `std_error = RR \cdot \mathrm{SE}[\log(\mathrm{RR})]`, `NaN` when `TP == 0` or `FP == 0` (or RR undefined). With `fisher`, that analytic `std_error` is `NaN` but **p-value** is still Fisher.
  - **95% Wald CI on the RR scale** (same $\mathrm{SE}[\log(\mathrm{RR})]$ as above): $\exp(\log(\mathrm{RR})\pm 1.96\cdot \mathrm{SE}[\log(\mathrm{RR})])$ when $\mathrm{TP}>0$, $\mathrm{FP}>0$, and $\mathrm{RR}>0$; stored as `rate_ratio_ci_lower` / `rate_ratio_ci_upper`. These are cleared to `NaN` under `--bootstrap` (same as enrichment CIs).
- `auc`: `p_value` is a two-sided test of $H_0:\mathrm{AUC}=0.5$ using the Hanley–McNeil variance for the AUC estimate and a normal approximation (`z=(\mathrm{AUC}-0.5)/\mathrm{SE}`). `auprc` still has `p_value = NaN`.
- `--bootstrap N` must use `N >= 2` when bootstrap is enabled.

Output paths are derived from `--out-fname`:

- main TSV: `<schema>.tsv`
- log JSON: `<schema>_log.json`
- missing TSV: `<schema>_missing.tsv` (when `--write-missing` is `all` or `any`)
- gene variant coverage TSV: `<schema>_gene_variant_coverage.tsv` (when `--write-gene-variant-coverage` is set)

## Threshold behavior

- Thresholds are percentile-based **fractions in `[0,1]`**
- A row counts as **above** threshold when `score >= t` (ties at `t` are included).
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
- `enrichment_ci_lower`, `enrichment_ci_upper` (95% CI on enrichment ratio; `NaN` except for `stat=enrichment` without `--bootstrap`, and always `NaN` for other stats)
- `rate_ratio_ci_lower`, `rate_ratio_ci_upper` (95% Wald CI on rate ratio; `NaN` except for `stat=rate_ratio` when analytic CI is defined and without `--bootstrap`, and always `NaN` for other stats)
- `p_value`
- `tp`, `fp`, `tn`, `fn`
- `rows_used`
- `total_eval_rows`

For gene-averaged stats (`gene_avg_enrichment`, `gene_avg_rate_ratio`, `gene_avg_auc`, `gene_avg_auprc`), additional columns:

- `n_genes_used` — number of genes with valid per-gene values
- `n_genes_excluded` — number of genes excluded from the average

For pairwise stats (`pairwise_enrichment`, `pairwise_rate_ratio`), additional columns:

- `anchor_value` - baseline value from anchor VSM on full set
- `adjustment_ratio` - ratio of VSM performance to anchor performance on pairwise intersection

For `pairwise_auc`, `p_value` is a **paired DeLong** two-sided test of whether the anchor and VSM AUCs differ on the **pairwise intersection** cohort (Sun–Xu fast DeLong covariance; same binary labels, two score vectors). The anchor-only baseline row keeps `p_value = NaN`. `pairwise_auprc` still uses `p_value = NaN`.

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
  --case-total-by-eval "is_pos_dd:1000" \
  --ctrl-total-by-eval "is_pos_dd:5000" \
  --out-fname ../results/VSM_pairwise
```

## Gene-averaged statistics

Gene-averaged stats compute each metric per gene and then average across genes, giving each gene equal weight regardless of variant count (macro-averaging).

### Available stats

- `gene_avg_enrichment` — per-gene enrichment averaged; null = 1.0
- `gene_avg_rate_ratio` — per-gene rate ratio averaged (cohort totals per eval via `--case-total-by-eval` / `--ctrl-total-by-eval` or resources JSON); null = 1.0
- `gene_avg_auc` — per-gene AUC averaged; null = 0.5
- `gene_avg_auprc` — per-gene AUPRC averaged

### Requirements

- `--eval-level variant` (not compatible with gene eval level)
- A gene identifier column must be present in the parquet (default: `ensg`, override with `--gene-col`)

### Output columns

Gene-averaged rows in the main TSV include two additional columns:

- `n_genes_used` — genes with a computable (non-NaN) per-gene value
- `n_genes_excluded` — genes excluded (insufficient labels, zero denominators, etc.)
- `tp/fp/tn/fn` are `NaN` (no single contingency table)

### Inference

- `std_error` = SEM across per-gene values: $\mathrm{SD}/\sqrt{n}$
- `p_value` = one-sample t-test against the null (enrichment/RR: 1.0; AUC: 0.5; AUPRC: NaN)

### Gene variant coverage report

When `--write-gene-variant-coverage` is enabled, a separate TSV is written to `<prefix>_gene_variant_coverage.tsv` with columns:

- `eval_name`, `filter_name`, `score_name`, `gene`
- `n_variants_used` — variants with non-null score in this gene
- `n_variants_excluded` — variants with null score in this gene
- `n_variants_total` — sum of the above

This is computed per (eval, filter, score) combination and is independent of thresholds.

### Example usage

```bash
python -m biostat_cli.cli \
  --resources-json ../files/vsm_all.json \
  --table-name VSM_all_inner_per_v1 \
  --eval-level variant \
  --stat "gene_avg_enrichment,gene_avg_auc" \
  --thresholds 0.90,0.95,0.98,0.99 \
  --gene-col ensg \
  --write-gene-variant-coverage \
  --out-fname ../results/VSM_gene_avg
```

## Pairwise VSM comparison (exact Poisson design options)

When adding Poisson-based pairwise VSM comparison (`p_greater`, `p_less`) for
`rate_ratio` interpretation, there are multiple valid exact formulations.
For this project, the default should be `full_2x2_exact_poisson`.

`lock-exact-formula` means selecting exactly one formulation before coding so
implementation, tests, and interpretation all match.

### Candidate formulas (default first)

- `full_2x2_exact_poisson` (default):
  - Use all four counts `(TP_i, FP_i, TP_j, FP_j)` in a single exact test
    targeting the relative rate-ratio contrast.
- `tp_only_conditional_binomial`:
  - Use `X = TP_i`, `n = TP_i + TP_j`, null `p0 = 0.5`.
  - `p_greater = P(X >= TP_i | Binomial(n, p0))`
  - `p_less = P(X <= TP_i | Binomial(n, p0))`
- `fp_only_conditional_binomial`:
  - Analogous conditional-binomial test using FP counts.

### Input-count definitions (used by all options)

- Thresholding uses `score >= t` (ties at the cutoff count as above threshold).
- For boolean evals:
  - `TP = count(above_threshold and eval == True)`
  - `FP = count(above_threshold and eval == False)`
- In gene `sum_variants` mode:
  - `TP = sum(n_case)` above threshold
  - `FP = sum(n_ctrl)` above threshold

### Important caveat

Current `vsm_comparison` uses per-model score-non-null row sets, so model `i`
and model `j` counts can come from different subsets when missingness differs.
This is a marginal comparison of model-level contingency summaries, not a
strict paired-intersection test.

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
  --pvalue-method fisher \
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
