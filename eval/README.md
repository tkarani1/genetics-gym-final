# BioStat-CLI

Command-line toolkit for **Genetics Gym** evaluation tables: enrichment, rate ratios, ROC/PR metrics, pairwise-adjusted scores, gene-level summaries, and publication figures — all on large Parquet files via Polars lazy execution.

**Location in the repo:** `genetics-gym-final/eval/`

---

## What it does

You point the CLI at a merged evaluation Parquet (local path or `gs://`) and a small JSON config that names score columns, boolean eval labels, and optional filters. It writes TSV metrics plus a JSON run log. Typical workflows:

| Goal | Entry point |
|------|-------------|
| One-off metrics on a table | `biostat-cli` (or `python -m biostat_cli.cli`) |
| Many eval/filter combos in parallel | `python -m biostat_cli.cli_parallel` |
| Paper Figure 1 panels (compute + QC + plots) | `figure1-pipeline` |

---

## Quick start

### 1. Install

Requires **Python ≥ 3.10**.

```bash
cd genetics-gym-final/eval
pip install -e .
```

This registers two commands: `biostat-cli` and `figure1-pipeline`.

### 2. Configure a table

Copy or edit `resources.json`. Each table under `Table_info` needs at least `Path`, `Level`, `Score_cols`, and usually `evals`:

```json
{
  "Table_info": {
    "my_variant_table": {
      "Path": "/path/to/merged.parquet",
      "Level": "variant",
      "Score_cols": ["AM_percentile", "mpc_score_percentile"],
      "Filters": { "ordered": "filter_ordered" },
      "evals": ["is_pos_schema", "is_pos_dd"],
      "Case_totals": { "is_pos_schema": 1000 },
      "Ctrl_totals": { "is_pos_schema": 5000 }
    }
  }
}
```

Keys are case-flexible where noted in code (`evals` / `Evals`, `Case_totals` / `case_totals`).

### 3. Run

```bash
biostat-cli \
  --resources-json resources.json \
  --table-name my_variant_table \
  --eval-level variant \
  --stat enrichment,auc \
  --thresholds 0.90,0.95,0.98,0.99 \
  --out-fname results/my_run
```

Outputs (prefix = `--out-fname` without extension):

| File | When |
|------|------|
| `results/my_run.tsv` | Always |
| `results/my_run_log.json` | Always |
| `results/my_run_curves.tsv` | Continuous stats requested (`auc`, `auprc`, …) |
| `results/my_run_vsm_comparison.tsv` | `--stat` includes `vsm_comparison` |
| `results/my_run_missing.tsv` | `--write-missing all` or `any` |
| `results/my_run_gene_variant_coverage.tsv` | `--write-gene-variant-coverage` |

---

## Core concepts

### Evaluation level (`--eval-level`)

- **`variant`** — one row per variant; boolean eval column marks case vs control (or use observed/expected columns; see below).
- **`gene`** — one row per gene. Either a boolean gene label **or** weighted burden via `n_case` / `n_ctrl` when no `evals` are listed in JSON (`sum_variants` mode). See [GENE_EVAL_EXPLANATION.md](GENE_EVAL_EXPLANATION.md).

### Thresholds

Thresholds are **percentile fractions in `[0, 1]`**, not 0–100. A variant/gene counts as “above” when `score >= t` (ties included).

- Default if omitted: `0.90, 0.95, 0.98, 0.99, 0.995`
- Any value `> 1.0` exits with code **22**

### Rate-ratio denominators

`rate_ratio`, `pairwise_rate_ratio`, and `gene_avg_rate_ratio` need cohort sizes per eval. Resolution order:

1. CLI: `--case-total-by-eval` / `--ctrl-total-by-eval` (`eval_name:value,...`)
2. Resources JSON: `Case_totals` / `Ctrl_totals` on the table entry

### Input data

- Parquet with score columns (often precomputed percentiles) and eval/filter columns.
- **Pairwise-adjusted** stats expect precomputed anchor/VSM percentile columns (patterns below).
- **Observed/expected** stats expect `{eval_stem}_observed` and `{eval_stem}_expected` for that eval stem.

---

## Statistics reference

Use `--stat all` or a comma-separated subset. Names are case-insensitive.

### Variant / gene (standard)

| Stat | Description |
|------|-------------|
| `enrichment` | Case rate / control rate above threshold (2×2 from labels) |
| `rate_ratio` | (TP / case_total) / (FP / ctrl_total) |
| `auc`, `auprc` | ROC / PR area |
| `tpr_at_threshold`, `fpr_at_threshold`, `precision_at_threshold`, `recall_at_threshold` | Point on ROC/PR at each threshold |
| `auc_trunc`, `auprc_trunc` | Same, restricted to variants with score ≥ threshold |
| `obs_exp_ratio` | Σ observed / Σ expected above threshold (requires `{stem}_observed`, `{stem}_expected`) |
| `vsm_comparison` | All-pairs comparison of score columns (separate TSV; see below) |

### Pairwise-adjusted (precomputed columns)

| Stat | Description |
|------|-------------|
| `pairwise_enrichment`, `pairwise_rate_ratio` | Anchor-adjusted enrichment / RR |
| `pairwise_auc`, `pairwise_auprc` | AUC / AUPRC on pairwise intersection cohort |
| `pairwise_tpr_at_threshold`, … | Threshold metrics on pairwise scores |
| `pairwise_auc_trunc`, `pairwise_auprc_trunc` | Truncated continuous metrics |
| `pairwise_obs_exp_ratio` | Adjusted O/E ratio (same column rules as `obs_exp_ratio`) |

`pairwise_auc` uses a **paired DeLong** test vs the anchor on the intersection (anchor row: `p_value = NaN`). Other pairwise continuous stats keep `p_value = NaN`.

### Gene-averaged (macro over genes)

Requires `--eval-level variant` and a gene ID column (`--gene-col`, default `ensg`). **Not available in `cli_parallel`.**

| Stat | Null for t-test |
|------|-----------------|
| `gene_avg_enrichment`, `gene_avg_rate_ratio` | 1.0 |
| `gene_avg_auc` | 0.5 |
| `gene_avg_auprc` | (p-value NaN) |
| `gene_avg_tpr_at_threshold`, `gene_avg_fpr_at_threshold`, … | Same structure as variant-level counterparts |

Extra output columns: `n_genes_used`, `n_genes_excluded`. Contingency `tp/fp/tn/fn` are `NaN`.

---

## CLI reference (`biostat-cli`)

### Required

| Flag | Description |
|------|-------------|
| `--table-name` | Key under `Table_info` in resources JSON |
| `--eval-level` | `variant` or `gene` |
| `--out-fname` | Output prefix (see [Output paths](#output-paths)) |

### Common options

| Flag | Default | Description |
|------|---------|-------------|
| `--resources-json` | `resources.json` | Table metadata and paths |
| `--stat` | `all` | Comma-separated stat names or `all` |
| `--eval-set` | all table `evals` | Comma-separated eval columns |
| `--filters` | all table `Filters` | Comma-separated filter names; `none` always runs |
| `--thresholds` | `0.90,…,0.995` | Comma-separated percentile cutoffs |
| `--case-total-by-eval` | — | `eval:value,...` for rate ratios |
| `--ctrl-total-by-eval` | — | Same for controls |
| `--bootstrap [N]` | off | Nonparametric row bootstrap for `std_error` (default **N=100**; need **N ≥ 2**) |
| `--pvalue-method` | `fisher` | `fisher` or `poisson` for enrichment / rate ratio |
| `--vsm-comparison-method` | `fisher` | `fisher` or `poisson` for `vsm_comparison` |
| `--within-gene-percentile` | off | Rank scores within `ensg` before thresholding (variant level only) |
| `--gene-col` | `ensg` | Gene ID for gene-averaged stats |
| `--write-gene-variant-coverage` | off | Per-gene variant coverage TSV |
| `--chromosomes` | — | e.g. `chr1,chrX,chrM` (column must be `chrom` or `CHROM`) |
| `--write-missing` | `none` | `none`, `all`, or `any` missing-score report |

Point estimates (`value`, `p_value`) are always computed on the full data. With `--bootstrap`, `std_error` is the SD across bootstrap replicates; analytic confidence intervals for enrichment/rate ratio are set to `NaN`.

### Uncertainty without bootstrap

- **`enrichment`**: analytic SE on ln(LR⁺) and 95% CI on the ratio scale (`enrichment_ci_lower` / `enrichment_ci_upper`); +0.5 smoothing only when raw-cell formula is undefined.
- **`rate_ratio`**: Wald CI on RR scale when Poisson log-RR SE is defined (`rate_ratio_ci_*`); with `fisher`, analytic `std_error` is `NaN` but Fisher p-value still computed.
- **`auc`**: two-sided test of AUC = 0.5 (Hanley–McNeil SE).
- **`gene_avg_*`**: SEM across genes; one-sample t-test vs null above.

Details: [biostat_cli/BIOSTAT_PVALUE_STDERR_SLIDES.md](biostat_cli/BIOSTAT_PVALUE_STDERR_SLIDES.md).

---

## Pairwise column layout

Pairwise stats auto-detect columns from names (no extra JSON):

| Pattern | Role |
|---------|------|
| `{anchor}_anchor_percentile` | Anchor on full variant set S* |
| `{vsm}_percentile_with_anchor` | VSM on S_i ∩ S* |
| `{anchor}_anchor_percentile_with_{vsm_short}` | Anchor on S_i ∩ S* |

Example for anchor `mpc_score` and VSM `esm1b_score`:

```
mpc_score_anchor_percentile
esm1b_score_percentile_with_anchor
mpc_score_anchor_percentile_with_esm1b
```

Adjusted enrichment:

$$
\mathrm{Enr}(\mathrm{VSM}_i) = \mathrm{Enr}(\mathrm{VSM}^*, S^* \cap S_e) \times \frac{\mathrm{Enr}(\mathrm{VSM}_i, S_i \cap S^* \cap S_e)}{\mathrm{Enr}(\mathrm{VSM}^*, S_i \cap S^* \cap S_e)}
$$

```bash
biostat-cli \
  --resources-json resources.json \
  --table-name VSM_pairwise_table \
  --eval-level variant \
  --stat pairwise_enrichment,pairwise_rate_ratio \
  --thresholds 0.90,0.95,0.98,0.99 \
  --case-total-by-eval "is_pos_dd:1000" \
  --ctrl-total-by-eval "is_pos_dd:5000" \
  --out-fname results/pairwise_run
```

---

## Observed / expected evals

If the table has `{eval_stem}_observed` and `{eval_stem}_expected` for an eval stem, requesting `obs_exp_ratio` or `pairwise_obs_exp_ratio` computes

$$
\text{O/E} = \frac{\sum \text{observed}}{\sum \text{expected}}
$$

over rows with `score >= threshold`. Boolean TP/FP/TN/FN stats are not mixed into the same eval row set.

---

## VSM comparison output

`--stat vsm_comparison` compares **every pair of score columns** at each threshold on the **intersection of non-null scores** for that pair. Results go to `<prefix>_vsm_comparison.tsv`, not the main TSV.

Columns include: `eval_name`, `filter_name`, `vsm_i`, `vsm_j`, `threshold`, `odds_ratio`, `p_greater`, `p_less`, `log_odds_ratio`, `standard_error`, confidence interval bounds, `rows_used_pair`.

Methods (`--vsm-comparison-method`):

- **`fisher`** (default) — Fisher exact on the 2×2 table built from each model’s TP/FP at the threshold.
- **`poisson`** — exact Poisson formulation on the same contingency summaries.

This is a **marginal** comparison of per-model summaries (row sets can differ when missingness differs), not a strict paired variant-level test.

---

## Main output TSV schema

| Column | Notes |
|--------|--------|
| `eval_name`, `filter_name`, `score_name`, `threshold`, `stat` | Keys |
| `value`, `std_error`, `p_value` | Estimate and uncertainty |
| `enrichment_ci_lower`, `enrichment_ci_upper` | Enrichment only; `NaN` under bootstrap |
| `rate_ratio_ci_lower`, `rate_ratio_ci_upper` | Rate ratio only; `NaN` under bootstrap |
| `tp`, `fp`, `tn`, `fn` | `NaN` for gene-averaged stats |
| `rows_used`, `total_eval_rows` | Rows contributing to the stat |
| `rows_retained`, `n_pos_retained`, `n_neg_retained` | Truncated / threshold-point stats |
| `anchor_value`, `adjustment_ratio` | Pairwise enrichment / RR only |
| `n_genes_used`, `n_genes_excluded` | Gene-averaged stats only |

**Curves** (`<prefix>_curves.tsv`): `curve_type` (`roc` / `pr`), `point_idx`, `score_threshold`, `fpr`, `tpr`, `precision`, `recall`.

**Missing report** (`<prefix>_missing.tsv`): optional; sorted by eval, filter, category, then locus. ID columns auto-detected (variant: `chrom,pos,ref,alt` or `locus,alleles`; gene: `ensg`, `gene_symbol`, etc.).

**Run log** (`<prefix>_log.json`): `run_args`, `table_path`, `output_files`, `elapsed_seconds`, `eval_filter_elapsed_seconds`.

---

## Parallel runner

```bash
python -m biostat_cli.cli_parallel \
  --resources-json resources.json \
  --table-name my_variant_table \
  --eval-level variant \
  --stat enrichment,auc \
  --out-fname results/my_run_parallel
```

Runs eval/filter combinations concurrently (same outputs as the serial CLI for supported flags).

**Use the serial `biostat-cli` when you need:**

- `--bootstrap`
- Gene-averaged stats (`gene_avg_*`)
- `--write-gene-variant-coverage`
- `--gene-col` override (parallel uses default evaluator paths only for gene-averaged — gene_avg not supported anyway)

---

## Figure 1 pipeline

End-to-end workflow: raw + pairwise metrics → panel table → QC → plots. Config (`figure1_pipeline_config.json`) holds **panel layout and column names**, not Parquet paths — paths are passed on the CLI.

```bash
figure1-pipeline run \
  --config figure1_pipeline_config.json \
  --mode both \
  --raw-parquet /path/to/vsm_all_figure1_ready.parquet \
  --pairwise-parquet /path/to/vsm_pairwise_anchor_figure1_ready.parquet \
  --outdir results/figure1_run
```

### Subcommands

| Command | Purpose |
|---------|---------|
| `run` | Compute, QC, and plot |
| `compute` | Metrics + panel table + QC only |
| `plot` | Plot from existing `panel_table.tsv` or metrics TSVs |

### Useful flags

| Flag | Description |
|------|-------------|
| `--mode` | `raw`, `pairwise`, or `both` |
| `--profile` | `paper_figure1` (fixed panels) or `all_variant` (auto `is_pos_*` evals) |
| `--output-layout` | `combined`, `per_eval`, or `both` |
| `--threshold` | Panel/plot threshold (e.g. `0.95`) |
| `--thresholds` | Override compute thresholds |
| `--bootstrap N` | Pass through to biostat for `std_error` |
| `--eval-set` | Override eval list |
| `--paper-strict` | Fail if required paper evals are missing |
| `--dry-run` | Print plan + validate parquets |
| `--overwrite` | Replace existing outdir |

### Typical outputs in `--outdir`

- `figure1_run_resources.json` — ephemeral resources pointing at your parquets
- `metrics_raw.tsv`, `metrics_pairwise.tsv`
- `panel_table.tsv`, `qc_summary.tsv`, `qc_report.md`
- `run_manifest.json`
- `figure1_raw.png/pdf`, `figure1_pairwise.png/pdf`, `figure1_combined.png/pdf` (depending on `--mode`)

Other preset configs in this directory: `figure1_pipeline_config_*.json`, `figure2_pipeline_config*.json`.

---

## Examples

### Enrichment + bootstrap

```bash
biostat-cli \
  --resources-json resources.json \
  --table-name VSM_all_inner_per_v1 \
  --eval-level variant \
  --stat enrichment \
  --thresholds 0.95 \
  --bootstrap 20 \
  --out-fname results/enrichment_boot20
```

### Gene-averaged + coverage report

```bash
biostat-cli \
  --resources-json resources.json \
  --table-name VSM_all_inner_per_v1 \
  --eval-level variant \
  --stat gene_avg_enrichment,gene_avg_auc \
  --gene-col ensg \
  --write-gene-variant-coverage \
  --out-fname results/gene_avg
```

### Within-gene percentiles

```bash
biostat-cli \
  --resources-json resources.json \
  --table-name VSM_all_inner_per_v1 \
  --eval-level variant \
  --stat enrichment \
  --within-gene-percentile \
  --out-fname results/within_gene
```

### Chromosome subset

```bash
biostat-cli \
  --resources-json resources.json \
  --table-name VSM_all_inner_per_v1 \
  --eval-level variant \
  --stat auc \
  --chromosomes chr1,chr2,chrX \
  --out-fname results/chr_subset
```

---

## Project layout

```
eval/
├── README.md                 # This file
├── pyproject.toml            # Package biostat-cli 0.2.0
├── resources.json            # Example table registry
├── biostat_cli/              # Library + CLI
│   ├── cli.py                # Serial runner (full feature set)
│   ├── cli_parallel.py       # Parallel eval/filter runner
│   ├── pipeline/             # figure1-pipeline
│   └── stats/                # Stat implementations
├── figure1_pipeline_config.json
├── GENE_EVAL_EXPLANATION.md  # Gene-level eval modes
└── tests/                    # pytest suite
```

Repo-specific driver scripts (e.g. `run_gene_eval.py`, `run_gsm_missense_pairwise_gene_avg.py`) wrap `biostat_cli` for particular Genetics Gym tables; treat them as examples once you understand the flags above.

---

## Development

```bash
cd genetics-gym-final/eval
pip install -e ".[dev]"
pytest
```

When changing CLI behavior, update this README in the same change (see `.cursor/rules/readme-sync.mdc`).

---

## Further reading

- [GENE_EVAL_EXPLANATION.md](GENE_EVAL_EXPLANATION.md) — boolean vs `sum_variants` gene evals
- [biostat_cli/BIOSTAT_PVALUE_STDERR_SLIDES.md](biostat_cli/BIOSTAT_PVALUE_STDERR_SLIDES.md) — p-values, bootstrap, and confidence intervals
