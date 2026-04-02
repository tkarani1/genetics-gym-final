# Gene-Level Evaluation: Statistical Tests

## Overview

Gene evaluation is activated with `--eval-level gene` and supports two modes depending on the eval column structure.

## Two Modes of Gene Evaluation

### Mode 1: Boolean gene eval

The input table has a boolean eval column per gene (e.g., `is_pos_schema`). Each gene is one observation — positive or negative.

### Mode 2: `sum_variants` (weighted gene-burden)

No boolean eval column. Instead, each gene row has pre-computed **`n_case`** and **`n_ctrl`** columns (variant counts). Selected automatically when the resources JSON provides no `evals`.

---

## Step 1: Contingency Table Construction

For a given score threshold `t`:

### Boolean mode

| | score > t | score ≤ t |
|---|---|---|
| **eval = True** | TP = count of positive genes above | FN = count of positive genes below |
| **eval = False** | FP = count of negative genes above | TN = count of negative genes below |

Every gene contributes exactly **1** to one cell.

### `sum_variants` mode

| | score > t | score ≤ t |
|---|---|---|
| **case variants** | TP = Σ `n_case` for genes above | FN = Σ `n_case` for genes below |
| **control variants** | FP = Σ `n_ctrl` for genes above | TN = Σ `n_ctrl` for genes below |

A gene with `n_case=50, n_ctrl=200` scoring above threshold contributes **50** to TP and **200** to FP. Cells are variant-weighted.

---

## Step 2: Enrichment

Same formula in both modes — only TP/FP/TN/FN values differ.

```
case_rate  = TP / (TP + FN)
ctrl_rate  = FP / (FP + TN)
enrichment = case_rate / ctrl_rate
```

- **Boolean mode**: fraction of positive genes above threshold ÷ fraction of negative genes above threshold. "Are positive genes disproportionately high-scoring?"
- **`sum_variants` mode**: fraction of total case variants in high-scoring genes ÷ fraction of total control variants in those genes. "Do high-scoring genes carry a disproportionate share of case variants?"

---

## Step 3: Rate Ratio

Uses externally-supplied denominators instead of in-table totals.

```
case_rate  = TP / case_total
ctrl_rate  = FP / ctrl_total
rate_ratio = case_rate / ctrl_rate
```

`case_total` and `ctrl_total` come from CLI args or the resources JSON — they represent total population sizes (e.g., total de novo variants in cases/controls across all genes).

- **Boolean mode**: rate = (positive genes above threshold) / (total positive genes in population).
- **`sum_variants` mode**: rate = (case variants in genes above threshold) / (total case variants in cohort).

---

## Step 4: P-value

Both enrichment and rate ratio use the same p-value test, selectable via `--pvalue-method` (`fisher` default, `poisson` legacy). Variables:

```
n1 = TP + FP    (total variants/genes above threshold)
k1 = TP         (case variants/positive genes above threshold)
n2 = FN + TN    (total variants/genes below threshold)
k2 = FN         (case variants/positive genes below threshold)
```

Calculation:

```
null_rate = k2 / n2          (case rate in low-scoring genes)
expected  = null_rate × n1   (expected cases in high-scoring genes under null)
p_value   = P(X ≥ k1)       where X ~ Poisson(expected)
```

**Poisson method** (`--pvalue-method poisson`): treats `k1` as Poisson-distributed. Good approximation when the above-threshold group is small.

**Fisher method** (`--pvalue-method fisher`, default): Fisher's exact test on the 2×2 contingency table (`alternative='greater'`). Exact, makes no distributional assumption.

**Null hypothesis** (both): the case variant rate (or positive gene rate) in high-scoring genes equals the rate in low-scoring genes. If k1 significantly exceeds the expected count, the p-value is small.

---

## Summary: What Differs Between Modes

| | Boolean mode | `sum_variants` mode |
|---|---|---|
| Unit of observation | 1 gene = 1 count | 1 gene = `n_case` + `n_ctrl` counts |
| TP meaning | # positive genes above threshold | # case variants in genes above threshold |
| Enrichment asks | Are positive genes enriched above threshold? | Are case variants enriched in high-scoring genes? |
| P-value null | Gene-level positive rate is uniform | Variant-level case rate is uniform across score strata |
| AUC / AUPRC | Available | Disabled (no binary labels) |

---

## Available Stats by Mode

| Stat | Boolean gene eval | `sum_variants` mode |
|---|:-:|:-:|
| Enrichment | ✓ | ✓ (weighted) |
| Rate ratio | ✓ | ✓ (weighted) |
| P-value (fisher/poisson) | ✓ | ✓ |
| AUC | ✓ | ✗ |
| AUPRC | ✓ | ✗ |
| Pairwise enrichment | ✓ (if columns exist) | ✓ (if columns exist) |
| Pairwise rate ratio | ✓ (if columns exist) | ✓ (if columns exist) |
