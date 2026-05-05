---
marp: true
theme: default
paginate: true
header: "biostat_cli — p-values & uncertainty"
footer: "GenGym / biostat_cli"
style: |
  section { font-size: 28px; }
  section.lead h1 { font-size: 48px; }
  .small { font-size: 22px; }
  table { font-size: 24px; }
---

<!-- _class: lead -->

# biostat_cli
## P-values, standard errors & bootstrap

Variant-level enrichment / rate ratio at percentile thresholds

---

## Setup: contingency at threshold $t$

- **Case** = `is_pos == true`, **control** = `is_pos == false`
- **Above** = score $\ge t$ (percentile threshold $t \in [0,1]$)

$$
\begin{bmatrix}
\mathrm{TP} & \mathrm{FP}\\
\mathrm{FN} & \mathrm{TN}
\end{bmatrix}
$$

Rows: above / below $t$. Columns: case / control.

---

## Point estimates (context)

**Enrichment**

$$
\mathrm{Enr}=\frac{\mathrm{TP}/(\mathrm{TP}+\mathrm{FN})}{\mathrm{FP}/(\mathrm{FP}+\mathrm{TN})}
$$

**Rate ratio** (with per-eval cohort totals $N_1$, $N_2$ from CLI `--*-total-by-eval` or resources JSON)

$$
\mathrm{RR}=\frac{\mathrm{TP}/N_1}{\mathrm{FP}/N_2}
$$

---

## Fisher — p-value

$2\times2$ table as above. **Fisher’s exact test**, `alternative="greater"` (one-sided).

Tests whether **cases are over-represented** among scores above $t$ vs independence.

Same contingency drives **enrichment** and **rate_ratio** $p$-values in this mode.

---

## Fisher — `std_error` in output

**Enrichment:** reported `std_error` is **not** SE of $\mathrm{Enr}$ on the ratio scale.

It is the **Wald SE of $\log \mathrm{OR}$** with **Haldane–Anscombe +0.5** on each cell:

$$
\mathrm{SE}_{\log\mathrm{OR}}=\sqrt{\sum_{c\in\{\mathrm{TP},\ldots\}}\frac{1}{c+0.5}}
$$

**Rate ratio:** `std_error` = **`nan`** when $p$-value method is Fisher.

---

## Poisson — p-value (approximate one-sided test)

Among rows **below** $t$: $\pi_{\text{below}}=\dfrac{\mathrm{FN}}{\mathrm{FN}+\mathrm{TN}}$.

Expected **case count above** $t$ if that rate applied to all **above** rows ($m=\mathrm{TP}+\mathrm{FP}$):

$$
\lambda=\pi_{\text{below}}\cdot m
$$

$$
p=\Pr\{X\ge \mathrm{TP}\},\quad X\sim\mathrm{Poisson}(\lambda)
$$

---

## Poisson — `std_error`

**Enrichment:** **`nan`** (no paired analytic SE in code for Enr in Poisson mode).

**Rate ratio:** delta method on $\log \mathrm{RR}$:

$$
\mathrm{SE}_{\log\mathrm{RR}}=\sqrt{\frac{1}{\mathrm{TP}}+\frac{1}{\mathrm{FP}}},\qquad
\mathrm{SE}_{\mathrm{RR}}\approx \mathrm{RR}\cdot \mathrm{SE}_{\log\mathrm{RR}}
$$

(requires $\mathrm{TP},\mathrm{FP}>0$ and finite $\mathrm{RR}$).

---

## Bootstrap (`--bootstrap N`)

1. Resample **rows with replacement** ($N$ times), same table height.
2. Recompute requested stats each replicate.

**Unchanged on full data:** `value`, `p_value` (still Fisher or Poisson from **original** table).

**Overwritten:** `std_error` ← sample SD of replicate **`value`**s ($n-1$ denominator).

---

## Bootstrap — formula

For each output row (same eval / filter / score / threshold / stat):

$$
\bar\theta=\frac{1}{N}\sum_{b=1}^{N}\hat\theta_b,\qquad
s^2=\frac{1}{N-1}\sum_{b=1}^{N}(\hat\theta_b-\bar\theta)^2
$$

`std_error` $=s$, or `nan` if $N<2$ / too few finite replicates.

---

## Side-by-side

| | **Fisher** | **Poisson** |
|---|------------|-------------|
| **$p$-value** | Exact Fisher, one-sided greater | Poisson tail vs $\lambda$ from below-$t$ rate |
| **Enr `std_error`** | $\mathrm{SE}(\log\mathrm{OR})$, +0.5 cells | `nan` |
| **RR `std_error`** | `nan` | $\mathrm{RR}\cdot\sqrt{1/\mathrm{TP}+1/\mathrm{FP}}$ |

**Bootstrap:** empirical SD of $\hat\theta$; does **not** replace $p$.

---

## CLI reminder

- `--pvalue-method fisher` | `poisson`
- `--bootstrap` or `--bootstrap N` (default $N{=}100$ in standalone CLI)

`figure1-pipeline` passes `--pvalue-method` and `--bootstrap` into the same `biostat_cli.run` logic.

---

## Caveats (one slide)

- Fisher vs Poisson changes **definition** of the contingency $p$-value, not the point estimate formula for Enr/RR.
- Enrichment Fisher `std_error` is on **log-odds** scale, not $\mathrm{Enr}$ itself.
- Bootstrap SD reflects **row resampling** of the prepared frame, not genes/variants as clusters unless the table already encodes that structure.

---

<!-- _class: lead -->

# Questions

`biostat_cli/stats/binary.py` · `biostat_cli/cli.py`

