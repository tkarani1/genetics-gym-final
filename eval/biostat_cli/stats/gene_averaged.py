from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import ttest_1samp


@dataclass(frozen=True)
class GeneAvgResult:
    value: float
    p_value: float
    std_error: float
    n_genes_used: int
    n_genes_excluded: int


def _aggregate(values: list[float], n_total_genes: int, null_value: float) -> GeneAvgResult:
    clean = np.array([v for v in values if not (isinstance(v, float) and math.isnan(v))], dtype=np.float64)
    n_used = int(clean.size)
    n_excluded = n_total_genes - n_used
    if n_used == 0:
        return GeneAvgResult(
            value=math.nan, p_value=math.nan, std_error=math.nan,
            n_genes_used=0, n_genes_excluded=n_total_genes,
        )
    mean_val = float(np.mean(clean))
    if n_used == 1:
        return GeneAvgResult(
            value=mean_val, p_value=math.nan, std_error=math.nan,
            n_genes_used=1, n_genes_excluded=n_excluded,
        )
    sem = float(np.std(clean, ddof=1) / np.sqrt(n_used))
    if sem < 1e-15:
        p_val = 1.0 if abs(mean_val - null_value) < 1e-12 else 0.0
    else:
        _, p_val = ttest_1samp(clean, null_value)
        p_val = float(p_val)
    return GeneAvgResult(
        value=mean_val, p_value=p_val, std_error=sem,
        n_genes_used=n_used, n_genes_excluded=n_excluded,
    )


def gene_avg_enrichment(per_gene_values: list[float], n_total_genes: int) -> GeneAvgResult:
    return _aggregate(per_gene_values, n_total_genes, null_value=1.0)


def gene_avg_rate_ratio(per_gene_values: list[float], n_total_genes: int) -> GeneAvgResult:
    return _aggregate(per_gene_values, n_total_genes, null_value=1.0)


def gene_avg_auc(per_gene_values: list[float], n_total_genes: int) -> GeneAvgResult:
    return _aggregate(per_gene_values, n_total_genes, null_value=0.5)


def gene_avg_auprc(per_gene_values: list[float], n_total_genes: int) -> GeneAvgResult:
    return _aggregate(per_gene_values, n_total_genes, null_value=math.nan)
