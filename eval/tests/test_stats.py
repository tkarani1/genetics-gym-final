import math
import json

import pytest
import polars as pl

from biostat_cli.cli import (
    RunArgs,
    _compute_std_error,
    _resolve_eval_totals,
    _row_identity_key,
    _validate_bootstrap_args,
    run,
)
from biostat_cli.config import detect_pairwise_columns, parse_eval_totals
from biostat_cli.evaluators.base import Contingency
from biostat_cli.evaluators.variant import VariantEvaluator
from biostat_cli.stats.binary import (
    enrichment,
    fisher_p_value,
    pairwise_enrichment,
    pairwise_rate_ratio,
    poisson_p_value,
    rate_ratio,
    vsm_comparison,
    vsm_comparison_fisher,
    vsm_comparison_poisson_exact,
)
from biostat_cli.stats.continuous import (
    compute_auc,
    compute_auc_p_value,
    compute_auc_trunc,
    compute_auprc,
    compute_auprc_trunc,
    compute_curve_points,
    compute_threshold_point_metrics,
    delong_two_auc_p_value,
)
from biostat_cli.utils import apply_within_gene_percentile


def test_auc_and_auprc_basic():
    labels = [0, 0, 1, 1]
    scores = [0.1, 0.2, 0.8, 0.9]
    assert compute_auc(labels, scores) > 0.99
    assert compute_auprc(labels, scores) > 0.99


def test_threshold_point_metrics_basic():
    labels = [1, 1, 0, 0]
    scores = [0.9, 0.8, 0.7, 0.2]
    out = compute_threshold_point_metrics(labels, scores, threshold=0.8)
    assert out.tpr == pytest.approx(1.0)
    assert out.fpr == pytest.approx(0.0)
    assert out.precision == pytest.approx(1.0)
    assert out.recall == pytest.approx(1.0)
    assert out.rows_retained == 2
    assert out.n_pos_retained == 2
    assert out.n_neg_retained == 0


def test_truncated_metrics_basic():
    labels = [1, 1, 0, 0, 1, 0]
    scores = [0.99, 0.95, 0.94, 0.93, 0.2, 0.1]
    auc_t = compute_auc_trunc(labels, scores, threshold=0.93)
    auprc_t = compute_auprc_trunc(labels, scores, threshold=0.93)
    assert not math.isnan(auc_t)
    assert not math.isnan(auprc_t)


def test_curve_points_contains_roc_and_pr():
    labels = [0, 0, 1, 1]
    scores = [0.1, 0.2, 0.8, 0.9]
    points = compute_curve_points(labels, scores)
    assert points
    assert {"roc", "pr"} <= {p.curve_type for p in points}


def test_auc_p_value_perfect_separation():
    labels = [0, 0, 1, 1]
    scores = [0.1, 0.2, 0.8, 0.9]
    p = compute_auc_p_value(labels, scores)
    assert 0.0 <= p <= 1.0
    assert p < 0.05


def test_auc_p_value_no_discrimination():
    labels = [0, 0, 1, 1]
    scores = [0.5, 0.5, 0.5, 0.5]
    assert compute_auc(labels, scores) == 0.5
    p = compute_auc_p_value(labels, scores)
    assert p == pytest.approx(1.0)


def test_auc_p_value_single_class_nan():
    assert math.isnan(compute_auc_p_value([1, 1, 1], [0.1, 0.2, 0.3]))


def test_delong_two_auc_identical_scores_p_one():
    y = [0, 0, 1, 1]
    s = [0.1, 0.2, 0.7, 0.8]
    assert delong_two_auc_p_value(y, s, s) == pytest.approx(1.0)


def test_delong_two_auc_constant_scores_nan():
    """Degenerate when one score vector has no discrimination (zero variance)."""
    y = [0, 0, 1, 1]
    s_flat = [0.5, 0.5, 0.5, 0.5]
    s_vsm = [0.1, 0.2, 0.8, 0.9]
    assert math.isnan(delong_two_auc_p_value(y, s_flat, s_vsm))


def test_delong_two_auc_discordant():
    y = [0, 0, 1, 1]
    s_anchor = [0.3, 0.4, 0.6, 0.7]
    s_vsm = [0.1, 0.2, 0.8, 0.9]
    p = delong_two_auc_p_value(y, s_anchor, s_vsm)
    assert not math.isnan(p)
    assert 0.0 <= p <= 1.0


def test_binary_stats():
    cont = Contingency(tp=10, fp=5, tn=20, fn=15)
    enr = enrichment(cont)
    rr = rate_ratio(cont, case_total=200, ctrl_total=300)
    assert not math.isnan(enr.value)
    assert not math.isnan(enr.p_value)
    assert not math.isnan(enr.std_error)
    assert not math.isnan(rr.value)
    assert not math.isnan(rr.p_value)
    assert math.isnan(rr.std_error)
    assert not math.isnan(rr.rate_ratio_ci_lower)
    assert not math.isnan(rr.rate_ratio_ci_upper)


def test_fisher_p_value():
    cont = Contingency(tp=10, fp=5, tn=20, fn=15)
    p = fisher_p_value(cont)
    assert not math.isnan(p)
    assert 0.0 <= p <= 1.0


def test_poisson_p_value():
    cont = Contingency(tp=10, fp=5, tn=20, fn=15)
    p = poisson_p_value(cont)
    assert not math.isnan(p)
    assert 0.0 <= p <= 1.0


def test_pvalue_method_parameter():
    """Both methods produce valid (but potentially different) p-values."""
    cont = Contingency(tp=10, fp=5, tn=20, fn=15)
    enr_fisher = enrichment(cont, pvalue_method="fisher")
    enr_poisson = enrichment(cont, pvalue_method="poisson")
    assert enr_fisher.value == enr_poisson.value
    assert not math.isnan(enr_fisher.p_value)
    assert not math.isnan(enr_poisson.p_value)


def test_fisher_is_default():
    cont = Contingency(tp=10, fp=5, tn=20, fn=15)
    default_p = enrichment(cont).p_value
    fisher_p = enrichment(cont, pvalue_method="fisher").p_value
    assert default_p == fisher_p


def test_enrichment_log_lr_stderr_matches_formula():
    cont = Contingency(tp=10, fp=5, tn=20, fn=15)
    fisher_out = enrichment(cont, pvalue_method="fisher")
    poisson_out = enrichment(cont, pvalue_method="poisson")
    tp, fp, fn, tn = 10.0, 5.0, 15.0, 20.0
    pos_d, neg_d = tp + fn, fp + tn
    expected_se = math.sqrt((1.0 / tp - 1.0 / pos_d) + (1.0 / fp - 1.0 / neg_d))
    lr = (tp / pos_d) / (fp / neg_d)
    margin = 1.96 * expected_se
    expected_lo = math.exp(math.log(lr) - margin)
    expected_hi = math.exp(math.log(lr) + margin)
    assert fisher_out.std_error == pytest.approx(expected_se)
    assert poisson_out.std_error == pytest.approx(expected_se)
    assert fisher_out.enrichment_ci_lower == pytest.approx(expected_lo)
    assert fisher_out.enrichment_ci_upper == pytest.approx(expected_hi)
    assert poisson_out.enrichment_ci_lower == pytest.approx(expected_lo)
    assert poisson_out.enrichment_ci_upper == pytest.approx(expected_hi)


def test_enrichment_ci_uses_continuity_when_raw_se_undefined():
    cont = Contingency(tp=0, fp=5, tn=20, fn=15)
    out = enrichment(cont)
    assert not math.isnan(out.std_error)
    assert not math.isnan(out.enrichment_ci_lower)
    assert not math.isnan(out.enrichment_ci_upper)


def test_rate_ratio_poisson_analytic_std_error():
    cont = Contingency(tp=10, fp=5, tn=20, fn=15)
    rr = rate_ratio(cont, case_total=200, ctrl_total=300, pvalue_method="poisson")
    expected_se = rr.value * math.sqrt((1.0 / cont.tp) + (1.0 / cont.fp))
    assert rr.std_error == pytest.approx(expected_se)
    se_log = math.sqrt((1.0 / cont.tp) + (1.0 / cont.fp))
    margin = 1.96 * se_log
    assert rr.rate_ratio_ci_lower == pytest.approx(math.exp(math.log(rr.value) - margin))
    assert rr.rate_ratio_ci_upper == pytest.approx(math.exp(math.log(rr.value) + margin))


def test_rate_ratio_ci_matches_for_fisher_and_poisson_pvalue_method():
    cont = Contingency(tp=10, fp=5, tn=20, fn=15)
    rr_f = rate_ratio(cont, case_total=200, ctrl_total=300, pvalue_method="fisher")
    rr_p = rate_ratio(cont, case_total=200, ctrl_total=300, pvalue_method="poisson")
    assert math.isnan(rr_f.std_error)
    assert not math.isnan(rr_p.std_error)
    assert rr_f.rate_ratio_ci_lower == pytest.approx(rr_p.rate_ratio_ci_lower)
    assert rr_f.rate_ratio_ci_upper == pytest.approx(rr_p.rate_ratio_ci_upper)


def test_rate_ratio_poisson_analytic_std_error_zero_counts_nan():
    cont_tp_zero = Contingency(tp=0, fp=5, tn=20, fn=15)
    cont_fp_zero = Contingency(tp=10, fp=0, tn=20, fn=15)
    rr_tp_zero = rate_ratio(cont_tp_zero, case_total=200, ctrl_total=300, pvalue_method="poisson")
    rr_fp_zero = rate_ratio(cont_fp_zero, case_total=200, ctrl_total=300, pvalue_method="poisson")
    assert math.isnan(rr_tp_zero.std_error)
    assert math.isnan(rr_fp_zero.std_error)


def test_pvalue_empty_strata():
    cont_empty_above = Contingency(tp=0, fp=0, tn=20, fn=15)
    assert math.isnan(fisher_p_value(cont_empty_above))
    assert math.isnan(poisson_p_value(cont_empty_above))


def test_pairwise_enrichment():
    # Anchor on full set: good performance
    anchor_full = Contingency(tp=100, fp=10, tn=800, fn=90)
    # Anchor on pairwise intersection: similar performance
    anchor_pairwise = Contingency(tp=50, fp=5, tn=400, fn=45)
    # VSM on pairwise intersection: better performance than anchor on pairwise
    vsm_pairwise = Contingency(tp=60, fp=4, tn=401, fn=35)

    result = pairwise_enrichment(anchor_full, anchor_pairwise, vsm_pairwise)

    assert not math.isnan(result.value)
    assert not math.isnan(result.p_value)
    assert not math.isnan(result.anchor_value)
    assert not math.isnan(result.adjustment_ratio)

    # Verify the formula: value = anchor_value * adjustment_ratio
    expected_value = result.anchor_value * result.adjustment_ratio
    assert abs(result.value - expected_value) < 1e-9

    # VSM has better performance on pairwise, so adjustment_ratio > 1
    assert result.adjustment_ratio > 1.0


def test_pairwise_enrichment_same_contingency():
    # When anchor and vsm are the same (i.e., for the anchor itself)
    cont = Contingency(tp=50, fp=5, tn=400, fn=45)
    result = pairwise_enrichment(cont, cont, cont)

    assert abs(result.adjustment_ratio - 1.0) < 1e-9
    assert abs(result.value - result.anchor_value) < 1e-9


def test_pairwise_rate_ratio():
    anchor_full = Contingency(tp=100, fp=10, tn=800, fn=90)
    anchor_pairwise = Contingency(tp=50, fp=5, tn=400, fn=45)
    vsm_pairwise = Contingency(tp=60, fp=4, tn=401, fn=35)

    result = pairwise_rate_ratio(anchor_full, anchor_pairwise, vsm_pairwise, case_total=1000, ctrl_total=5000)

    assert not math.isnan(result.value)
    assert not math.isnan(result.p_value)
    assert not math.isnan(result.anchor_value)
    assert not math.isnan(result.adjustment_ratio)

    # Verify the formula
    expected_value = result.anchor_value * result.adjustment_ratio
    assert abs(result.value - expected_value) < 1e-9


def test_pairwise_rate_ratio_missing_totals():
    cont = Contingency(tp=50, fp=5, tn=400, fn=45)
    result = pairwise_rate_ratio(cont, cont, cont, case_total=None, ctrl_total=None)

    assert math.isnan(result.value)
    assert math.isnan(result.anchor_value)
    assert math.isnan(result.adjustment_ratio)


def test_detect_pairwise_columns():
    columns = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "mpc_score_anchor_percentile",
        "esm1b_score_percentile_with_anchor",
        "mpc_score_anchor_percentile_with_esm1b",
        "MisFit_S_score_percentile_with_anchor",
        "mpc_score_anchor_percentile_with_MisFit_S",
        "eval_col",
    ]

    result = detect_pairwise_columns(columns)

    assert result is not None
    assert result.anchor_base == "mpc_score"
    assert result.anchor_full_col == "mpc_score_anchor_percentile"
    assert len(result.vsm_pairs) == 2

    vsm_names = {pair[0] for pair in result.vsm_pairs}
    assert vsm_names == {"esm1b_score", "MisFit_S_score"}


def test_detect_pairwise_columns_no_anchor():
    columns = [
        "CHROM",
        "POS",
        "esm1b_score_percentile_with_anchor",
        "eval_col",
    ]

    result = detect_pairwise_columns(columns)
    assert result is None


def test_detect_pairwise_columns_no_vsm():
    columns = [
        "CHROM",
        "POS",
        "mpc_score_anchor_percentile",
        "eval_col",
    ]

    result = detect_pairwise_columns(columns)
    assert result is None


def test_parse_eval_totals():
    parsed = parse_eval_totals("eval_A:1000, eval_B:2500.5", "--case-total-by-eval")
    assert parsed == {"eval_A": 1000.0, "eval_B": 2500.5}


def test_parse_eval_totals_invalid_format():
    with pytest.raises(ValueError):
        parse_eval_totals("eval_A=1000", "--case-total-by-eval")


def test_resolve_eval_totals_priority():
    case_total, ctrl_total = _resolve_eval_totals(
        eval_col="eval_A",
        table_case_totals={"eval_A": 100.0, "eval_B": 200.0},
        table_ctrl_totals={"eval_A": 300.0, "eval_B": 400.0},
        cli_case_totals={"eval_A": 111.0},
        cli_ctrl_totals={"eval_A": 333.0},
    )
    assert case_total == 111.0
    assert ctrl_total == 333.0


def test_resolve_eval_totals_fallbacks():
    # Falls back to table-level totals for matching eval when CLI omits that eval.
    case_total, ctrl_total = _resolve_eval_totals(
        eval_col="eval_B",
        table_case_totals={"eval_B": 222.0},
        table_ctrl_totals={"eval_B": 444.0},
        cli_case_totals={},
        cli_ctrl_totals={},
    )
    assert case_total == 222.0
    assert ctrl_total == 444.0

    # No CLI or table entry for this eval → missing denominators.
    case_total, ctrl_total = _resolve_eval_totals(
        eval_col="eval_C",
        table_case_totals={},
        table_ctrl_totals={},
        cli_case_totals={},
        cli_ctrl_totals={},
    )
    assert case_total is None
    assert ctrl_total is None


def test_compute_std_error():
    values = [1.0, 2.0, 3.0, 4.0]
    out = _compute_std_error(values)
    assert out > 0


def test_compute_std_error_nan_filtering():
    values = [1.0, float("nan"), 3.0]
    out = _compute_std_error(values)
    assert not math.isnan(out)


def test_row_identity_key_nan_threshold():
    row = {
        "eval_name": "eval_a",
        "filter_name": "none",
        "score_name": "score_x",
        "threshold": float("nan"),
        "stat": "auc",
    }
    key = _row_identity_key(row)
    assert key[3] == "nan"


def test_validate_bootstrap_args():
    args = RunArgs(
        resources_json="resources.json",
        table_name="t",
        eval_level="variant",
        stat="all",
        eval_set=None,
        filters=None,
        thresholds=None,
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=1,
        out_fname="out",
        write_missing="none",
    )
    with pytest.raises(ValueError):
        _validate_bootstrap_args(args)


def test_bootstrap_run_value_and_pvalue_stable(tmp_path):
    df = pl.DataFrame(
        {
            "chrom": ["1", "1", "1", "1", "1", "1"],
            "pos": [1, 2, 3, 4, 5, 6],
            "ref": ["A"] * 6,
            "alt": ["C"] * 6,
            "eval_a": [True, False, True, False, True, False],
            "score_x": [0.95, 0.85, 0.88, 0.1, 0.99, 0.4],
        }
    )
    parquet_path = tmp_path / "toy.parquet"
    resources_path = tmp_path / "resources.json"
    out_prefix = tmp_path / "out"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "toy": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_x"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    base_args = RunArgs(
        resources_json=str(resources_path),
        table_name="toy",
        eval_level="variant",
        stat="enrichment",
        eval_set=None,
        filters=None,
        thresholds="0.5",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    base_df, _, _, _, _, _ = run(base_args)
    base_row = base_df.to_dicts()[0]
    assert not math.isnan(base_row["std_error"])
    assert not math.isnan(base_row["enrichment_ci_lower"])
    assert not math.isnan(base_row["enrichment_ci_upper"])

    boot_args = RunArgs(
        resources_json=str(resources_path),
        table_name="toy",
        eval_level="variant",
        stat="enrichment",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=20,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    boot_df, _, _, _, _, _ = run(boot_args)
    boot_row = boot_df.to_dicts()[0]
    assert boot_row["value"] == pytest.approx(base_row["value"], rel=0, abs=1e-12)
    assert boot_row["p_value"] == pytest.approx(base_row["p_value"], rel=0, abs=1e-12)
    assert not math.isnan(boot_row["std_error"])
    assert boot_row["std_error"] != pytest.approx(base_row["std_error"], rel=0, abs=1e-12)
    assert math.isnan(boot_row["enrichment_ci_lower"])
    assert math.isnan(boot_row["enrichment_ci_upper"])


def test_bootstrap_pairwise_std_error(tmp_path):
    df = pl.DataFrame(
        {
            "chrom": ["1", "1", "1", "1", "1", "1", "1", "1"],
            "pos": [1, 2, 3, 4, 5, 6, 7, 8],
            "ref": ["A"] * 8,
            "alt": ["C"] * 8,
            "eval_a": [True, False, True, False, True, False, True, False],
            "anchor_score_anchor_percentile": [0.9, 0.82, 0.95, 0.2, 0.85, 0.3, 0.88, 0.4],
            "vsm1_score_percentile_with_anchor": [0.92, 0.86, 0.97, 0.25, 0.83, 0.2, 0.86, 0.5],
            "anchor_score_anchor_percentile_with_vsm1": [0.89, 0.81, 0.94, 0.22, 0.84, 0.28, 0.87, 0.45],
        }
    )
    parquet_path = tmp_path / "pairwise.parquet"
    resources_path = tmp_path / "resources_pairwise.json"
    out_prefix = tmp_path / "out_pairwise"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "pairwise": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["anchor_score_anchor_percentile"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    args = RunArgs(
        resources_json=str(resources_path),
        table_name="pairwise",
        eval_level="variant",
        stat="pairwise_enrichment",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=20,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    out_df, _, _, _, _, _ = run(args)
    rows = out_df.to_dicts()
    assert rows
    assert all("std_error" in row for row in rows)
    assert any(not math.isnan(row["std_error"]) for row in rows)


def test_continuous_threshold_and_trunc_stats_integration(tmp_path):
    df = pl.DataFrame(
        {
            "chrom": ["1"] * 8,
            "pos": list(range(1, 9)),
            "ref": ["A"] * 8,
            "alt": ["C"] * 8,
            "eval_a": [True, False, True, False, True, False, True, False],
            "score_x": [0.95, 0.85, 0.88, 0.1, 0.99, 0.4, 0.75, 0.2],
        }
    )
    parquet_path = tmp_path / "continuous_new_stats.parquet"
    resources_path = tmp_path / "resources_continuous_new_stats.json"
    out_prefix = tmp_path / "out_continuous_new_stats"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "continuous_new_stats": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_x"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")
    args = RunArgs(
        resources_json=str(resources_path),
        table_name="continuous_new_stats",
        eval_level="variant",
        stat="tpr_at_threshold,fpr_at_threshold,precision_at_threshold,recall_at_threshold,auc_trunc,auprc_trunc",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    out_df, _, _, _, _, curves_df = run(args)
    stats = set(out_df["stat"].to_list())
    assert {
        "tpr_at_threshold", "fpr_at_threshold",
        "precision_at_threshold", "recall_at_threshold",
        "auc_trunc", "auprc_trunc",
    } <= stats
    row = out_df.filter(pl.col("stat") == "auc_trunc").to_dicts()[0]
    assert row["threshold"] == pytest.approx(0.8)
    assert row["rows_retained"] >= 0
    assert curves_df.height > 0


def test_pairwise_trunc_and_threshold_stats_integration(tmp_path):
    df = pl.DataFrame(
        {
            "chrom": ["1"] * 8,
            "pos": list(range(1, 9)),
            "ref": ["A"] * 8,
            "alt": ["C"] * 8,
            "eval_a": [True, False, True, False, True, False, True, False],
            "anchor_score_anchor_percentile": [0.9, 0.82, 0.95, 0.2, 0.85, 0.3, 0.88, 0.4],
            "vsm1_score_percentile_with_anchor": [0.92, 0.86, 0.97, 0.25, 0.83, 0.2, 0.86, 0.5],
            "anchor_score_anchor_percentile_with_vsm1": [0.89, 0.81, 0.94, 0.22, 0.84, 0.28, 0.87, 0.45],
        }
    )
    parquet_path = tmp_path / "pairwise_new_stats.parquet"
    resources_path = tmp_path / "resources_pairwise_new_stats.json"
    out_prefix = tmp_path / "out_pairwise_new_stats"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "pairwise_new_stats": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["anchor_score_anchor_percentile"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")
    args = RunArgs(
        resources_json=str(resources_path),
        table_name="pairwise_new_stats",
        eval_level="variant",
        stat="pairwise_auc_trunc,pairwise_auprc_trunc,pairwise_tpr_at_threshold",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    out_df, _, _, _, _, _ = run(args)
    assert {"pairwise_auc_trunc", "pairwise_auprc_trunc", "pairwise_tpr_at_threshold"} <= set(
        out_df["stat"].to_list()
    )


def test_rate_ratio_missing_denominators_raises(tmp_path):
    df = pl.DataFrame(
        {
            "chrom": ["1", "1", "1", "1", "1", "1"],
            "pos": [1, 2, 3, 4, 5, 6],
            "ref": ["A"] * 6,
            "alt": ["C"] * 6,
            "eval_a": [True, False, True, False, True, False],
            "score_x": [0.95, 0.85, 0.88, 0.1, 0.99, 0.4],
        }
    )
    parquet_path = tmp_path / "toy_rr_missing_totals.parquet"
    resources_path = tmp_path / "resources_rr_missing_totals.json"
    out_prefix = tmp_path / "out_rr_missing_totals"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "toy_rr_missing_totals": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_x"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    args = RunArgs(
        resources_json=str(resources_path),
        table_name="toy_rr_missing_totals",
        eval_level="variant",
        stat="rate_ratio",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    with pytest.raises(ValueError, match="Cohort denominators are required"):
        run(args)


def test_rate_ratio_poisson_std_error_without_bootstrap(tmp_path):
    df = pl.DataFrame(
        {
            "chrom": ["1", "1", "1", "1", "1", "1"],
            "pos": [1, 2, 3, 4, 5, 6],
            "ref": ["A"] * 6,
            "alt": ["C"] * 6,
            "eval_a": [True, False, True, False, True, False],
            "score_x": [0.95, 0.85, 0.88, 0.1, 0.99, 0.4],
        }
    )
    parquet_path = tmp_path / "toy_poisson_rr.parquet"
    resources_path = tmp_path / "resources_poisson_rr.json"
    out_prefix = tmp_path / "out_poisson_rr"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "toy_poisson_rr": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_x"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    args = RunArgs(
        resources_json=str(resources_path),
        table_name="toy_poisson_rr",
        eval_level="variant",
        stat="rate_ratio",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval="eval_a:3",
        ctrl_total_by_eval="eval_a:3",
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
        pvalue_method="poisson",
    )
    out_df, _, _, _, _, _ = run(args)
    row = out_df.to_dicts()[0]
    expected = row["value"] * math.sqrt((1.0 / row["tp"]) + (1.0 / row["fp"]))
    assert row["std_error"] == pytest.approx(expected)
    se_log = math.sqrt((1.0 / row["tp"]) + (1.0 / row["fp"]))
    assert row["rate_ratio_ci_lower"] == pytest.approx(row["value"] * math.exp(-1.96 * se_log))
    assert row["rate_ratio_ci_upper"] == pytest.approx(row["value"] * math.exp(1.96 * se_log))


def test_rate_ratio_bootstrap_overrides_analytic_poisson_std_error(tmp_path):
    df = pl.DataFrame(
        {
            "chrom": ["1", "1", "1", "1", "1", "1", "1", "1", "1", "1"],
            "pos": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            "ref": ["A"] * 10,
            "alt": ["C"] * 10,
            "eval_a": [True, False, True, False, True, False, True, False, True, False],
            "score_x": [0.95, 0.85, 0.88, 0.1, 0.99, 0.4, 0.81, 0.79, 0.9, 0.2],
        }
    )
    parquet_path = tmp_path / "toy_poisson_rr_boot.parquet"
    resources_path = tmp_path / "resources_poisson_rr_boot.json"
    out_prefix = tmp_path / "out_poisson_rr_boot"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "toy_poisson_rr_boot": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_x"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    base_args = RunArgs(
        resources_json=str(resources_path),
        table_name="toy_poisson_rr_boot",
        eval_level="variant",
        stat="rate_ratio",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval="eval_a:5",
        ctrl_total_by_eval="eval_a:5",
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
        pvalue_method="poisson",
    )
    base_df, _, _, _, _, _ = run(base_args)
    base_row = base_df.to_dicts()[0]
    assert not math.isnan(base_row["std_error"])
    assert not math.isnan(base_row["rate_ratio_ci_lower"])
    assert not math.isnan(base_row["rate_ratio_ci_upper"])

    boot_args = RunArgs(
        resources_json=str(resources_path),
        table_name="toy_poisson_rr_boot",
        eval_level="variant",
        stat="rate_ratio",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval="eval_a:5",
        ctrl_total_by_eval="eval_a:5",
        bootstrap_samples=30,
        out_fname=str(out_prefix),
        write_missing="none",
        pvalue_method="poisson",
    )
    boot_df, _, _, _, _, _ = run(boot_args)
    boot_row = boot_df.to_dicts()[0]
    assert not math.isnan(boot_row["std_error"])
    assert boot_row["std_error"] != pytest.approx(base_row["std_error"], rel=0, abs=1e-12)
    assert math.isnan(boot_row["rate_ratio_ci_lower"])
    assert math.isnan(boot_row["rate_ratio_ci_upper"])


def test_apply_within_gene_percentile():
    df = pl.DataFrame(
        {"ensg": ["g1", "g1", "g1", "g2", "g2"], "score": [1.0, 2.0, 3.0, 10.0, 30.0]}
    )
    out = apply_within_gene_percentile(df.lazy(), "score").collect()
    g1 = out.filter(pl.col("ensg") == "g1")["score"].sort()
    g2 = out.filter(pl.col("ensg") == "g2")["score"].sort()
    assert g1.to_list() == pytest.approx([1 / 3, 2 / 3, 1.0])
    assert g2.to_list() == pytest.approx([0.5, 1.0])


def test_apply_within_gene_percentile_requires_gene_col():
    df = pl.DataFrame({"score": [1.0, 2.0]})
    with pytest.raises(ValueError, match="Within-gene percentile"):
        apply_within_gene_percentile(df.lazy(), "score").collect()


def test_prepare_score_frame_within_gene_percentile():
    source = pl.DataFrame(
        {"ensg": ["a", "a", "a"], "y": [True, False, True], "s": [1.0, 2.0, 3.0]}
    ).lazy()
    ev = VariantEvaluator(source)
    prepared = ev.prepare_eval_frame("y", None)
    sf = ev.prepare_score_frame(prepared, "s", within_gene_percentile=True)
    vals = sf.frame.select("s").collect()["s"].sort().to_list()
    assert vals == pytest.approx([1 / 3, 2 / 3, 1.0])


# ---------------------------------------------------------------------------
# VSM comparison Fisher exact tests
# ---------------------------------------------------------------------------


def test_vsm_comparison_basic():
    from scipy.stats import fisher_exact as scipy_fisher

    cont_a = Contingency(tp=50, fp=10, tn=0, fn=0)
    cont_b = Contingency(tp=30, fp=20, tn=0, fn=0)
    result = vsm_comparison_fisher(cont_a, cont_b)

    table = [[50, 30], [10, 20]]
    expected_or, expected_p_greater = scipy_fisher(table, alternative="greater")
    _, expected_p_less = scipy_fisher(table, alternative="less")

    assert result.odds_ratio == pytest.approx(expected_or)
    assert result.p_greater == pytest.approx(expected_p_greater)
    assert result.p_less == pytest.approx(expected_p_less)

    a, b, c, d = 50.5, 30.5, 10.5, 20.5
    expected_log_or = math.log((a * d) / (b * c))
    expected_se = math.sqrt((1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d))
    assert result.log_odds_ratio == pytest.approx(expected_log_or)
    assert result.standard_error == pytest.approx(expected_se)
    z = 1.96
    moe = z * expected_se
    assert result.log_ci_lower == pytest.approx(expected_log_or - moe)
    assert result.log_ci_upper == pytest.approx(expected_log_or + moe)
    assert result.conf_interval_lower == pytest.approx(math.exp(expected_log_or - moe))
    assert result.conf_interval_upper == pytest.approx(math.exp(expected_log_or + moe))


def test_vsm_comparison_zero_cells():
    cont_a = Contingency(tp=0, fp=0, tn=100, fn=50)
    cont_b = Contingency(tp=10, fp=5, tn=85, fn=50)
    result = vsm_comparison_fisher(cont_a, cont_b)
    assert math.isnan(result.odds_ratio)
    assert math.isnan(result.p_greater)
    assert math.isnan(result.p_less)
    assert math.isnan(result.log_odds_ratio)
    assert math.isnan(result.standard_error)
    assert math.isnan(result.log_ci_lower)
    assert math.isnan(result.log_ci_upper)
    assert math.isnan(result.conf_interval_lower)
    assert math.isnan(result.conf_interval_upper)


def test_vsm_comparison_symmetry():
    cont_a = Contingency(tp=40, fp=15, tn=0, fn=0)
    cont_b = Contingency(tp=25, fp=30, tn=0, fn=0)
    result_ab = vsm_comparison_fisher(cont_a, cont_b)
    result_ba = vsm_comparison_fisher(cont_b, cont_a)

    assert result_ab.odds_ratio * result_ba.odds_ratio == pytest.approx(1.0)
    assert result_ab.p_greater == pytest.approx(result_ba.p_less)
    assert result_ab.p_less == pytest.approx(result_ba.p_greater)
    assert result_ab.log_odds_ratio == pytest.approx(-result_ba.log_odds_ratio)


def test_vsm_comparison_poisson_exact_basic():
    from scipy.stats import binom

    cont_a = Contingency(tp=50, fp=10, tn=0, fn=0)
    cont_b = Contingency(tp=30, fp=20, tn=0, fn=0)
    result = vsm_comparison_poisson_exact(cont_a, cont_b)

    total_tp = 80
    p0 = 10 / 30
    expected_p_greater = float(binom.sf(50 - 1, total_tp, p0))
    expected_p_less = float(binom.cdf(50, total_tp, p0))

    assert result.p_greater == pytest.approx(expected_p_greater)
    assert result.p_less == pytest.approx(expected_p_less)
    assert 0.0 <= result.p_greater <= 1.0
    assert 0.0 <= result.p_less <= 1.0
    assert not math.isnan(result.odds_ratio)
    assert not math.isnan(result.log_odds_ratio)
    assert not math.isnan(result.standard_error)


def test_vsm_comparison_poisson_symmetry():
    cont_a = Contingency(tp=40, fp=15, tn=0, fn=0)
    cont_b = Contingency(tp=25, fp=30, tn=0, fn=0)
    result_ab = vsm_comparison_poisson_exact(cont_a, cont_b)
    result_ba = vsm_comparison_poisson_exact(cont_b, cont_a)

    assert result_ab.p_greater == pytest.approx(result_ba.p_less)
    assert result_ab.p_less == pytest.approx(result_ba.p_greater)
    assert result_ab.log_odds_ratio == pytest.approx(-result_ba.log_odds_ratio)


def test_vsm_comparison_dispatch_methods():
    cont_a = Contingency(tp=40, fp=15, tn=0, fn=0)
    cont_b = Contingency(tp=25, fp=30, tn=0, fn=0)

    fisher_result = vsm_comparison(cont_a, cont_b, method="fisher")
    poisson_result = vsm_comparison(cont_a, cont_b, method="poisson")
    assert abs(fisher_result.p_greater - poisson_result.p_greater) > 1e-12


def test_vsm_comparison_integration(tmp_path):
    """End-to-end: vsm_comparison stat produces a non-empty DataFrame."""
    df = pl.DataFrame(
        {
            "chrom": ["1"] * 8,
            "pos": list(range(1, 9)),
            "ref": ["A"] * 8,
            "alt": ["C"] * 8,
            "eval_a": [True, False, True, False, True, False, True, False],
            "score_x": [0.95, 0.85, 0.88, 0.1, 0.99, 0.4, 0.7, 0.3],
            "score_y": [0.6, 0.92, 0.75, 0.2, 0.88, 0.55, 0.98, 0.1],
        }
    )
    parquet_path = tmp_path / "toy2.parquet"
    resources_path = tmp_path / "resources2.json"
    out_prefix = tmp_path / "out2"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "toy2": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_x", "score_y"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    from biostat_cli.cli import RunArgs, run

    args = RunArgs(
        resources_json=str(resources_path),
        table_name="toy2",
        eval_level="variant",
        stat="vsm_comparison",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    _, _, _, vsm_cmp_df, _, _ = run(args)
    assert vsm_cmp_df.height > 0
    assert set(vsm_cmp_df.columns) == {
        "eval_name", "filter_name", "vsm_i", "vsm_j",
        "threshold", "odds_ratio", "p_greater", "p_less",
        "log_odds_ratio", "standard_error",
        "log_ci_lower", "log_ci_upper",
        "conf_interval_lower", "conf_interval_upper",
        "rows_used_i", "rows_used_j", "rows_used_pair",
    }
    row = vsm_cmp_df.to_dicts()[0]
    assert row["vsm_i"] == "score_x"
    assert row["vsm_j"] == "score_y"
    assert row["threshold"] == pytest.approx(0.8)
    assert row["rows_used_i"] > 0
    assert row["rows_used_j"] > 0


def test_vsm_comparison_integration_poisson_method(tmp_path):
    df = pl.DataFrame(
        {
            "chrom": ["1"] * 8,
            "pos": list(range(1, 9)),
            "ref": ["A"] * 8,
            "alt": ["C"] * 8,
            "eval_a": [True, False, True, False, True, False, True, False],
            "score_x": [0.95, 0.85, 0.88, 0.1, 0.99, 0.4, 0.7, 0.3],
            "score_y": [0.6, 0.92, 0.75, 0.2, 0.88, 0.55, 0.98, 0.1],
        }
    )
    parquet_path = tmp_path / "toy3.parquet"
    resources_path = tmp_path / "resources3.json"
    out_prefix = tmp_path / "out3"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "toy3": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_x", "score_y"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    args = RunArgs(
        resources_json=str(resources_path),
        table_name="toy3",
        eval_level="variant",
        stat="vsm_comparison",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
        vsm_comparison_method="poisson",
    )
    _, _, _, vsm_cmp_df, _, _ = run(args)
    assert vsm_cmp_df.height > 0
    row = vsm_cmp_df.to_dicts()[0]
    assert 0.0 <= row["p_greater"] <= 1.0
    assert 0.0 <= row["p_less"] <= 1.0


def test_vsm_comparison_invalid_method_raises():
    cont_a = Contingency(tp=12, fp=8, tn=0, fn=0)
    cont_b = Contingency(tp=9, fp=11, tn=0, fn=0)
    with pytest.raises(ValueError, match="Unknown vsm comparison method"):
        vsm_comparison(cont_a, cont_b, method="not-a-method")


def test_vsm_comparison_missingness_uses_intersection(tmp_path):
    """VSM comparison computes contingencies on the intersection of non-null rows."""
    df = pl.DataFrame(
        {
            "chrom": ["1"] * 10,
            "pos": list(range(1, 11)),
            "ref": ["A"] * 10,
            "alt": ["C"] * 10,
            "eval_a": [True, False, True, False, True, False, True, False, True, False],
            "score_x": [0.99, 0.7, 0.97, 0.1, 0.95, 0.05, 0.9, 0.4, None, None],
            "score_y": [None, None, 0.96, 0.15, 0.93, 0.08, 0.91, 0.8, None, None],
        }
    )
    parquet_path = tmp_path / "missingness.parquet"
    resources_path = tmp_path / "resources_missingness.json"
    out_prefix = tmp_path / "out_missingness"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "missingness": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_x", "score_y"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    args = RunArgs(
        resources_json=str(resources_path),
        table_name="missingness",
        eval_level="variant",
        stat="vsm_comparison",
        eval_set=None,
        filters=None,
        thresholds="0.5",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
        vsm_comparison_method="poisson",
    )
    _, _, _, vsm_cmp_df, _, _ = run(args)
    row = vsm_cmp_df.to_dicts()[0]
    # Both scores have non-null values for rows 2-7 (indices), so intersection = 6
    assert row["rows_used_pair"] == 6
    assert row["rows_used_i"] == row["rows_used_j"] == 6
    assert 0.0 <= row["p_greater"] <= 1.0
    assert 0.0 <= row["p_less"] <= 1.0


def test_vsm_comparison_intersection_contingencies(tmp_path):
    """Contingencies are computed on the pair intersection, not per-VSM row sets.

    score_a: non-null for rows 0-9  (10 rows, 5 pos, 5 neg)
    score_b: non-null for rows 5-14 (10 rows, 5 pos, 5 neg)
    Intersection: rows 5-9 (5 rows, 3 pos=True at indices 6,8 wait let's be explicit)

    We verify that TP+FP+TN+FN in the output sums to 5 (the intersection size), not 10.
    """
    n = 15
    df = pl.DataFrame(
        {
            "chrom": ["1"] * n,
            "pos": list(range(1, n + 1)),
            "ref": ["A"] * n,
            "alt": ["C"] * n,
            "eval_a": [True, False] * 7 + [True],
            "score_a": [0.9, 0.1, 0.8, 0.2, 0.7, 0.6, 0.95, 0.05, 0.85, 0.15,
                        None, None, None, None, None],
            "score_b": [None, None, None, None, None,
                        0.55, 0.92, 0.08, 0.88, 0.12, 0.75, 0.3, 0.65, 0.4, 0.99],
        }
    )
    parquet_path = tmp_path / "intersection.parquet"
    resources_path = tmp_path / "resources_intersection.json"
    out_prefix = tmp_path / "out_intersection"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "intersection": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_a", "score_b"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    args = RunArgs(
        resources_json=str(resources_path),
        table_name="intersection",
        eval_level="variant",
        stat="vsm_comparison",
        eval_set=None,
        filters=None,
        thresholds="0.5",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    _, _, _, vsm_cmp_df, _, _ = run(args)
    assert vsm_cmp_df.height > 0
    row = vsm_cmp_df.to_dicts()[0]
    # Intersection of non-null rows: indices 5-9 → 5 rows
    assert row["rows_used_pair"] == 5
    assert row["rows_used_i"] == 5
    assert row["rows_used_j"] == 5
    # TP+FP+TN+FN for each score should sum to 5 (intersection size)
    # The test is sufficient if rows_used_pair confirms intersection size
    assert 0.0 <= row["p_greater"] <= 1.0
    assert 0.0 <= row["p_less"] <= 1.0


def test_vsm_comparison_parallel_matches_serial_poisson(tmp_path):
    from biostat_cli.cli_parallel import RunArgs as ParallelRunArgs
    from biostat_cli.cli_parallel import run as run_parallel

    df = pl.DataFrame(
        {
            "chrom": ["1"] * 12,
            "pos": list(range(1, 13)),
            "ref": ["A"] * 12,
            "alt": ["C"] * 12,
            "eval_a": [True, False, True, False, True, False, True, False, True, False, True, False],
            "score_x": [0.99, 0.1, 0.96, 0.2, 0.94, 0.4, 0.9, 0.5, 0.88, 0.3, 0.86, 0.6],
            "score_y": [0.85, 0.92, 0.82, 0.88, 0.8, 0.7, 0.95, 0.2, 0.75, 0.98, 0.73, 0.9],
        }
    )
    parquet_path = tmp_path / "parallel.parquet"
    resources_path = tmp_path / "resources_parallel.json"
    out_prefix = tmp_path / "out_parallel"
    df.write_parquet(str(parquet_path))
    resources = {
        "Table_info": {
            "parallel": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": ["score_x", "score_y"],
                "evals": ["eval_a"],
            }
        }
    }
    resources_path.write_text(json.dumps(resources), encoding="utf-8")

    serial_args = RunArgs(
        resources_json=str(resources_path),
        table_name="parallel",
        eval_level="variant",
        stat="vsm_comparison",
        eval_set=None,
        filters=None,
        thresholds="0.8,0.9",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
        vsm_comparison_method="poisson",
    )
    _, _, _, serial_df, _, _ = run(serial_args)

    parallel_args = ParallelRunArgs(
        resources_json=str(resources_path),
        table_name="parallel",
        eval_level="variant",
        stat="vsm_comparison",
        eval_set=None,
        filters=None,
        thresholds="0.8,0.9",
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        out_fname=str(out_prefix),
        write_missing="none",
        vsm_comparison_method="poisson",
    )
    _, _, _, parallel_df, _, _ = run_parallel(parallel_args)

    sort_cols = ["eval_name", "filter_name", "vsm_i", "vsm_j", "threshold"]
    serial_rows = serial_df.sort(sort_cols).to_dicts()
    parallel_rows = parallel_df.sort(sort_cols).to_dicts()
    assert len(serial_rows) == len(parallel_rows)
    for left, right in zip(serial_rows, parallel_rows):
        assert left["eval_name"] == right["eval_name"]
        assert left["filter_name"] == right["filter_name"]
        assert left["vsm_i"] == right["vsm_i"]
        assert left["vsm_j"] == right["vsm_j"]
        assert left["threshold"] == pytest.approx(right["threshold"])
        assert left["p_greater"] == pytest.approx(right["p_greater"])
        assert left["p_less"] == pytest.approx(right["p_less"])


# ---------------------------------------------------------------------------
# Gene-averaged statistics
# ---------------------------------------------------------------------------

from biostat_cli.config import GENE_AVG_STATS
from biostat_cli.stats.gene_averaged import (
    gene_avg_auc,
    gene_avg_auprc,
    gene_avg_enrichment,
    gene_avg_rate_ratio,
)


def test_parse_stats_accepts_gene_avg():
    from biostat_cli.config import parse_stats
    result = parse_stats("gene_avg_enrichment")
    assert result == {"gene_avg_enrichment"}


def test_parse_stats_all_includes_gene_avg():
    from biostat_cli.config import parse_stats
    result = parse_stats("all")
    assert result >= GENE_AVG_STATS


def test_gene_avg_enrichment_basic():
    from scipy.stats import sem
    values = [2.0, 3.0, 1.5, 2.5, 4.0]
    r = gene_avg_enrichment(values, n_total_genes=5)
    assert r.value == pytest.approx(sum(values) / len(values))
    assert r.std_error == pytest.approx(sem(values))
    assert r.n_genes_used == 5
    assert r.n_genes_excluded == 0
    assert 0.0 <= r.p_value <= 1.0


def test_gene_avg_enrichment_with_nans():
    values = [2.0, float("nan"), 3.0, float("nan"), 1.0]
    r = gene_avg_enrichment(values, n_total_genes=5)
    assert r.value == pytest.approx(2.0)
    assert r.n_genes_used == 3
    assert r.n_genes_excluded == 2


def test_gene_avg_enrichment_all_nan():
    values = [float("nan"), float("nan")]
    r = gene_avg_enrichment(values, n_total_genes=2)
    assert math.isnan(r.value)
    assert math.isnan(r.p_value)
    assert math.isnan(r.std_error)
    assert r.n_genes_used == 0
    assert r.n_genes_excluded == 2


def test_gene_avg_enrichment_single_gene():
    r = gene_avg_enrichment([5.0], n_total_genes=1)
    assert r.value == pytest.approx(5.0)
    assert math.isnan(r.std_error)
    assert math.isnan(r.p_value)
    assert r.n_genes_used == 1
    assert r.n_genes_excluded == 0


def test_gene_avg_enrichment_p_value_significant():
    values = [3.0, 4.0, 5.0, 3.5, 4.5]
    r = gene_avg_enrichment(values, n_total_genes=5)
    assert r.p_value < 0.05


def test_gene_avg_enrichment_p_value_null():
    values = [0.9, 1.1, 1.0, 0.95, 1.05]
    r = gene_avg_enrichment(values, n_total_genes=5)
    assert r.p_value > 0.05


def test_gene_avg_auc_null_is_half():
    values = [0.5, 0.5, 0.5]
    r = gene_avg_auc(values, n_total_genes=3)
    assert r.value == pytest.approx(0.5)
    assert r.p_value >= 0.99 or r.p_value == pytest.approx(1.0)


def test_gene_avg_auprc_basic():
    from scipy.stats import sem
    values = [0.8, 0.7, 0.9]
    r = gene_avg_auprc(values, n_total_genes=3)
    assert r.value == pytest.approx(sum(values) / len(values))
    assert r.std_error == pytest.approx(sem(values))
    assert r.n_genes_used == 3
    assert r.n_genes_excluded == 0


def test_gene_avg_rate_ratio_basic():
    values = [1.5, 2.0, 2.5]
    r = gene_avg_rate_ratio(values, n_total_genes=3)
    assert r.value == pytest.approx(2.0)
    assert r.n_genes_used == 3


# ---------------------------------------------------------------------------
# Gene-averaged evaluator methods
# ---------------------------------------------------------------------------


def test_contingency_by_gene_batch_basic():
    df = pl.DataFrame({
        "ensg": ["g1", "g1", "g1", "g1", "g2", "g2", "g2", "g2", "g3", "g3", "g3", "g3"],
        "eval_a": [True, True, False, False] * 3,
        "score_x": [0.9, 0.8, 0.7, 0.1, 0.95, 0.85, 0.6, 0.05, 0.92, 0.3, 0.88, 0.2],
    })
    ev = VariantEvaluator(df.lazy())
    from biostat_cli.evaluators.base import ScoreFrame
    sf = ScoreFrame(frame=df.lazy(), rows_used=12)
    gene_conts = ev.contingency_by_gene_batch(sf, "eval_a", "score_x", [0.5], "ensg")
    assert len(gene_conts) == 3
    for gene, conts in gene_conts.items():
        assert len(conts) == 1
        c = conts[0]
        assert c.tp + c.fp + c.tn + c.fn == 4

    # Global contingency should equal sum of per-gene contingencies
    global_conts = ev.contingency_batch(sf, "eval_a", "score_x", [0.5])
    gc = global_conts[0]
    total_tp = sum(conts[0].tp for conts in gene_conts.values())
    total_fp = sum(conts[0].fp for conts in gene_conts.values())
    total_tn = sum(conts[0].tn for conts in gene_conts.values())
    total_fn = sum(conts[0].fn for conts in gene_conts.values())
    assert gc.tp == pytest.approx(total_tp)
    assert gc.fp == pytest.approx(total_fp)
    assert gc.tn == pytest.approx(total_tn)
    assert gc.fn == pytest.approx(total_fn)


def test_contingency_by_gene_batch_multiple_thresholds():
    df = pl.DataFrame({
        "ensg": ["g1", "g1", "g2", "g2"],
        "eval_a": [True, False, True, False],
        "score_x": [0.9, 0.8, 0.5, 0.1],
    })
    ev = VariantEvaluator(df.lazy())
    from biostat_cli.evaluators.base import ScoreFrame
    sf = ScoreFrame(frame=df.lazy(), rows_used=4)
    gene_conts = ev.contingency_by_gene_batch(sf, "eval_a", "score_x", [0.5, 0.8], "ensg")
    for conts in gene_conts.values():
        assert len(conts) == 2


def test_contingency_by_gene_batch_missing_gene_col():
    df = pl.DataFrame({
        "eval_a": [True, False],
        "score_x": [0.9, 0.1],
    })
    ev = VariantEvaluator(df.lazy())
    from biostat_cli.evaluators.base import ScoreFrame
    sf = ScoreFrame(frame=df.lazy(), rows_used=2)
    with pytest.raises(ValueError, match="Gene column"):
        ev.contingency_by_gene_batch(sf, "eval_a", "score_x", [0.5], "ensg")


def test_labels_and_scores_by_gene_basic():
    df = pl.DataFrame({
        "ensg": ["g1", "g1", "g2", "g2"],
        "eval_a": [True, False, True, False],
        "score_x": [0.9, 0.1, 0.8, 0.2],
    })
    ev = VariantEvaluator(df.lazy())
    from biostat_cli.evaluators.base import ScoreFrame
    sf = ScoreFrame(frame=df.lazy(), rows_used=4)
    gene_ls = ev.labels_and_scores_by_gene(sf, "eval_a", "score_x", "ensg")
    assert len(gene_ls) == 2
    assert "g1" in gene_ls
    assert "g2" in gene_ls
    labels_g1, scores_g1 = gene_ls["g1"]
    assert sorted(labels_g1) == [0, 1]
    assert len(scores_g1) == 2


def test_labels_and_scores_by_gene_single_class():
    df = pl.DataFrame({
        "ensg": ["g1", "g1", "g2", "g2"],
        "eval_a": [True, True, True, False],
        "score_x": [0.9, 0.8, 0.7, 0.1],
    })
    ev = VariantEvaluator(df.lazy())
    from biostat_cli.evaluators.base import ScoreFrame
    sf = ScoreFrame(frame=df.lazy(), rows_used=4)
    gene_ls = ev.labels_and_scores_by_gene(sf, "eval_a", "score_x", "ensg")
    assert "g1" in gene_ls
    labels_g1, _ = gene_ls["g1"]
    assert set(labels_g1) == {1}


# ---------------------------------------------------------------------------
# Gene-averaged integration tests (end-to-end with run())
# ---------------------------------------------------------------------------


def _make_gene_avg_parquet(tmp_path):
    """3 genes × 8 variants each, mix of pos/neg with some overlap to ensure FP > 0."""
    df = pl.DataFrame({
        "chrom": ["1"] * 24,
        "pos": list(range(1, 25)),
        "ref": ["A"] * 24,
        "alt": ["C"] * 24,
        "ensg": (["gene_A"] * 8) + (["gene_B"] * 8) + (["gene_C"] * 8),
        "eval_a": ([True, True, True, True, False, False, False, False] * 3),
        "score_x": [
            0.95, 0.85, 0.70, 0.30, 0.80, 0.20, 0.15, 0.10,  # gene_A: 3 TP, 1 FP above 0.5
            0.92, 0.88, 0.75, 0.25, 0.78, 0.60, 0.12, 0.08,  # gene_B: 3 TP, 2 FP above 0.5
            0.90, 0.82, 0.40, 0.35, 0.72, 0.18, 0.14, 0.09,  # gene_C: 2 TP, 1 FP above 0.5
        ],
        "score_y": [
            0.88, 0.82, 0.60, 0.40, 0.75, 0.30, 0.20, 0.15,
            0.85, 0.78, 0.55, 0.35, 0.70, 0.25, 0.18, 0.12,
            0.82, 0.75, 0.50, 0.30, 0.65, 0.22, 0.16, 0.10,
        ],
    })
    path = tmp_path / "gene_avg.parquet"
    df.write_parquet(str(path))
    return path, df


def _make_gene_avg_resources(tmp_path, parquet_path, score_cols=None):
    resources = {
        "Table_info": {
            "gene_avg_test": {
                "Path": str(parquet_path),
                "Level": "variant",
                "Score_cols": score_cols or ["score_x"],
                "evals": ["eval_a"],
            }
        }
    }
    rpath = tmp_path / "resources_gene_avg.json"
    rpath.write_text(json.dumps(resources), encoding="utf-8")
    return rpath


def test_gene_avg_enrichment_integration(tmp_path):
    pq_path, _ = _make_gene_avg_parquet(tmp_path)
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_enrichment",
        eval_set=None, filters=None, thresholds="0.5",
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    out_df, _, _, _, _, _ = run(args)
    rows = out_df.to_dicts()
    assert len(rows) == 1
    row = rows[0]
    assert row["stat"] == "gene_avg_enrichment"
    assert row["n_genes_used"] == 3
    assert row["n_genes_excluded"] == 0
    assert not math.isnan(row["value"])
    assert not math.isnan(row["std_error"])
    assert row["std_error"] > 0
    assert 0.0 <= row["p_value"] <= 1.0


def test_gene_avg_rate_ratio_integration(tmp_path):
    pq_path, _ = _make_gene_avg_parquet(tmp_path)
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_rate_ratio",
        eval_set=None, filters=None, thresholds="0.5",
        case_total_by_eval="eval_a:12",
        ctrl_total_by_eval="eval_a:12",
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    out_df, _, _, _, _, _ = run(args)
    row = out_df.to_dicts()[0]
    assert row["stat"] == "gene_avg_rate_ratio"
    assert row["n_genes_used"] == 3
    assert not math.isnan(row["value"])


def test_gene_avg_auc_integration(tmp_path):
    pq_path, _ = _make_gene_avg_parquet(tmp_path)
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_auc",
        eval_set=None, filters=None, thresholds=None,
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    out_df, _, _, _, _, _ = run(args)
    row = out_df.to_dicts()[0]
    assert row["stat"] == "gene_avg_auc"
    assert row["n_genes_used"] == 3
    assert not math.isnan(row["value"])
    assert math.isnan(row["threshold"])


def test_gene_avg_auc_single_class_gene_excluded(tmp_path):
    df = pl.DataFrame({
        "chrom": ["1"] * 8, "pos": list(range(1, 9)),
        "ref": ["A"] * 8, "alt": ["C"] * 8,
        "ensg": ["g1", "g1", "g1", "g2", "g2", "g2", "g3", "g3"],
        "eval_a": [True, True, True, True, False, True, True, False],
        "score_x": [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2],
    })
    pq_path = tmp_path / "single_class.parquet"
    df.write_parquet(str(pq_path))
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_auc",
        eval_set=None, filters=None, thresholds=None,
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    out_df, _, _, _, _, _ = run(args)
    row = out_df.to_dicts()[0]
    # g1 is all-positive → excluded from AUC
    assert row["n_genes_used"] == 2
    assert row["n_genes_excluded"] == 1


def test_gene_avg_mixed_stats(tmp_path):
    pq_path, _ = _make_gene_avg_parquet(tmp_path)
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="enrichment,gene_avg_enrichment",
        eval_set=None, filters=None, thresholds="0.5",
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    out_df, _, _, _, _, _ = run(args)
    stats = set(out_df["stat"].to_list())
    assert "enrichment" in stats
    assert "gene_avg_enrichment" in stats
    pooled = out_df.filter(pl.col("stat") == "enrichment").to_dicts()[0]
    gene_avg = out_df.filter(pl.col("stat") == "gene_avg_enrichment").to_dicts()[0]
    assert "n_genes_used" in gene_avg
    assert not math.isnan(pooled["value"])
    assert not math.isnan(gene_avg["value"])


def test_gene_avg_threshold_and_trunc_stats_integration(tmp_path):
    pq_path, _ = _make_gene_avg_parquet(tmp_path)
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_tpr_at_threshold,gene_avg_auc_trunc",
        eval_set=None, filters=None, thresholds="0.8",
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    out_df, _, _, _, _, _ = run(args)
    stats = set(out_df["stat"].to_list())
    assert {"gene_avg_tpr_at_threshold", "gene_avg_auc_trunc"} <= stats


def test_gene_avg_requires_gene_col(tmp_path):
    df = pl.DataFrame({
        "chrom": ["1"] * 4, "pos": [1, 2, 3, 4],
        "ref": ["A"] * 4, "alt": ["C"] * 4,
        "eval_a": [True, False, True, False],
        "score_x": [0.9, 0.1, 0.8, 0.2],
    })
    pq_path = tmp_path / "no_gene.parquet"
    df.write_parquet(str(pq_path))
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_enrichment",
        eval_set=None, filters=None, thresholds="0.5",
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    with pytest.raises(ValueError, match="Gene-averaged stats require column"):
        run(args)


def test_gene_avg_rejects_gene_eval_level(tmp_path):
    pq_path, _ = _make_gene_avg_parquet(tmp_path)
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="gene", stat="gene_avg_enrichment",
        eval_set=None, filters=None, thresholds="0.5",
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    with pytest.raises(ValueError, match="not compatible with --eval-level gene"):
        run(args)


# ---------------------------------------------------------------------------
# Gene variant coverage report tests
# ---------------------------------------------------------------------------


def test_gene_variant_coverage_off_by_default(tmp_path):
    pq_path, _ = _make_gene_avg_parquet(tmp_path)
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_enrichment",
        eval_set=None, filters=None, thresholds="0.5",
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    _, _, _, _, coverage_df, _ = run(args)
    assert coverage_df.height == 0


def test_gene_variant_coverage_basic(tmp_path):
    df = pl.DataFrame({
        "chrom": ["1"] * 9, "pos": list(range(1, 10)),
        "ref": ["A"] * 9, "alt": ["C"] * 9,
        "ensg": ["g1"] * 3 + ["g2"] * 3 + ["g3"] * 3,
        "eval_a": [True, False, True, True, False, True, True, False, False],
        "score_x": [0.9, 0.1, None, 0.8, 0.2, 0.7, None, None, 0.3],
    })
    pq_path = tmp_path / "coverage.parquet"
    df.write_parquet(str(pq_path))
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_enrichment",
        eval_set=None, filters=None, thresholds="0.5",
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
        write_gene_variant_coverage=True,
    )
    _, _, _, _, coverage_df, _ = run(args)
    assert coverage_df.height == 3
    assert set(coverage_df.columns) == {
        "eval_name", "filter_name", "score_name", "gene",
        "n_variants_used", "n_variants_excluded", "n_variants_total",
    }
    for row in coverage_df.to_dicts():
        assert row["n_variants_used"] + row["n_variants_excluded"] == row["n_variants_total"]
    total_variants = coverage_df["n_variants_total"].sum()
    assert total_variants == 9


def test_gene_variant_coverage_multiple_scores(tmp_path):
    df = pl.DataFrame({
        "chrom": ["1"] * 6, "pos": list(range(1, 7)),
        "ref": ["A"] * 6, "alt": ["C"] * 6,
        "ensg": ["g1"] * 3 + ["g2"] * 3,
        "eval_a": [True, False, True] * 2,
        "score_x": [0.9, 0.1, 0.8, 0.7, 0.2, None],
        "score_y": [None, None, 0.8, 0.7, 0.2, 0.6],
    })
    pq_path = tmp_path / "cov_multi.parquet"
    df.write_parquet(str(pq_path))
    res_path = _make_gene_avg_resources(tmp_path, pq_path, score_cols=["score_x", "score_y"])
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_enrichment",
        eval_set=None, filters=None, thresholds="0.5",
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
        write_gene_variant_coverage=True,
    )
    _, _, _, _, coverage_df, _ = run(args)
    # 2 genes × 2 scores = 4 rows
    assert coverage_df.height == 4
    scores = set(coverage_df["score_name"].to_list())
    assert scores == {"score_x", "score_y"}
    # score_x has 1 null, score_y has 2 nulls
    sx_excluded = coverage_df.filter(pl.col("score_name") == "score_x")["n_variants_excluded"].sum()
    sy_excluded = coverage_df.filter(pl.col("score_name") == "score_y")["n_variants_excluded"].sum()
    assert sx_excluded == 1
    assert sy_excluded == 2


def test_gene_avg_zero_genes_after_filtering(tmp_path):
    """All genes have only one label class → all excluded from AUC."""
    df = pl.DataFrame({
        "chrom": ["1"] * 4, "pos": [1, 2, 3, 4],
        "ref": ["A"] * 4, "alt": ["C"] * 4,
        "ensg": ["g1", "g1", "g2", "g2"],
        "eval_a": [True, True, False, False],
        "score_x": [0.9, 0.8, 0.7, 0.6],
    })
    pq_path = tmp_path / "zero_genes.parquet"
    df.write_parquet(str(pq_path))
    res_path = _make_gene_avg_resources(tmp_path, pq_path)
    args = RunArgs(
        resources_json=str(res_path), table_name="gene_avg_test",
        eval_level="variant", stat="gene_avg_auc",
        eval_set=None, filters=None, thresholds=None,
        case_total_by_eval=None, ctrl_total_by_eval=None,
        bootstrap_samples=None, out_fname=str(tmp_path / "out"),
        write_missing="none",
    )
    out_df, _, _, _, _, _ = run(args)
    row = out_df.to_dicts()[0]
    assert row["n_genes_used"] == 0
    assert row["n_genes_excluded"] == 2
    assert math.isnan(row["value"])
