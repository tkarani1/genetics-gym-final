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
from biostat_cli.stats.binary import enrichment, pairwise_enrichment, pairwise_rate_ratio, rate_ratio
from biostat_cli.stats.continuous import compute_auc, compute_auprc


def test_auc_and_auprc_basic():
    labels = [0, 0, 1, 1]
    scores = [0.1, 0.2, 0.8, 0.9]
    assert compute_auc(labels, scores) > 0.99
    assert compute_auprc(labels, scores) > 0.99


def test_binary_stats():
    cont = Contingency(tp=10, fp=5, tn=20, fn=15)
    enr = enrichment(cont)
    rr = rate_ratio(cont, case_total=200, ctrl_total=300)
    assert not math.isnan(enr.value)
    assert not math.isnan(enr.p_value)
    assert not math.isnan(rr.value)
    assert not math.isnan(rr.p_value)


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
        global_case_total=999.0,
        global_ctrl_total=888.0,
    )
    assert case_total == 111.0
    assert ctrl_total == 333.0


def test_resolve_eval_totals_fallbacks():
    # Falls back to table-level totals for matching eval.
    case_total, ctrl_total = _resolve_eval_totals(
        eval_col="eval_B",
        table_case_totals={"eval_B": 222.0},
        table_ctrl_totals={"eval_B": 444.0},
        cli_case_totals={},
        cli_ctrl_totals={},
        global_case_total=999.0,
        global_ctrl_total=888.0,
    )
    assert case_total == 222.0
    assert ctrl_total == 444.0

    # Falls back to global totals when eval-specific totals are absent.
    case_total, ctrl_total = _resolve_eval_totals(
        eval_col="eval_C",
        table_case_totals={},
        table_ctrl_totals={},
        cli_case_totals={},
        cli_ctrl_totals={},
        global_case_total=999.0,
        global_ctrl_total=888.0,
    )
    assert case_total == 999.0
    assert ctrl_total == 888.0


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
        case_total=None,
        ctrl_total=None,
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
        thresholds="0.8",
        case_total=None,
        ctrl_total=None,
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=None,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    base_df, _, _ = run(base_args)
    base_row = base_df.to_dicts()[0]
    assert math.isnan(base_row["std_error"])

    boot_args = RunArgs(
        resources_json=str(resources_path),
        table_name="toy",
        eval_level="variant",
        stat="enrichment",
        eval_set=None,
        filters=None,
        thresholds="0.8",
        case_total=None,
        ctrl_total=None,
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=20,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    boot_df, _, _ = run(boot_args)
    boot_row = boot_df.to_dicts()[0]
    assert boot_row["value"] == pytest.approx(base_row["value"], rel=0, abs=1e-12)
    assert boot_row["p_value"] == pytest.approx(base_row["p_value"], rel=0, abs=1e-12)
    assert not math.isnan(boot_row["std_error"])


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
        case_total=None,
        ctrl_total=None,
        case_total_by_eval=None,
        ctrl_total_by_eval=None,
        bootstrap_samples=20,
        out_fname=str(out_prefix),
        write_missing="none",
    )
    out_df, _, _ = run(args)
    rows = out_df.to_dicts()
    assert rows
    assert all("std_error" in row for row in rows)
    assert any(not math.isnan(row["std_error"]) for row in rows)
