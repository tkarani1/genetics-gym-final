#!/usr/bin/env python3
"""
Run biostat_cli pairwise metrics for is_pos_schema on the polyphen-anchor pairwise parquet.

Intended for `hailctl dataproc submit` with `--pyfiles` containing the `biostat_cli` package.
Writes TSV locally then `gsutil cp` to the requested GCS output prefix (directory or .tsv path).

Dataproc Hail images may not include Polars; we pip-install a pinned wheel before importing biostat_cli.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def _ensure_polars() -> None:
    try:
        import polars as pl  # noqa: F401
    except ModuleNotFoundError:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "-q", "polars==1.26.0"],
        )


_ensure_polars()

# Deploy layout: repo root contains eval/ with this file; pyfiles adds biostat_cli on PYTHONPATH.
_HERE = Path(__file__).resolve().parent
_REPO_EVAL = _HERE
if str(_REPO_EVAL) not in sys.path:
    sys.path.insert(0, str(_REPO_EVAL))

from biostat_cli import cli as bcli  # noqa: E402
from biostat_cli.io import write_tsv  # noqa: E402

DEFAULT_PAIRWISE_GS = (
    "gs://grohlicek/genetics_gym_vsm_all_content/vsm_all_2026_04_21/"
    "vsm_ensg_pairwise_polyphen_score_anchor_post.parquet"
)


def pairwise_score_columns_from_parquet(path: str) -> list[str]:
    """VSM percentile columns present in this build (Table_info Score_cols; schema varies by table version)."""
    import polars as pl

    names = pl.scan_parquet(path).collect_schema().names()
    return sorted(n for n in names if n.endswith("_percentile_with_anchor"))


def main() -> None:
    p = argparse.ArgumentParser(description="Dataproc one-off: pairwise stats for is_pos_schema")
    p.add_argument("--pairwise-parquet", default=os.environ.get("PAIRWISE_PARQUET_GS", DEFAULT_PAIRWISE_GS))
    p.add_argument(
        "--out-gs",
        default=(
            "gs://grohlicek/genetics_gym_vsm_all_content/vsm_all_2026_04_15/"
            "results/metrics_pairwise_is_pos_schema.tsv"
        ),
        help="gs://... path for final TSV (uploaded via gsutil cp)",
    )
    p.add_argument("--thresholds", default="0.90,0.95,0.98,0.99,0.995")
    p.add_argument("--bootstrap", type=int, default=None, help="Bootstrap replicates (>=2), omit to disable")
    p.add_argument("--stat", default="pairwise_enrichment,pairwise_rate_ratio,pairwise_auc,pairwise_auprc")
    p.add_argument("--eval-set", default="is_pos_schema")
    p.add_argument("--pvalue-method", default="fisher")
    args = p.parse_args()

    boot = args.bootstrap if args.bootstrap is not None and args.bootstrap >= 2 else None

    score_cols = pairwise_score_columns_from_parquet(args.pairwise_parquet)
    if not score_cols:
        raise SystemExit(f"No *_percentile_with_anchor columns found in {args.pairwise_parquet!r}")

    with tempfile.TemporaryDirectory() as td:
        resources_path = Path(td) / "resources_pairwise.json"
        out_prefix = str(Path(td) / "metrics_pairwise_schema")
        resources_path.write_text(
            json.dumps(
                {
                    "Table_info": {
                        "PAIRWISE": {
                            "Path": args.pairwise_parquet,
                            "Level": "variant",
                            "Score_cols": score_cols,
                            "evals": [args.eval_set],
                        }
                    }
                },
                indent=2,
            ),
            encoding="utf-8",
        )

        run_args = bcli.RunArgs(
            resources_json=str(resources_path),
            table_name="PAIRWISE",
            eval_level="variant",
            stat=args.stat,
            eval_set=args.eval_set,
            filters="none",
            thresholds=args.thresholds,
            case_total_by_eval=None,
            ctrl_total_by_eval=None,
            bootstrap_samples=boot,
            out_fname=out_prefix,
            write_missing="none",
            within_gene_percentile=False,
            pvalue_method=args.pvalue_method,
            vsm_comparison_method="fisher",
            gene_col="ensg",
            write_gene_variant_coverage=False,
        )

        print(f"[dataproc_pairwise_schema] parquet={args.pairwise_parquet}", flush=True)
        print(f"[dataproc_pairwise_schema] eval_set={args.eval_set} bootstrap={boot}", flush=True)

        out_df, timings, missing_df, vsm_cmp_df, cov_df = bcli.run(run_args)

        local_tsv = f"{out_prefix}.tsv"
        write_tsv(out_df, local_tsv)
        print(f"[dataproc_pairwise_schema] rows={out_df.height} wrote {local_tsv}", flush=True)
        for t in timings:
            print(f"  timing eval={t.get('eval_name')} filter={t.get('filter_name')} s={t.get('elapsed_seconds')}", flush=True)

        dest = args.out_gs
        if dest.endswith("/"):
            dest = dest + "metrics_pairwise_is_pos_schema.tsv"
        subprocess.check_call(["gsutil", "-q", "cp", local_tsv, dest])
        print(f"[dataproc_pairwise_schema] uploaded -> {dest}", flush=True)

        if vsm_cmp_df.height > 0:
            vsm_local = f"{out_prefix}_vsm_comparison.tsv"
            write_tsv(vsm_cmp_df, vsm_local)
            vsm_dest = dest.replace(".tsv", "_vsm_comparison.tsv")
            subprocess.check_call(["gsutil", "-q", "cp", vsm_local, vsm_dest])
            print(f"[dataproc_pairwise_schema] uploaded -> {vsm_dest}", flush=True)


if __name__ == "__main__":
    main()
