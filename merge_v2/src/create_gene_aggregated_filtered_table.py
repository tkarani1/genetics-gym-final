#!/usr/bin/env python3
"""Attach the variant- and gene-level filters to a gene-aggregated analysis table.

This is the gene-aggregated analogue of ``created_variant_scores_filtered_tables.py``.
That script filters the wide *variant-level* score table and therefore has to
join the variant->gene linker first (to obtain ``ensg``). The gene-aggregated
``full_analysis_tables/gene_aggregated/*_eval.parquet`` tables produced by
``create_analysis_tables.py`` **already** carry ``ensg`` (one row per
``(chrom, pos, ref, alt, ensg)``), so here the linker step is skipped and only
the filter columns are appended:

1. **Every variant-level filter** under ``../data/raw_data/filters/variant``
   (the ``variant_filters_*.tsv.gz`` files). Each contributes a single BOOLEAN
   membership column named after the file (``variant_filters_<name>.tsv.gz`` ->
   ``<name>``): TRUE when the variant is in that filter, FALSE otherwise.
2. **The gene/ensg-level filter** (``../data/raw_data/filters/ensg/ensg_filters.tsv``)
   joined on ``ensg``; its non-key columns are carried through with an ``ensg_``
   prefix (``uniprot_id`` stays VARCHAR, every other column is cast to BOOLEAN
   and is TRUE/NULL).

All joins are **LEFT OUTER** joins onto the (already gene-fanned) analysis
table, so the output has the **same row count** as the input -- no variant is
added or dropped, the filters only add columns. The output column order is the
input's columns, then the variant membership columns (sorted by filter file
name), then the ``ensg_``-prefixed gene-filter columns -- byte-for-byte the same
layout as the corresponding ``*_eval_filtered.parquet`` already published on GCS.

The heavy lifting (per-join SQL, the pipelined/checkpointed streaming ``COPY``
strategy, resume support) is reused verbatim from
``created_variant_scores_filtered_tables.py``; this module only changes which
steps make up the join chain (no linker) and forces a LEFT join.

Run from the ``src/`` directory::

    python create_gene_aggregated_filtered_table.py            # default outer table
    python create_gene_aggregated_filtered_table.py --input <path>_eval.parquet
    python create_gene_aggregated_filtered_table.py --memory-limit 12GB --threads 6
    python create_gene_aggregated_filtered_table.py --dry-run
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import duckdb

# Reuse the proven builders / runners from the variant-level filter script so
# the per-join SQL and the streaming/checkpoint engine are identical.
import created_variant_scores_filtered_tables as vf  # noqa: E402

SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # merge_v2/
GENE_AGG_DIR = (
    PROJECT_DIR / "data" / "processed_data" / "full_analysis_tables" / "gene_aggregated"
)
DEFAULT_INPUT = GENE_AGG_DIR / "variant_scores_all_outer_ensg_stats_eval.parquet"

FILTERED_SUFFIX = "_filtered"


def default_output(input_path: Path) -> Path:
    """``<stem>_eval.parquet`` -> ``<stem>_eval_filtered.parquet`` (same dir)."""
    return input_path.with_name(f"{input_path.stem}{FILTERED_SUFFIX}.parquet")


def build_filter_steps(
    con: duckdb.DuckDBPyConnection,
    variant_filters: list[Path],
    gene_filter: Path | None,
    compact: bool = False,
) -> list[vf.Step]:
    """The variant- and gene-filter join steps only (the linker step dropped).

    ``vf.build_steps`` returns ``[linker, variant_filter..., gene_filter]`` and
    also validates membership-column-name collisions; we reuse it and discard the
    leading linker step because the gene-aggregated input already has ``ensg``.
    """
    steps = vf.build_steps(con, vf.LINKER_PARQUET, variant_filters, gene_filter, compact)
    if steps and steps[0].kind == "linker":
        steps = steps[1:]
    return steps


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT,
        help=f"Gene-aggregated *_eval analysis table to filter (default: {DEFAULT_INPUT}).",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Output parquet path (default: <input stem>_filtered.parquet beside the input).",
    )
    parser.add_argument(
        "--variant-filters-dir", "--variant_filters_dir", dest="variant_filters_dir",
        type=Path, default=vf.VARIANT_FILTERS_DIR,
        help=f"Directory of variant_filters_*.tsv.gz files (default: {vf.VARIANT_FILTERS_DIR}).",
    )
    parser.add_argument(
        "--ensg-filters-dir", "--ensg_filters_dir", dest="ensg_filters_dir",
        type=Path, default=vf.ENSG_FILTERS_DIR,
        help=f"Directory of the gene/ensg-level filter TSV (default: {vf.ENSG_FILTERS_DIR}).",
    )
    parser.add_argument(
        "--compact-dtypes", action=argparse.BooleanOptionalAction, default=False,
        help="Match the compact numeric key profile of the input analysis table: "
             "encode each variant filter's chrom/pos/ref/alt the same way (so the "
             "joins match). Must match the flag the input score table was built "
             "with (the default gene-aggregated tables use the non-compact profile).",
    )
    parser.add_argument(
        "--checkpoint-every", type=int, default=10, metavar="N",
        help="Persist a Parquet checkpoint (and recycle the build DB) every N "
             "joins (default: 10). A crash loses at most N joins of work.",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Discard any existing output and checkpoints and start clean.",
    )
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Keep the checkpoint directory after a successful run (for debugging).",
    )
    parser.add_argument(
        "--memory-limit", default=None,
        help="DuckDB memory limit, e.g. '12GB' (default: DuckDB's ~80%% of RAM).",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="DuckDB worker threads (default: DuckDB auto).",
    )
    parser.add_argument(
        "--temp-dir", type=Path, default=None,
        help="Directory for DuckDB spill files during each batch "
             "(default: <output dir>/.duckdb_spill).",
    )
    parser.add_argument(
        "--compression", default="zstd",
        help="Parquet compression codec (default: zstd).",
    )
    parser.add_argument(
        "--row-group-size", type=int, default=512_000,
        help="Parquet row group size (default: 512000).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the join plan (and representative SQL) and exit.",
    )
    args = parser.parse_args()

    # The gene-aggregated table already carries ``ensg``; the filters only add
    # columns, so a LEFT join is the only correct flavor (preserves every row).
    args.join_type = "left"

    if args.checkpoint_every < 1:
        sys.exit("ERROR: --checkpoint-every must be >= 1.")

    if not args.input.exists():
        sys.exit(
            f"ERROR: input not found: {args.input}\n"
            f"       Build it with create_analysis_tables.py --groups gene first."
        )

    args.output = args.output or default_output(args.input)

    variant_filters = vf.variant_filter_files(args.variant_filters_dir)
    gene_filter = vf.gene_filter_file(args.ensg_filters_dir)
    if not variant_filters and gene_filter is None:
        sys.exit(
            "ERROR: no filter tables found under\n"
            f"  {args.variant_filters_dir}\n  {args.ensg_filters_dir}\n"
            "Run: python download_source_data.py --filter-tables all"
        )

    out_dir = args.output.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or (out_dir / ".duckdb_spill")
    ckpt_dir = vf.checkpoint_dir_for(args.output)

    # A short-lived connection only for schema introspection while planning.
    plan_con = duckdb.connect()
    try:
        steps = build_filter_steps(
            plan_con, variant_filters, gene_filter, args.compact_dtypes
        )
    finally:
        plan_con.close()
    n_steps = len(steps)
    last_idx = n_steps - 1

    print(f"Input:  {args.input}")
    print(f"Output: {args.output}")
    print(f"Variant-level filters: {len(variant_filters)} "
          f"(-> {len(variant_filters)} boolean membership columns)")
    print(f"Gene-level filter: {gene_filter.name if gene_filter else 'none'} "
          f"(-> '{vf.ENSG_FILTER_PREFIX}'-prefixed columns)")
    print(f"Join type: LEFT OUTER on {vf.VARIANT_KEYS} (variant filters) / "
          f"[{vf.ENSG_COL}] (gene filter); linker step skipped (ensg already present)")
    print(f"Total joins: {n_steps}  | checkpoint every {args.checkpoint_every}")
    print(f"Dtype profile: {'compact numeric' if args.compact_dtypes else 'default'}")
    print(f"Checkpoints: {ckpt_dir}")
    print(f"Spill dir: {temp_dir}\n")

    if args.dry_run:
        print("Join plan:")
        for i, s in enumerate(steps):
            print(f"  [{i}] {s.name}")
        print("\nRepresentative SQL:")
        if steps and steps[0].kind == "variant_filter":
            print("-- first variant filter join (onto the input table) --")
            print(steps[0].join_sql(vf.parquet_reader(args.input), args.join_type) + ";\n")
        if steps and steps[-1].kind == "gene_filter":
            print("-- gene filter join --")
            print(steps[-1].join_sql('"stage_prev"', args.join_type) + ";")
        return 0

    if args.overwrite:
        shutil.rmtree(ckpt_dir, ignore_errors=True)
        args.output.unlink(missing_ok=True)

    if args.output.exists():
        print("Output already exists; nothing to do (use --overwrite to rebuild).")
        return 0

    vf.validate_resume(ckpt_dir, args.input, args.join_type, steps, args.compact_dtypes)
    vf.write_plan(ckpt_dir, args.input, args.join_type, steps, args.compact_dtypes)

    next_idx, ckpt_source = vf.latest_checkpoint(ckpt_dir, n_steps)
    current_source = ckpt_source if ckpt_source is not None else args.input
    if ckpt_source is not None:
        print(f"Resuming from checkpoint {ckpt_source.name} (next join: step {next_idx}).\n")

    n_rows = 0
    idx = next_idx
    while idx <= last_idx:
        batch_end = min(idx + args.checkpoint_every - 1, last_idx)
        is_final = batch_end == last_idx
        out_path = args.output if is_final else vf.checkpoint_path(ckpt_dir, batch_end)
        print(f"  == batch: steps {idx}..{batch_end} "
              f"-> {'OUTPUT' if is_final else out_path.name} ==", flush=True)
        n_rows = vf.run_batch(
            temp_dir, Path(current_source), steps, idx, batch_end, out_path, args
        )
        print(f"     checkpoint written: {n_rows:,} rows -> {out_path.name}")

        prev = Path(current_source)
        if prev != args.input and prev.parent == ckpt_dir and prev.exists():
            prev.unlink()
        current_source = out_path
        idx = batch_end + 1

    if not args.keep_intermediates:
        shutil.rmtree(ckpt_dir, ignore_errors=True)
    shutil.rmtree(temp_dir, ignore_errors=True)

    final_con = duckdb.connect()
    try:
        cols = vf.column_names(final_con, vf.parquet_reader(args.output))
    finally:
        final_con.close()
    print(f"\nDone. Wrote {n_rows:,} rows x {len(cols)} columns to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
