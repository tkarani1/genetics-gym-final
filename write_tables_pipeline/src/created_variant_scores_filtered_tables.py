#!/usr/bin/env python3
"""Attach the variant->gene linker and every filter table to the wide score table.

This builds ``variant_scores_outer_pre_percentile_filtered.parquet`` from
``../data/processed_data/scores/variant_scores_outer_pre_percentile.parquet`` by
a long chain of **FULL OUTER JOIN**s, in this order:

1. **Linker** (``../data/raw_data/linker/linker_all.parquet``) on the variant key
   ``(chrom, pos, ref, alt)``, in the *same manner* as
   ``create_variant_scores_gene_aggregation.py`` -- it adds the ``ensg`` column
   and, because a variant can map to several genes, fans the table out to one
   row per ``(variant, gene)``. (That script INNER-joins; here it is a full
   outer join so no variant is dropped.)
2. **Every variant-level filter** under ``../data/raw_data/filters/variant`` (the
   ``variant_filters_*.tsv.gz`` files downloaded by ``download_source_data.py``).
   Each such file *is* a set of variants -- membership in the file is the
   filter -- so each one contributes a single BOOLEAN membership column named
   after the file (``variant_filters_<name>.tsv.gz`` -> ``<name>``): TRUE when
   the variant is in that filter, FALSE otherwise.
3. **The gene/ensg-level filter** (``../data/raw_data/filters/ensg/ensg_filters.tsv``)
   on ``ensg``. Its non-key columns are carried through with an ``ensg_`` prefix
   (e.g. ``CATH_class_1`` -> ``ensg_CATH_class_1``) so they cannot collide with
   the case-insensitive variant-level filter names; ``uniprot_id`` stays a
   string, every other column is cast to BOOLEAN.

The result has the variant key, ``ensg``, every original score, one boolean per
variant-level filter (~65), and the prefixed gene-level filter columns.

Why this is structured so defensively
--------------------------------------
65+ sequential joins over a ~80M-row (post fan-out) table is extremely taxing.
The key observation is that every right-hand input here is *small*: each variant
filter is reduced to its distinct variant keys, the gene filter is ~one row per
gene, and only the linker is large (one big hash build, built once). So unlike a
join of many large tables, a batch of these joins can be **pipelined and streamed
in a single pass** without ever materializing the wide table. This script:

* Streams each batch of joins as **one ``COPY (...) TO`` straight from the
  previous checkpoint Parquet to the next** -- the wide table is never stored in
  a DuckDB table (which is what exhausted the disk: DuckDB's single build file
  does not release dropped-stage space back to the OS mid-connection). Hash
  tables for the small build sides stay in memory; anything large spills to
  ``--temp-dir`` under ``--memory-limit``.
* **Checkpoints to Parquet every ``--checkpoint-every`` joins** with a fresh
  connection per batch, so a crash only loses the joins since the last
  checkpoint. On restart the script finds the latest checkpoint and **resumes**
  from the next join.
* Keeps **at most two checkpoint files on disk** (the new one is written and
  its row count verified before the previous one is deleted), so disk use is
  bounded to roughly two copies of the (Parquet-compressed) wide table.

A ``plan.json`` next to the checkpoints records the exact join order; a resume
whose plan does not match the current inputs is refused (re-run with
``--overwrite`` to start clean).

Run from the ``src/`` directory (requires the score table, the linker, and the
filters to be present locally first)::

    python download_source_data.py --linker-table all --filter-tables all
    python create_percentile_score_tables.py --variant   # builds the input

    python created_variant_scores_filtered_tables.py
    python created_variant_scores_filtered_tables.py --memory-limit 20GB --threads 4
    python created_variant_scores_filtered_tables.py --checkpoint-every 5
    python created_variant_scores_filtered_tables.py --join-type left
    python created_variant_scores_filtered_tables.py --dry-run    # print the plan
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import duckdb

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: write_tables_pipeline/src/...)
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # write_tables_pipeline/
DATA_DIR = PROJECT_DIR / "data"
SCORES_DIR = DATA_DIR / "processed_data" / "scores"
LINKER_PARQUET = DATA_DIR / "raw_data" / "linker" / "linker_all.parquet"
VARIANT_FILTERS_DIR = DATA_DIR / "raw_data" / "filters" / "variant"
ENSG_FILTERS_DIR = DATA_DIR / "raw_data" / "filters" / "ensg"

DEFAULT_INPUT = SCORES_DIR / "variant_scores_outer_pre_percentile.parquet"
DEFAULT_OUTPUT_DIR = SCORES_DIR / "filtered"
DEFAULT_OUTPUT = DEFAULT_OUTPUT_DIR / "variant_scores_outer_pre_percentile_filtered.parquet"

# Variant key shared by the score table, the linker, and the variant filters.
VARIANT_KEYS = ["chrom", "pos", "ref", "alt"]
ENSG_COL = "ensg"

# Prefix applied to the gene-level filter's data columns so they cannot collide
# (case-insensitively, as DuckDB compares identifiers) with the variant filter
# membership columns -- e.g. variant ``cath_class_1`` vs gene ``CATH_class_1``.
ENSG_FILTER_PREFIX = "ensg_"
# Gene-level filter columns kept as VARCHAR; everything else is cast to BOOLEAN.
ENSG_STRING_COLS = {"uniprot_id"}

# Filenames of the variant filters carry this prefix; the membership column is
# the remainder once the prefix and the (double) extension are stripped.
VARIANT_FILTER_PREFIX = "variant_filters_"

# Null sentinels for the TSV sources (matches create_variant_scores_all_table.py).
CSV_NULL_VALUES = ["NA", "N/A", "N/a", "n/a", "na", "Na", "NaN", "nan", ""]

# Marker column used to detect filter membership across a FULL JOIN.
PRESENT_COL = "_filter_present"

PLAN_FILE = "plan.json"


def q(identifier: str) -> str:
    """Quote a SQL identifier (DuckDB double-quotes; e.g. 'ref' is reserved)."""
    return '"' + identifier.replace('"', '""') + '"'


def sql_str(value: str) -> str:
    """Quote a SQL string literal."""
    return "'" + value.replace("'", "''") + "'"


def parquet_reader(path: Path) -> str:
    """``read_parquet(...)`` for a file or a partitioned directory."""
    pattern = str(path / "**" / "*.parquet") if path.is_dir() else str(path)
    return f"read_parquet({sql_str(pattern)})"


def csv_reader(path: Path) -> str:
    """``read_csv(...)`` for a (optionally gzipped) TSV, read as all-VARCHAR.

    Quoting is disabled and every column is read as text; the key columns are
    re-cast explicitly downstream, so we never depend on type sniffing.
    """
    compression = "gzip" if path.suffix in (".bgz", ".gz") else "auto"
    nullstr = "[" + ", ".join(sql_str(v) for v in CSV_NULL_VALUES) + "]"
    return (
        "read_csv("
        f"{sql_str(str(path))}, "
        "delim='\t', header=true, quote='', "
        f"nullstr={nullstr}, "
        f"compression={sql_str(compression)}, "
        "all_varchar=true)"
    )


def column_names(con: duckdb.DuckDBPyConnection, reader: str) -> list[str]:
    """Column names exposed by a reader expression (metadata only, no scan)."""
    return [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {reader}").fetchall()]


# ---------------------------------------------------------------------------
# Discovery of the inputs that make up the join chain.
# ---------------------------------------------------------------------------
def variant_filter_files(directory: Path) -> list[Path]:
    """Sorted list of ``variant_filters_*.tsv.gz`` files (deterministic order)."""
    return sorted(directory.glob(f"{VARIANT_FILTER_PREFIX}*.tsv.gz"))


def gene_filter_file(directory: Path) -> Path | None:
    """The single gene/ensg-level filter TSV, if present."""
    candidates = sorted(directory.glob("*.tsv"))
    return candidates[0] if candidates else None


def variant_filter_colname(path: Path) -> str:
    """Membership column name for a variant filter file.

    ``variant_filters_cath_arch_3_90.tsv.gz`` -> ``cath_arch_3_90``.
    """
    name = path.name
    if name.startswith(VARIANT_FILTER_PREFIX):
        name = name[len(VARIANT_FILTER_PREFIX):]
    for suffix in (".tsv.gz", ".tsv.bgz", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


# ---------------------------------------------------------------------------
# Join SQL builders. Each returns the SELECT that turns ``prev_from`` (a table
# name or a read_*() expression) into the next wide table.
# ---------------------------------------------------------------------------
def keys_using() -> str:
    return "(" + ", ".join(q(k) for k in VARIANT_KEYS) + ")"


def linker_join_sql(prev_from: str, linker_reader: str, join_type: str) -> str:
    """FULL/LEFT OUTER JOIN the wide table to the linker, adding ``ensg``."""
    linker = (
        f"SELECT {', '.join(q(k) for k in VARIANT_KEYS)}, {q(ENSG_COL)} "
        f"FROM {linker_reader}"
    )
    return (
        f"SELECT * FROM {prev_from} "
        f"{join_type} JOIN (\n{linker}\n) USING {keys_using()}"
    )


def variant_filter_join_sql(
    prev_from: str, filter_reader: str, colname: str, join_type: str
) -> str:
    """FULL/LEFT OUTER JOIN one variant filter, adding a BOOLEAN membership column.

    The filter file is reduced to its DISTINCT variant keys (so a variant listed
    under several transcripts/genes cannot fan the join out) and tagged with a
    presence marker; non-members of an already-present row become FALSE.
    """
    casts = [
        f"TRY_CAST({q('chrom')} AS VARCHAR) AS {q('chrom')}",
        f"TRY_CAST({q('pos')} AS BIGINT) AS {q('pos')}",
        f"TRY_CAST({q('ref')} AS VARCHAR) AS {q('ref')}",
        f"TRY_CAST({q('alt')} AS VARCHAR) AS {q('alt')}",
    ]
    nonnull = " AND ".join(f"{q(k)} IS NOT NULL" for k in VARIANT_KEYS)
    keyed = (
        f"SELECT DISTINCT {', '.join(q(k) for k in VARIANT_KEYS)}, "
        f"TRUE AS {q(PRESENT_COL)} FROM (\n"
        f"  SELECT {', '.join(casts)} FROM {filter_reader}\n"
        f") WHERE {nonnull}"
    )
    return (
        f"SELECT * EXCLUDE ({q(PRESENT_COL)}), "
        f"COALESCE({q(PRESENT_COL)}, FALSE) AS {q(colname)} "
        f"FROM {prev_from} "
        f"{join_type} JOIN (\n{keyed}\n) USING {keys_using()}"
    )


def gene_filter_join_sql(
    prev_from: str,
    gene_reader: str,
    gene_value_cols: list[str],
    join_type: str,
) -> str:
    """FULL/LEFT OUTER JOIN the gene-level filter on ``ensg``.

    The filter is deduped to one row per ``ensg`` (``max`` per column, safe even
    if the source is already unique). Boolean-flag columns are cast to BOOLEAN;
    string columns are kept; every output column is ``ensg_``-prefixed.
    """
    agg_cols = []
    for c in gene_value_cols:
        out = q(ENSG_FILTER_PREFIX + c)
        if c in ENSG_STRING_COLS:
            agg_cols.append(f"max({q(c)}) AS {out}")
        else:
            agg_cols.append(f"max(TRY_CAST({q(c)} AS BOOLEAN)) AS {out}")
    deduped = (
        f"SELECT {q(ENSG_COL)}, {', '.join(agg_cols)} "
        f"FROM {gene_reader} WHERE {q(ENSG_COL)} IS NOT NULL "
        f"GROUP BY {q(ENSG_COL)}"
    )
    return (
        f"SELECT * FROM {prev_from} "
        f"{join_type} JOIN (\n{deduped}\n) USING ({q(ENSG_COL)})"
    )


# ---------------------------------------------------------------------------
# Step plan.
# ---------------------------------------------------------------------------
class Step:
    """One join in the chain."""

    def __init__(self, name: str, kind: str, reader: str, extra: dict | None = None):
        self.name = name        # short label, also used for plan validation
        self.kind = kind        # 'linker' | 'variant_filter' | 'gene_filter'
        self.reader = reader     # read_*() expression for the right-hand input
        self.extra = extra or {}

    def join_sql(self, prev_from: str, join_type: str) -> str:
        if self.kind == "linker":
            return linker_join_sql(prev_from, self.reader, join_type)
        if self.kind == "variant_filter":
            return variant_filter_join_sql(
                prev_from, self.reader, self.extra["colname"], join_type
            )
        if self.kind == "gene_filter":
            return gene_filter_join_sql(
                prev_from, self.reader, self.extra["value_cols"], join_type
            )
        raise ValueError(f"unknown step kind: {self.kind}")


def build_steps(
    con: duckdb.DuckDBPyConnection,
    linker: Path,
    variant_filters: list[Path],
    gene_filter: Path | None,
) -> list[Step]:
    """Assemble the ordered list of joins and validate column names up front."""
    steps: list[Step] = [Step("linker", "linker", parquet_reader(linker))]

    seen: dict[str, str] = {}
    for vf in variant_filters:
        col = variant_filter_colname(vf)
        lower = col.lower()
        if lower in seen:
            sys.exit(
                f"ERROR: variant filter membership column '{col}' (from {vf.name}) "
                f"collides with '{seen[lower]}'. Filter names must be unique."
            )
        seen[lower] = col
        steps.append(
            Step(f"variant:{col}", "variant_filter", csv_reader(vf), {"colname": col})
        )

    if gene_filter is not None:
        cols = column_names(con, csv_reader(gene_filter))
        value_cols = [c for c in cols if c != ENSG_COL]
        if ENSG_COL not in cols:
            sys.exit(
                f"ERROR: gene-level filter {gene_filter.name} has no '{ENSG_COL}' "
                f"column to join on. Found: {cols}"
            )
        steps.append(
            Step("gene", "gene_filter", csv_reader(gene_filter), {"value_cols": value_cols})
        )
    return steps


# ---------------------------------------------------------------------------
# Checkpointing / resume.
# ---------------------------------------------------------------------------
def checkpoint_dir_for(output: Path) -> Path:
    return output.parent / f".ckpt_{output.stem}"


def checkpoint_path(ckpt_dir: Path, step_idx: int) -> Path:
    return ckpt_dir / f"after_step_{step_idx:03d}.parquet"


def write_plan(ckpt_dir: Path, base_input: Path, join_type: str, steps: list[Step]) -> None:
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    plan = {
        "base_input": str(base_input),
        "join_type": join_type,
        "steps": [s.name for s in steps],
    }
    (ckpt_dir / PLAN_FILE).write_text(json.dumps(plan, indent=2))


def load_plan(ckpt_dir: Path) -> dict | None:
    p = ckpt_dir / PLAN_FILE
    if not p.exists():
        return None
    return json.loads(p.read_text())


def validate_resume(
    ckpt_dir: Path, base_input: Path, join_type: str, steps: list[Step]
) -> None:
    """Refuse to resume against a plan that no longer matches the inputs."""
    plan = load_plan(ckpt_dir)
    if plan is None:
        return
    current = {
        "base_input": str(base_input),
        "join_type": join_type,
        "steps": [s.name for s in steps],
    }
    if plan != current:
        sys.exit(
            "ERROR: existing checkpoints in\n"
            f"  {ckpt_dir}\n"
            "were built with a different plan (inputs / join order / join type "
            "changed). Re-run with --overwrite to start clean, or restore the "
            "original inputs."
        )


def latest_checkpoint(ckpt_dir: Path, n_steps: int) -> tuple[int, Path | None]:
    """Return (next_step_index, source_parquet) given existing checkpoints."""
    best = -1
    for idx in range(n_steps):
        if checkpoint_path(ckpt_dir, idx).exists():
            best = idx
    if best == -1:
        return 0, None
    return best + 1, checkpoint_path(ckpt_dir, best)


def configure(
    con: duckdb.DuckDBPyConnection,
    memory_limit: str | None,
    threads: int | None,
    temp_dir: Path,
) -> None:
    """Apply memory / spill settings so each large join stays out-of-core."""
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET temp_directory = {sql_str(str(temp_dir))}")
    con.execute("SET preserve_insertion_order = false")
    if memory_limit:
        con.execute(f"SET memory_limit = {sql_str(memory_limit)}")
    if threads:
        con.execute(f"SET threads = {int(threads)}")


def copy_sql(source_sql: str, out_path: Path, compression: str, row_group_size: int) -> str:
    return (
        f"COPY ({source_sql}) TO {sql_str(str(out_path))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(compression)}, "
        f"ROW_GROUP_SIZE {int(row_group_size)})"
    )


def run_batch(
    temp_dir: Path,
    source_parquet: Path,
    steps: list[Step],
    start_idx: int,
    end_idx: int,
    out_path: Path,
    args: argparse.Namespace,
) -> int:
    """Stream steps[start_idx..end_idx] as one ``COPY`` and return the row count.

    The batch's joins are nested into a single pipelined query that reads the
    previous checkpoint Parquet and writes the next one -- the wide table is
    never materialized in a DuckDB table, so the build file cannot grow. A fresh
    in-memory connection (spilling to ``temp_dir``) is used per batch.
    """
    shutil.rmtree(temp_dir, ignore_errors=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    # Nest each join as a subquery feeding the next; the build sides are small
    # (distinct filter keys / per-gene rows / the single linker), so DuckDB can
    # pipeline the whole batch in one streaming pass to Parquet.
    from_expr = parquet_reader(source_parquet)
    select_sql = f"SELECT * FROM {from_expr}"
    labels = []
    for j in range(start_idx, end_idx + 1):
        select_sql = steps[j].join_sql(from_expr, args.join_type)
        from_expr = f"(\n{select_sql}\n)"
        labels.append(f"[{j}] {steps[j].name}")
    print("    joins: " + ", ".join(labels), flush=True)

    con = duckdb.connect()
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        print(f"    streaming -> {out_path.name} ...", flush=True)
        con.execute(copy_sql(select_sql, out_path, args.compression, args.row_group_size))
        n_rows = con.execute(
            f"SELECT count(*) FROM {parquet_reader(out_path)}"
        ).fetchone()[0]
    finally:
        con.close()
    shutil.rmtree(temp_dir, ignore_errors=True)
    return n_rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT,
        help=f"Wide variant score table to filter (default: {DEFAULT_INPUT}).",
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help=f"Output parquet path (default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--linker", type=Path, default=LINKER_PARQUET,
        help=f"variant->gene linker parquet (default: {LINKER_PARQUET}).",
    )
    parser.add_argument(
        "--variant-filters-dir", "--variant_filters_dir", dest="variant_filters_dir",
        type=Path, default=VARIANT_FILTERS_DIR,
        help=f"Directory of variant_filters_*.tsv.gz files (default: {VARIANT_FILTERS_DIR}).",
    )
    parser.add_argument(
        "--ensg-filters-dir", "--ensg_filters_dir", dest="ensg_filters_dir",
        type=Path, default=ENSG_FILTERS_DIR,
        help=f"Directory of the gene/ensg-level filter TSV (default: {ENSG_FILTERS_DIR}).",
    )
    parser.add_argument(
        "--join-type", choices=["full", "left"], default="full",
        help="Outer join flavor. 'full' (default) keeps variants present only in "
             "a filter/linker; 'left' keeps only the score table's variants.",
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
        help="DuckDB memory limit, e.g. '20GB' (default: DuckDB's ~80%% of RAM).",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="DuckDB worker threads (default: DuckDB auto).",
    )
    parser.add_argument(
        "--temp-dir", type=Path, default=None,
        help="Directory for DuckDB spill files during each batch "
             "(default: <output dir>/.duckdb_spill). Needs free space for the "
             "largest join's spill (the linker step).",
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
        help="Print the join plan (and a few representative SQL statements) and exit.",
    )
    args = parser.parse_args()

    if args.checkpoint_every < 1:
        sys.exit("ERROR: --checkpoint-every must be >= 1.")

    if not args.input.exists():
        sys.exit(
            f"ERROR: input not found: {args.input}\n"
            f"       Build it with create_variant_scores_all_table.py first."
        )
    if not args.linker.exists():
        sys.exit(
            f"ERROR: linker not found: {args.linker}\n"
            f"       Run: python download_source_data.py --linker-table all"
        )

    variant_filters = variant_filter_files(args.variant_filters_dir)
    gene_filter = gene_filter_file(args.ensg_filters_dir)
    if not variant_filters and gene_filter is None:
        sys.exit(
            "ERROR: no filter tables found under\n"
            f"  {args.variant_filters_dir}\n  {args.ensg_filters_dir}\n"
            "Run: python download_source_data.py --filter-tables all"
        )

    out_dir = args.output.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or (out_dir / ".duckdb_spill")
    ckpt_dir = checkpoint_dir_for(args.output)

    # A short-lived connection just for schema introspection while planning.
    plan_con = duckdb.connect()
    try:
        steps = build_steps(plan_con, args.linker, variant_filters, gene_filter)
    finally:
        plan_con.close()
    n_steps = len(steps)
    last_idx = n_steps - 1

    print(f"Input:  {args.input}")
    print(f"Output: {args.output}")
    print(f"Linker: {args.linker}")
    print(f"Variant-level filters: {len(variant_filters)} "
          f"(-> {len(variant_filters)} boolean membership columns)")
    print(f"Gene-level filter: {gene_filter.name if gene_filter else 'none'}")
    print(f"Join type: {args.join_type.upper()} OUTER on "
          f"{VARIANT_KEYS} (filters) / [{ENSG_COL}] (gene filter)")
    print(f"Total joins: {n_steps}  | checkpoint every {args.checkpoint_every}")
    print(f"Checkpoints: {ckpt_dir}")
    print(f"Spill dir: {temp_dir}\n")

    if args.dry_run:
        print("Join plan:")
        for i, s in enumerate(steps):
            print(f"  [{i}] {s.name}")
        print("\nRepresentative SQL:")
        print("-- linker join --")
        print(steps[0].join_sql(parquet_reader(args.input), args.join_type) + ";\n")
        if len(steps) > 1 and steps[1].kind == "variant_filter":
            print("-- a variant filter join --")
            print(steps[1].join_sql('"stage_0"', args.join_type) + ";\n")
        if steps[-1].kind == "gene_filter":
            print("-- gene filter join --")
            print(steps[-1].join_sql('"stage_prev"', args.join_type) + ";")
        return 0

    if args.overwrite:
        shutil.rmtree(ckpt_dir, ignore_errors=True)
        args.output.unlink(missing_ok=True)

    if args.output.exists():
        print("Output already exists; nothing to do (use --overwrite to rebuild).")
        return 0

    validate_resume(ckpt_dir, args.input, args.join_type, steps)
    write_plan(ckpt_dir, args.input, args.join_type, steps)

    next_idx, ckpt_source = latest_checkpoint(ckpt_dir, n_steps)
    current_source = ckpt_source if ckpt_source is not None else args.input
    if ckpt_source is not None:
        print(f"Resuming from checkpoint {ckpt_source.name} (next join: step {next_idx}).\n")

    n_rows = 0
    idx = next_idx
    while idx <= last_idx:
        batch_end = min(idx + args.checkpoint_every - 1, last_idx)
        is_final = batch_end == last_idx
        out_path = args.output if is_final else checkpoint_path(ckpt_dir, batch_end)
        print(f"  == batch: steps {idx}..{batch_end} "
              f"-> {'OUTPUT' if is_final else out_path.name} ==", flush=True)
        n_rows = run_batch(
            temp_dir, Path(current_source), steps, idx, batch_end, out_path, args
        )
        print(f"     checkpoint written: {n_rows:,} rows -> {out_path.name}")

        # Now that the new checkpoint/output exists, drop the previous checkpoint
        # parquet so at most two checkpoints are ever on disk at once.
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
        cols = column_names(final_con, parquet_reader(args.output))
    finally:
        final_con.close()
    print(f"\nDone. Wrote {n_rows:,} rows x {len(cols)} columns to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
