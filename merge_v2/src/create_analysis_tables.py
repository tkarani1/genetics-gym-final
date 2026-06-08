#!/usr/bin/env python3
"""Join the score tables to the eval tables into ``full_analysis_tables``.

This is the final assembly step of the merge_v2 pipeline. It takes most of the
score Parquet files under ``../data/processed_data/scores`` and LEFT-JOINs them
onto the merged eval tables under ``../data/processed_data/evals``, writing the
results to ``../data/processed_data/full_analysis_tables``. Every output is
named ``{original score file stem}_eval.parquet``.

What gets joined
----------------
Three groups of score tables are processed (select a subset with ``--groups``):

* ``variant`` -- each top-level ``scores/variant_scores_*.parquet`` table
  (the ``_all_*`` and ``_*_percentile`` wide tables) LEFT-JOINed onto
  ``evals/variant_evals_all.parquet`` on ``(chrom, pos, ref, alt)``. Outputs are
  written flat in ``full_analysis_tables``.

* ``pairwise`` -- every file under ``scores/pairwise/{pairwise_raw,pairwise_pre,
  pairwise_post}`` LEFT-JOINed onto ``variant_evals_all.parquet`` on the same
  variant key. The ``pairwise/<flavor>/`` sub-structure is mirrored under
  ``full_analysis_tables/pairwise/<flavor>/``.

* ``gene`` -- the two ``scores/gene_aggregated/*_ensg_stats.parquet`` tables.
  These carry **both** the variant key ``(chrom, pos, ref, alt)`` and the gene
  key ``ensg``, so they are LEFT-JOINed onto ``variant_evals_all.parquet`` (on
  the variant key) **and** onto ``evals/ensg_evals_all.parquet`` (on ``ensg``).
  Outputs go under ``full_analysis_tables/gene_aggregated/``. For each gene
  table two outputs are written: ``{stem}_eval.parquet`` (the join as-is) and
  ``{stem}_eval_deduped.parquet`` (see below).

The ``filtered`` score tables are intentionally **not** processed.

Join semantics: why LEFT JOIN is safe
-------------------------------------
Both eval tables are built (by ``create_variant_eval_all_table.py`` /
``create_ensg_eval_all_table.py``) with exactly **one row per key** -- one per
``(chrom, pos, ref, alt)`` for the variant evals and one per ``ensg`` for the
gene evals. Joining a unique-keyed table onto the score table therefore cannot
duplicate ("fan out") any score row, and the LEFT direction keeps every row of
the score table (eval columns are simply NULL where there is no match). So no
data is lost from the original score tables and no extra rows are introduced.
As a guard, each output's row count is compared to its input's; any increase
(which would mean an eval table was not actually unique on its key) is reported
loudly.

Gene-aggregated dedup output
----------------------------
A gene-aggregated score table has one row per ``(variant, gene)`` pair, so a
variant that maps to several genes appears in several rows that share the same
``(chrom, pos, ref, alt)``. The ``_eval_deduped`` output keeps only the rows
whose ``(chrom, pos, ref, alt)`` is **unique** in the table: any key that occurs
more than once has *all* of its rows dropped (an unambiguous resolution, since
there is no principled single row to keep among the gene copies). It is computed
from the just-written ``_eval`` file via a per-key count + semi-join, so the
join is not recomputed and only the four key columns are scanned to find the
duplicates.

Engine / memory strategy
------------------------
Same DuckDB out-of-core philosophy as the upstream scripts. Each ``_eval``
output is a single streaming ``COPY (...) TO`` over a ``read_parquet`` LEFT JOIN
(DuckDB builds a hash table on the small eval side and probes it with the score
table, spilling to ``--temp-dir`` under ``--memory-limit``). The dedup step
first materializes the set of unique keys (``GROUP BY ... HAVING count(*) = 1``
over the four key columns only) into a TEMP table, then writes the deduped file
with a streaming semi-join. Existing outputs are skipped (resumable) unless
``--overwrite`` is set.

Run from the ``src/`` directory::

    python create_analysis_tables.py
    python create_analysis_tables.py --groups variant pairwise
    python create_analysis_tables.py --groups gene --memory-limit 20GB
    python create_analysis_tables.py --dry-run
    python create_analysis_tables.py --overwrite
"""

from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

import duckdb

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: merge_v2/src/...)
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # merge_v2/
PROC_DIR = PROJECT_DIR / "data" / "processed_data"
SCORES_DIR = PROC_DIR / "scores"
EVALS_DIR = PROC_DIR / "evals"
DEFAULT_OUTPUT_DIR = PROC_DIR / "full_analysis_tables"

VARIANT_EVALS = EVALS_DIR / "variant_evals_all.parquet"
ENSG_EVALS = EVALS_DIR / "ensg_evals_all.parquet"

# Join keys.
VARIANT_KEYS = ["chrom", "pos", "ref", "alt"]
ENSG_KEYS = ["ensg"]

# Pairwise flavor subdirectories (mirrored into the output tree).
PAIRWISE_FLAVORS = ["pairwise_raw", "pairwise_pre", "pairwise_post"]

# Only the per-gene *stats* tables are joined for the gene group.
GENE_INPUT_GLOB = "*_ensg_stats.parquet"

# Output naming.
EVAL_SUFFIX = "_eval"
DEDUP_SUFFIX = "_eval_deduped"

GROUPS = ("variant", "pairwise", "gene")


def q(identifier: str) -> str:
    """Quote a SQL identifier (DuckDB double-quotes; e.g. 'ref' is reserved)."""
    return '"' + identifier.replace('"', '""') + '"'


def sql_str(value: str) -> str:
    """Quote a SQL string literal."""
    return "'" + value.replace("'", "''") + "'"


def reader_sql(resolved: Path) -> str:
    """Build the DuckDB ``read_parquet(...)`` expression for a parquet input."""
    pattern = str(resolved / "**" / "*.parquet") if resolved.is_dir() else str(resolved)
    return f"read_parquet({sql_str(pattern)})"


def column_names(con: duckdb.DuckDBPyConnection, reader: str) -> list[str]:
    """Return the column names exposed by a reader expression (no full scan)."""
    return [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {reader}").fetchall()]


@dataclass
class EvalJoin:
    """One eval table to LEFT-JOIN onto a score table."""

    alias: str
    path: Path
    keys: list[str]
    value_cols: list[str] = field(default_factory=list)  # non-key columns


@dataclass
class Job:
    """One score table to join, with its output location."""

    input_path: Path
    out_dir: Path
    joins: list[EvalJoin]
    dedup: bool = False  # also emit a {stem}_eval_deduped.parquet


def copy_sql(source_sql: str, out_path: Path, compression: str, row_group_size: int) -> str:
    """Build a ``COPY (...) TO '<path>' (FORMAT PARQUET, ...)`` statement."""
    return (
        f"COPY ({source_sql}) TO {sql_str(str(out_path))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(compression)}, "
        f"ROW_GROUP_SIZE {int(row_group_size)})"
    )


def eval_value_columns(
    con: duckdb.DuckDBPyConnection, eval_path: Path, keys: list[str]
) -> list[str]:
    """Return an eval table's non-key columns (the columns it contributes)."""
    if not eval_path.exists():
        sys.exit(
            f"ERROR: eval table not found: {eval_path}\n"
            f"       Have you built it (create_variant_eval_all_table.py / "
            f"create_ensg_eval_all_table.py)?"
        )
    cols = column_names(con, reader_sql(eval_path.resolve()))
    missing = [k for k in keys if k not in cols]
    if missing:
        sys.exit(
            f"ERROR: eval table {eval_path.name} is missing key column(s) "
            f"{missing}. Found: {cols}"
        )
    return [c for c in cols if c not in keys]


def validate_job(con: duckdb.DuckDBPyConnection, job: Job) -> list[str]:
    """Validate a job's keys/columns; return the score table's column names.

    Ensures every join key exists in the score table and that no eval value
    column collides with a score column or with another eval's value column
    (a collision would yield duplicate output column names).
    """
    score_cols = column_names(con, reader_sql(job.input_path.resolve()))
    score_set = set(score_cols)

    for join in job.joins:
        missing = [k for k in join.keys if k not in score_set]
        if missing:
            sys.exit(
                f"ERROR: {job.input_path.name} is missing join key(s) {missing} "
                f"needed to join {join.path.name}. Found columns: {score_cols}"
            )

    seen: dict[str, str] = {c: "score table" for c in score_cols}
    for join in job.joins:
        for c in join.value_cols:
            if c in seen:
                sys.exit(
                    f"ERROR: column name collision on '{c}' between "
                    f"{join.path.name} and {seen[c]} when joining "
                    f"{job.input_path.name}."
                )
            seen[c] = join.path.name
    return score_cols


def join_select_sql(job: Job) -> str:
    """Build the LEFT-JOIN projection for a job.

    Keeps every column of the score table (``s.*``, original order) and appends
    each eval table's value columns (qualified by its alias). Explicit ``ON``
    predicates (rather than ``USING``) keep the score key columns first and the
    schema order stable.
    """
    score_reader = reader_sql(job.input_path.resolve())
    select_parts = ["s.*"]
    from_parts = [f"{score_reader} AS s"]
    for join in job.joins:
        reader = reader_sql(join.path.resolve())
        on = " AND ".join(f"s.{q(k)} = {join.alias}.{q(k)}" for k in join.keys)
        from_parts.append(f"LEFT JOIN {reader} AS {join.alias} ON {on}")
        select_parts += [f"{join.alias}.{q(c)} AS {q(c)}" for c in join.value_cols]
    return "SELECT " + ", ".join(select_parts) + "\nFROM " + "\n".join(from_parts)


def dedup_select_sql(eval_path: Path, keys: list[str], unique_keys_table: str) -> str:
    """Semi-join the written ``_eval`` file to the set of unique-key rows.

    ``unique_keys_table`` is expected to hold exactly the key tuples that occur
    once in ``eval_path``; the semi-join therefore keeps only rows whose key is
    unique (dropping *all* rows of any duplicated key).
    """
    reader = reader_sql(eval_path.resolve())
    on = " AND ".join(f"j.{q(k)} = u.{q(k)}" for k in keys)
    return (
        f"SELECT j.* FROM {reader} AS j "
        f"SEMI JOIN {q(unique_keys_table)} AS u ON {on}"
    )


def count_rows(con: duckdb.DuckDBPyConnection, path: Path) -> int:
    return con.execute(
        f"SELECT count(*) FROM read_parquet({sql_str(str(path))})"
    ).fetchone()[0]


def configure(
    con: duckdb.DuckDBPyConnection,
    memory_limit: str | None,
    threads: int | None,
    temp_dir: Path,
) -> None:
    """Apply memory / spill settings so large joins stay out-of-core."""
    temp_dir.mkdir(parents=True, exist_ok=True)
    # temp_directory MUST be set for an in-memory connection to spill to disk.
    con.execute(f"SET temp_directory = {sql_str(str(temp_dir))}")
    con.execute("SET preserve_insertion_order = false")
    if memory_limit:
        con.execute(f"SET memory_limit = {sql_str(memory_limit)}")
    if threads:
        con.execute(f"SET threads = {int(threads)}")


def build_jobs(
    groups: list[str],
    scores_dir: Path,
    out_dir: Path,
    variant_evals: Path,
    ensg_evals: Path,
) -> list[Job]:
    """Discover the score tables to process and assemble the join jobs."""
    jobs: list[Job] = []

    if "variant" in groups:
        # Top-level wide variant score tables (filtered/ lives in a subdir and
        # is excluded by globbing only the top level).
        for p in sorted(scores_dir.glob("variant_scores_*.parquet")):
            jobs.append(
                Job(
                    input_path=p,
                    out_dir=out_dir,
                    joins=[EvalJoin("ve", variant_evals, VARIANT_KEYS)],
                )
            )

    if "pairwise" in groups:
        for flavor in PAIRWISE_FLAVORS:
            flavor_dir = scores_dir / "pairwise" / flavor
            if not flavor_dir.is_dir():
                continue
            for p in sorted(flavor_dir.glob("*.parquet")):
                jobs.append(
                    Job(
                        input_path=p,
                        out_dir=out_dir / "pairwise" / flavor,
                        joins=[EvalJoin("ve", variant_evals, VARIANT_KEYS)],
                    )
                )

    if "gene" in groups:
        gene_dir = scores_dir / "gene_aggregated"
        for p in sorted(gene_dir.glob(GENE_INPUT_GLOB)):
            jobs.append(
                Job(
                    input_path=p,
                    out_dir=out_dir / "gene_aggregated",
                    joins=[
                        EvalJoin("ve", variant_evals, VARIANT_KEYS),
                        EvalJoin("ee", ensg_evals, ENSG_KEYS),
                    ],
                    dedup=True,
                )
            )

    return jobs


def process_job(con: duckdb.DuckDBPyConnection, job: Job, args: argparse.Namespace) -> None:
    """Run (or print, when --dry-run) one join job and its optional dedup."""
    stem = job.input_path.stem
    eval_out = job.out_dir / f"{stem}{EVAL_SUFFIX}.parquet"
    dedup_out = job.out_dir / f"{stem}{DEDUP_SUFFIX}.parquet"

    # Resolve eval value columns and validate keys / collisions.
    for join in job.joins:
        join.value_cols = eval_value_columns(con, join.path, join.keys)
    validate_job(con, job)

    join_sql = join_select_sql(job)

    print(f"\n=== {job.input_path.name} ===")
    print(f"Input:  {job.input_path}")
    for join in job.joins:
        print(f"Join:   LEFT JOIN {join.path.name} ON {join.keys} "
              f"(+{len(join.value_cols)} eval cols)")
    print(f"Output: {eval_out}")
    if job.dedup:
        print(f"Dedup:  {dedup_out}  (keep only unique {VARIANT_KEYS} rows)")

    if args.dry_run:
        print("  (dry run -- nothing written)")
        print(copy_sql(join_sql, eval_out, args.compression, args.row_group_size) + ";")
        if job.dedup:
            print("  [dedup] computed at run time from the written _eval file.")
        return

    job.out_dir.mkdir(parents=True, exist_ok=True)

    if eval_out.exists() and not args.overwrite:
        print("  skip _eval (output exists; use --overwrite to rebuild)")
    else:
        print("  joining + writing _eval parquet ...", flush=True)
        con.execute(copy_sql(join_sql, eval_out, args.compression, args.row_group_size))
        n_in = count_rows(con, job.input_path)
        n_out = count_rows(con, eval_out)
        print(f"  done: {n_in:,} input rows -> {n_out:,} output rows")
        if n_out != n_in:
            print(
                f"  WARNING: row count changed ({n_in:,} -> {n_out:,}). A LEFT "
                f"JOIN onto a unique-keyed eval table should preserve the row "
                f"count exactly; an increase means an eval table is not unique "
                f"on its key (unexpected fan-out).",
                flush=True,
            )

    if not job.dedup:
        return

    if dedup_out.exists() and not args.overwrite:
        print("  skip _eval_deduped (output exists; use --overwrite to rebuild)")
        return

    # Drop every row of any (chrom,pos,ref,alt) that occurs more than once.
    keys_sql = ", ".join(q(k) for k in VARIANT_KEYS)
    unique_keys_table = "_unique_keys"
    con.execute(
        f"CREATE OR REPLACE TEMP TABLE {q(unique_keys_table)} AS "
        f"SELECT {keys_sql} FROM {reader_sql(eval_out.resolve())} "
        f"GROUP BY {keys_sql} HAVING count(*) = 1"
    )
    print("  writing _eval_deduped parquet ...", flush=True)
    con.execute(
        copy_sql(
            dedup_select_sql(eval_out, VARIANT_KEYS, unique_keys_table),
            dedup_out,
            args.compression,
            args.row_group_size,
        )
    )
    con.execute(f"DROP TABLE IF EXISTS {q(unique_keys_table)}")
    n_eval = count_rows(con, eval_out)
    n_dedup = count_rows(con, dedup_out)
    print(
        f"  dedup done: {n_eval:,} rows -> {n_dedup:,} unique-key rows "
        f"({n_eval - n_dedup:,} dropped)"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--groups", nargs="+", choices=list(GROUPS), default=list(GROUPS),
        metavar="GROUP",
        help=f"Which score groups to process (default: {' '.join(GROUPS)}).",
    )
    parser.add_argument(
        "--scores-dir", type=Path, default=SCORES_DIR,
        help=f"Score tables root (default: {SCORES_DIR}).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help=f"Output root for the *_eval tables (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--variant-evals", "--variant_evals", dest="variant_evals",
        type=Path, default=VARIANT_EVALS,
        help=f"Variant-level eval table (default: {VARIANT_EVALS}).",
    )
    parser.add_argument(
        "--ensg-evals", "--ensg_evals", dest="ensg_evals",
        type=Path, default=ENSG_EVALS,
        help=f"Gene-level eval table (default: {ENSG_EVALS}).",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Rebuild outputs that already exist (default: skip them).",
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
        help="Directory for DuckDB spill files (default: <output dir>/.duckdb_spill).",
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
        help="Print the plan (inputs, joins, outputs, SQL) and exit.",
    )
    args = parser.parse_args()

    jobs = build_jobs(
        args.groups, args.scores_dir, args.output_dir, args.variant_evals, args.ensg_evals
    )
    if not jobs:
        sys.exit(
            f"ERROR: no score tables found for groups {args.groups} under "
            f"{args.scores_dir}."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or (args.output_dir / ".duckdb_spill")
    temp_dir.mkdir(parents=True, exist_ok=True)

    print(f"Groups: {args.groups}")
    print(f"Score tables to process: {len(jobs)}")
    print(f"Variant evals: {args.variant_evals}")
    print(f"Gene evals:    {args.ensg_evals}")
    print(f"Output root:   {args.output_dir}")
    print(f"Spill dir:     {temp_dir}")

    con = duckdb.connect()
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        for job in jobs:
            process_job(con, job, args)
    finally:
        con.close()
        shutil.rmtree(temp_dir, ignore_errors=True)

    print("\nAll done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
