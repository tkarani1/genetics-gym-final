#!/usr/bin/env python3
"""Derive percentile-transformed score tables from a wide ``*_scores_all.parquet``.

Given a wide score table (one row per key plus one FLOAT column per score, with
full-outer semantics so most cells are NULL) this script writes a set of sibling
Parquet files into the same directory. For an input named ``<base>_all_outer.parquet``
(e.g. ``variant_scores_all_outer.parquet`` -> ``<base>`` = ``variant_scores``):

* ``<base>_outer_pre_percentile.parquet``
  The original (outer) table with every score column replaced by its per-column
  percentile. NULL stays NULL (percentiles are computed over the non-NULL values
  of each column independently).

* ``<base>_inner_pre_percentile.parquet``
  The ``..._outer_pre_percentile`` table with every row that has any NULL score
  dropped. Because percentiles preserve NULLs, dropping NULL rows here is
  equivalent to an INNER JOIN across the score columns of the original outer
  table -- the percentile values are the *outer* percentiles, just restricted
  to the fully-populated rows.

* ``<base>_all_inner.parquet``
  The original table with every NA-score row dropped (an INNER JOIN / drop-NA of
  the input). This is both a standalone deliverable and the intermediate that
  the ``inner_post`` table below is built from.

* ``<base>_inner_post_percentile.parquet``
  Built from ``<base>_all_inner.parquet``: percentiles computed on the inner
  (fully-populated) rows. These differ from the ``inner_pre`` percentiles
  because the denominator/ranking is the inner set, not the outer set.

* ``<base>_percentile_thresholds.tsv`` (optional byproduct, on by default)
  For each score, the raw score value at each requested percentile threshold
  (default ``[0.85, 0.9, 0.95, 0.98, 0.99, 0.995]``), over that score's full
  non-NULL distribution (the same basis as ``outer_pre_percentile``). Disable
  with ``--no-tsv`` or override the thresholds with ``--tsv-thresholds``.

Datasets
--------
By default the script processes both ``variant_scores_all_outer.parquet`` and
``ensg_scores_all.parquet`` (gene level). ``--variant`` / ``--ensg`` restrict it
to one; ``--variant-score-path`` / ``--ensg-score-path`` override the input
locations. Each requested dataset's input must exist (the script errors if it
does not), so pass ``--variant`` while the gene-level table is not yet built.

Percentile definition (max tie-break)
-------------------------------------
For a column with ``N`` non-NULL values, a value ``v`` maps to

    pct(v) = (count of non-NULL values <= v) / N

i.e. ``cume_dist`` with **max** tie-breaking: every member of a tie group gets
the rank of the *last* element. This guarantees the maximum value maps to
exactly ``1.0`` even when it is heavily duplicated (an "average" tie-break
would pull a 4%-mass maximum down to ~0.98).

Engine / memory strategy
------------------------
This is a **DuckDB out-of-core** pipeline, like the table it consumes. The
percentile transform is applied **one column at a time into an on-disk DuckDB
table**: for each score column we build a small distinct-value -> percentile map
(``GROUP BY`` to distinct values, then a single ordered running ``SUM`` for the
max-tie cumulative fraction) and ``LEFT JOIN`` it back, re-materializing the
wide table each step and dropping the previous one. Only one column's ranking is
in flight at a time, so peak memory/disk is bounded. Runtime is intentionally
traded for that footprint.

Run from the ``src/`` directory::

    python create_percentile_score_tables.py
    python create_percentile_score_tables.py --variant --memory-limit 20GB
    python create_percentile_score_tables.py --tsv-thresholds 0.9 0.99 0.999
    python create_percentile_score_tables.py --no-tsv --dry-run
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import duckdb

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: write_tables_pipeline/src/...)
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # write_tables_pipeline/
SCORES_DIR = PROJECT_DIR / "data" / "processed_data" / "scores"
DEFAULT_VARIANT_INPUT = SCORES_DIR / "variant_scores_all_outer.parquet"
DEFAULT_ENSG_INPUT = SCORES_DIR / "ensg_scores_all.parquet"

# Per-dataset key columns carried through unchanged; everything else is a score.
VARIANT_KEYS = ["chrom", "pos", "ref", "alt"]
ENSG_KEYS = ["ensg"]

DEFAULT_TSV_THRESHOLDS = [0.85, 0.9, 0.95, 0.98, 0.99, 0.995]


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


def discover_score_columns(
    con: duckdb.DuckDBPyConnection, reader: str, key_cols: list[str]
) -> list[str]:
    """Return the score columns of the input (all columns except the keys).

    Order is preserved from the input schema. Errors if a key column is missing.
    """
    described = con.execute(f"DESCRIBE SELECT * FROM {reader}").fetchall()
    all_cols = [r[0] for r in described]
    missing = [k for k in key_cols if k not in all_cols]
    if missing:
        sys.exit(
            f"ERROR: input is missing expected key column(s) {missing}. "
            f"Found columns: {all_cols}"
        )
    score_cols = [c for c in all_cols if c not in key_cols]
    if not score_cols:
        sys.exit("ERROR: input has no score columns (only the key columns).")
    return score_cols


def notna(col: str) -> str:
    """SQL predicate: a score value is present (not NULL and not NaN)."""
    return f"({q(col)} IS NOT NULL AND NOT isnan({q(col)}))"


def percentile_map_sql(from_expr: str, col: str) -> str:
    """Subquery mapping each distinct non-NULL value of ``col`` to its percentile.

    ``GROUP BY`` collapses the column to its distinct values with counts (so the
    ordered running sum is over distinct values, not all rows). The running
    ``SUM(cnt) OVER (ORDER BY val)`` includes every peer of a tie group (RANGE
    framing), giving the **max**-tie cumulative count; dividing by the total
    non-NULL count ``SUM(cnt) OVER ()`` yields the percentile in ``(0, 1]``.
    """
    return (
        "SELECT val, "
        "CAST(SUM(cnt) OVER (ORDER BY val) AS DOUBLE) / SUM(cnt) OVER () AS pct "
        "FROM (SELECT "
        f"{q(col)} AS val, COUNT(*) AS cnt FROM {from_expr} "
        f"WHERE {notna(col)} GROUP BY {q(col)})"
    )


def build_percentile_table(
    con: duckdb.DuckDBPyConnection,
    base_from: str,
    score_cols: list[str],
    out_table: str,
    log_prefix: str = "",
) -> None:
    """Materialize ``out_table`` = ``base_from`` with each score col percentiled.

    Processes one column at a time, re-materializing the wide table each step
    and dropping the previous intermediate, so only one column's ranking is in
    flight at once. ``base_from`` is a FROM-clause expression (a quoted table
    name or a ``read_parquet(...)`` reader) and is never dropped. A NULL/NaN
    score does not match the percentile map (LEFT JOIN), so it is carried
    through as NULL.
    """
    prev_from = base_from
    prev_drop: str | None = None
    n = len(score_cols)
    for i, col in enumerate(score_cols):
        cur = out_table if i == n - 1 else f"{out_table}__s{i}"
        # The percentile is computed in DOUBLE (in percentile_map_sql) but stored
        # as single-precision FLOAT to roughly halve the on-disk footprint;
        # percentiles in (0, 1] retain ample precision at float32.
        stmt = (
            f"CREATE TABLE {q(cur)} AS "
            f"SELECT s.* REPLACE (CAST(m.pct AS FLOAT) AS {q(col)}) "
            f"FROM {prev_from} AS s "
            f"LEFT JOIN ({percentile_map_sql(prev_from, col)}) AS m "
            f"ON s.{q(col)} = m.val"
        )
        print(f"  {log_prefix}[{i + 1}/{n}] percentile of {col} ...", flush=True)
        con.execute(stmt)
        if prev_drop is not None:
            con.execute(f"DROP TABLE {q(prev_drop)}")
            # Reclaim the dropped table's blocks now: without an explicit
            # CHECKPOINT, DuckDB keeps freed pages allocated in the single build
            # DB file, so the file would grow ~N-fold across the N stages and
            # exhaust the disk. Checkpointing bounds it to ~2 stages.
            con.execute("CHECKPOINT")
        prev_from = q(cur)
        prev_drop = cur


def write_threshold_tsv(
    con: duckdb.DuckDBPyConnection,
    reader: str,
    score_cols: list[str],
    thresholds: list[float],
    tsv_path: Path,
    log_prefix: str = "",
) -> None:
    """Write a TSV of the raw score value at each percentile threshold.

    Rows are scores, columns are thresholds; each cell is the smallest raw value
    whose (max-tie) percentile is >= the threshold -- i.e. the threshold-th
    quantile over the score's full non-NULL distribution. ``NA`` is written for
    an all-NULL column.
    """
    # CAST to DOUBLE so the emitted value is a plain float regardless of the
    # source column's numeric type (e.g. DECIMAL), keeping the TSV clean.
    agg_exprs = ", ".join(
        f"min(CAST(val AS DOUBLE)) FILTER (WHERE pct >= {float(t)!r}) AS t{i}"
        for i, t in enumerate(thresholds)
    )
    rows: list[list[str]] = []
    for col in score_cols:
        result = con.execute(
            f"SELECT {agg_exprs} FROM ({percentile_map_sql(reader, col)})"
        ).fetchone()
        rows.append([col] + ["NA" if v is None else repr(v) for v in result])

    header = ["score"] + [repr(float(t)) for t in thresholds]
    with tsv_path.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print(f"  {log_prefix}wrote thresholds tsv -> {tsv_path.name}", flush=True)


def configure(
    con: duckdb.DuckDBPyConnection,
    memory_limit: str | None,
    threads: int | None,
    temp_dir: Path,
) -> None:
    """Apply memory / spill settings so large sorts/joins stay out-of-core."""
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET temp_directory = {sql_str(str(temp_dir))}")
    # Lower peak memory for streaming COPY by not preserving row order.
    con.execute("SET preserve_insertion_order = false")
    if memory_limit:
        con.execute(f"SET memory_limit = {sql_str(memory_limit)}")
    if threads:
        con.execute(f"SET threads = {int(threads)}")


def copy_sql(source_sql: str, out_path: Path, compression: str, row_group_size: int) -> str:
    """Build a ``COPY (...) TO '<path>' (FORMAT PARQUET, ...)`` statement."""
    return (
        f"COPY ({source_sql}) TO {sql_str(str(out_path))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(compression)}, "
        f"ROW_GROUP_SIZE {int(row_group_size)})"
    )


def renamed_projection(
    con: duckdb.DuckDBPyConnection, from_expr: str, score_cols: list[str], suffix: str
) -> str:
    """Projection that keeps key columns and suffixes each score column's name.

    Preserves the original column order (from the schema of ``from_expr``) and
    only renames the score columns (``<score>`` -> ``<score><suffix>``); key
    columns pass through unchanged. The percentile *values* are untouched -- this
    only relabels the output column names.
    """
    all_cols = [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {from_expr}").fetchall()]
    score_set = set(score_cols)
    parts = [
        f"{q(c)} AS {q(c + suffix)}" if c in score_set else q(c)
        for c in all_cols
    ]
    return ", ".join(parts)


def output_names(input_path: Path) -> dict[str, str]:
    """Derive the output basenames for an input ``<base>_all_outer.parquet``.

    The full-outer input is named ``<base>_all_outer`` (its inner sibling is
    ``<base>_all_inner``); a bare ``<base>_all`` is also accepted for backward
    compatibility. The ``_outer`` join-type marker is stripped before deriving
    the percentile-table base so the outputs keep their canonical names.
    """
    stem = input_path.stem  # e.g. "variant_scores_all_outer"
    core = stem[:-len("_outer")] if stem.endswith("_outer") else stem  # -> "variant_scores_all"
    base = core[:-4] if core.endswith("_all") else core  # -> "variant_scores"
    return {
        "outer_pre": f"{base}_outer_pre_percentile.parquet",
        "inner_pre": f"{base}_inner_pre_percentile.parquet",
        "all_inner": f"{core}_inner.parquet",
        "inner_post": f"{base}_inner_post_percentile.parquet",
        "tsv": f"{base}_percentile_thresholds.tsv",
    }


def process_scores_table(
    con: duckdb.DuckDBPyConnection,
    label: str,
    input_path: Path,
    key_cols: list[str],
    out_dir: Path,
    args: argparse.Namespace,
) -> None:
    """Build all percentile outputs (and the optional TSV) for one scores table."""
    if not input_path.exists():
        sys.exit(
            f"ERROR: {label} input not found: {input_path}\n"
            f"       Have you built the {input_path.name} table?"
        )

    reader = reader_sql(input_path.resolve())
    score_cols = discover_score_columns(con, reader, key_cols)
    names = output_names(input_path)
    paths = {k: out_dir / v for k, v in names.items()}

    print(f"\n=== {label} ===")
    print(f"Input: {input_path}")
    print(f"Keys carried through: {key_cols}")
    print(f"Score columns ({len(score_cols)}): {score_cols}")
    print("Outputs:")
    for key in ("outer_pre", "inner_pre", "all_inner", "inner_post"):
        print(f"  - {paths[key]}")
    if args.tsv:
        print(f"  - {paths['tsv']}  (thresholds {args.tsv_thresholds})")

    if args.dry_run:
        print("  (dry run -- nothing written)")
        return

    # Score columns in the percentile outputs are renamed to make the transform
    # explicit: pre-percentile tables get "_pre_percentile", the post table gets
    # "_post_percentile". Keys and the all_inner table keep their original names.
    pre_proj = renamed_projection(con, reader, score_cols, "_pre_percentile")
    post_proj = renamed_projection(con, reader, score_cols, "_post_percentile")

    # --- 1) outer_pre: percentile every column of the full outer table --------
    print("[1/4] Building outer pre-percentile table ...", flush=True)
    outer_tbl = f"{label}_outer_pre"
    build_percentile_table(con, reader, score_cols, outer_tbl, log_prefix=f"{label} outer ")
    print(f"  writing {names['outer_pre']} ...", flush=True)
    con.execute(
        copy_sql(
            f"SELECT {pre_proj} FROM {q(outer_tbl)}",
            paths["outer_pre"], args.compression, args.row_group_size,
        )
    )

    # --- 2) inner_pre: drop NULL-score rows from outer_pre (== INNER JOIN) -----
    print(f"[2/4] Writing {names['inner_pre']} (drop-NA of outer_pre) ...", flush=True)
    inner_pre_filter = " AND ".join(f"{q(c)} IS NOT NULL" for c in score_cols)
    con.execute(
        copy_sql(
            f"SELECT {pre_proj} FROM {q(outer_tbl)} WHERE {inner_pre_filter}",
            paths["inner_pre"], args.compression, args.row_group_size,
        )
    )
    con.execute(f"DROP TABLE {q(outer_tbl)}")
    con.execute("CHECKPOINT")  # reclaim the outer table's blocks before continuing

    # --- 3) all_inner: INNER JOIN / drop-NA of the ORIGINAL table --------------
    # Written as a standalone deliverable AND used as the intermediate for the
    # inner_post percentiles below.
    print(f"[3/4] Writing {names['all_inner']} (drop-NA of original) ...", flush=True)
    drop_na = " AND ".join(notna(c) for c in score_cols)
    con.execute(
        copy_sql(
            f"SELECT * FROM {reader} WHERE {drop_na}",
            paths["all_inner"], args.compression, args.row_group_size,
        )
    )

    # --- 4) inner_post: percentile the all_inner table -------------------------
    print("[4/4] Building inner post-percentile table (from all_inner) ...", flush=True)
    inner_reader = reader_sql(paths["all_inner"].resolve())
    inner_tbl = f"{label}_inner_post"
    build_percentile_table(con, inner_reader, score_cols, inner_tbl, log_prefix=f"{label} inner ")
    print(f"  writing {names['inner_post']} ...", flush=True)
    con.execute(
        copy_sql(
            f"SELECT {post_proj} FROM {q(inner_tbl)}",
            paths["inner_post"], args.compression, args.row_group_size,
        )
    )
    con.execute(f"DROP TABLE {q(inner_tbl)}")
    con.execute("CHECKPOINT")

    # --- optional TSV byproduct (outer / full-column basis) --------------------
    if args.tsv:
        write_threshold_tsv(
            con, reader, score_cols, args.tsv_thresholds, paths["tsv"],
            log_prefix=f"{label} ",
        )

    # --- summary ---------------------------------------------------------------
    def count(path: Path) -> int:
        return con.execute(
            f"SELECT count(*) FROM read_parquet({sql_str(str(path))})"
        ).fetchone()[0]

    n_outer = count(paths["outer_pre"])
    n_inner_pre = count(paths["inner_pre"])
    n_all_inner = count(paths["all_inner"])
    n_inner_post = count(paths["inner_post"])
    print(f"\n{label} done:")
    print(f"  {names['outer_pre']:<46} {n_outer:>14,} rows")
    print(f"  {names['inner_pre']:<46} {n_inner_pre:>14,} rows")
    print(f"  {names['all_inner']:<46} {n_all_inner:>14,} rows")
    print(f"  {names['inner_post']:<46} {n_inner_post:>14,} rows")
    if not (n_inner_pre == n_all_inner == n_inner_post):
        print(
            f"\nWARNING: inner row counts disagree (inner_pre={n_inner_pre:,}, "
            f"all_inner={n_all_inner:,}, inner_post={n_inner_post:,}); they "
            f"should all be the INNER set. Investigate NaN/NULL handling."
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--variant", action="store_true",
        help="Only produce percentile tables for the variant scores table "
             "(default: both variant and ensg).",
    )
    selection.add_argument(
        "--ensg", action="store_true",
        help="Only produce percentile tables for the ensg (gene-level) scores "
             "table (default: both).",
    )
    parser.add_argument(
        "--variant-score-path", "--variant_score_path", dest="variant_score_path",
        type=Path, default=DEFAULT_VARIANT_INPUT,
        help=f"Custom variant *_scores_all.parquet path (default: {DEFAULT_VARIANT_INPUT}).",
    )
    parser.add_argument(
        "--ensg-score-path", "--ensg_score_path", dest="ensg_score_path",
        type=Path, default=DEFAULT_ENSG_INPUT,
        help=f"Custom ensg *_scores_all.parquet path (default: {DEFAULT_ENSG_INPUT}).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="Directory for the output files (default: same directory as each input).",
    )
    parser.add_argument(
        "--tsv-thresholds", nargs="+", type=float, default=DEFAULT_TSV_THRESHOLDS,
        metavar="P",
        help="Percentile thresholds (in (0, 1]) for the raw-value TSV byproduct "
             f"(default: {DEFAULT_TSV_THRESHOLDS}).",
    )
    parser.add_argument(
        "--no-tsv", action="store_true",
        help="Do not produce the per-score raw-value-at-threshold TSV byproduct.",
    )
    parser.add_argument(
        "--memory-limit", default=None,
        help="DuckDB memory limit, e.g. '20GB' (default: DuckDB's ~80%% of RAM). "
             "Lower it to force earlier spilling on constrained machines.",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="DuckDB worker threads (default: DuckDB auto).",
    )
    parser.add_argument(
        "--temp-dir", type=Path, default=None,
        help="Directory for the on-disk build DB + DuckDB spill files "
             "(default: <output dir>/.duckdb_build). Needs free space.",
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
        help="Print the plan (datasets, columns, outputs) and exit without writing.",
    )
    args = parser.parse_args()

    args.tsv = not args.no_tsv
    bad = [t for t in args.tsv_thresholds if not (0.0 < t <= 1.0)]
    if bad:
        sys.exit(f"ERROR: --tsv-thresholds must be in (0, 1]; got out-of-range {bad}.")

    # Selection: default is both; --variant / --ensg restrict to one.
    want_variant = args.variant or not args.ensg
    want_ensg = args.ensg or not args.variant

    datasets: list[tuple[str, Path, list[str]]] = []
    if want_variant:
        datasets.append(("variant", args.variant_score_path, VARIANT_KEYS))
    if want_ensg:
        datasets.append(("ensg", args.ensg_score_path, ENSG_KEYS))

    out_dir = args.output_dir
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)

    temp_dir = args.temp_dir or ((out_dir or SCORES_DIR) / ".duckdb_build")
    temp_dir.mkdir(parents=True, exist_ok=True)
    db_path = temp_dir / "build.duckdb"

    con = duckdb.connect(str(db_path))
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        for label, input_path, key_cols in datasets:
            dataset_out_dir = out_dir or input_path.parent
            dataset_out_dir.mkdir(parents=True, exist_ok=True)
            process_scores_table(con, label, input_path, key_cols, dataset_out_dir, args)
    finally:
        con.close()
        shutil.rmtree(temp_dir, ignore_errors=True)

    print("\nAll done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
