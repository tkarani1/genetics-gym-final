#!/usr/bin/env python3
"""Attach gene (``ensg``) labels to the non-percentile variant score tables.

This is the **first step** of the gene-aggregation feature. For each
non-percentile wide variant score table (one row per ``(chrom, pos, ref, alt)``
plus one FLOAT column per score) it INNER-JOINs the table onto the variant ->
gene linker (``../data/raw_data/linker/linker_all.parquet``) and writes the
result to ``../data/processed_data/scores/gene_aggregated/{name}_ensg.parquet``.

The linker carries ``(chrom, pos, ref, alt, ensg)``. Because a single variant
can map to more than one gene, the INNER JOIN may **duplicate** rows of the
original score table (one copy per gene a variant belongs to). That is expected
and fine: downstream steps aggregate per ``ensg``, and it is acceptable for a
variant to contribute to two genes. Variants with no gene in the linker are
dropped (INNER JOIN semantics), which is the intended behavior for a
per-gene aggregation.

The output schema is the variant key, the ``ensg`` column, then every score
column from the input, in that order::

    chrom, pos, ref, alt, ensg, <score_1>, <score_2>, ...

Inputs (non-percentile variant score tables)
--------------------------------------------
By default both non-percentile variant score tables produced by the upstream
scripts are processed:

* ``variant_scores_all_outer.parquet``  (full outer table)
* ``variant_scores_all_inner.parquet``  (inner / drop-NA table)

The percentile tables (``*_pre_percentile`` / ``*_post_percentile``) are
intentionally **not** included here. Pass ``--input`` to process a specific
table (or set of tables) instead of the defaults.

Per-gene statistics (``--gene-agg-stats``)
------------------------------------------
With ``--gene-agg-stats``, after each ``{name}_ensg.parquet`` is written the
script computes per-``ensg`` statistics for every score column and writes a
second file ``{name}_ensg_stats.parquet`` with the original columns plus the
statistic columns appended (named ``{score}_{statistic}``). The table still has
one row per ``(variant, gene)``.

For each score ``c`` and gene ``g`` (only the non-NULL/non-NaN values of ``c``
within ``g`` contribute to ``g``'s statistics):

* ``c_mean``    -- mean of the gene's values (broadcast to every row of the gene).
* ``c_median``  -- median (continuous / interpolated; broadcast).
* ``c_max``     -- maximum (broadcast).
* ``c_p90`` / ``c_p95`` -- the 90th / 95th percentile score, using the **same
  discrete max-tie convention** as ``create_percentile_score_tables.py``: the
  smallest value ``v`` whose ``(count of values <= v) / N >= p`` (broadcast).
* ``c_zscore`` -- the **per-row** standard score ``(x - mean_g) / std_g`` where
  ``std_g`` is the **population** standard deviation of the gene's values. NULL
  when the gene has zero variance (``std_g = 0``, e.g. a single variant) or when
  this row's value is missing.
* ``c_modified_zscore`` -- the **per-row** modified z-score
  ``0.6745 * (x - median_g) / MAD_g`` where ``MAD_g = median(|x_i - median_g|)``
  over the gene's values (the Iglewicz-Hoaglin / statology definition). NULL when
  ``MAD_g = 0`` or when this row's value is missing.

``std`` and ``MAD`` are intermediates for the z-scores and are not emitted as
columns. The mean/median/max/p90/p95 columns are gene-level properties, so they
are populated for every row of a gene (including rows whose own score is NULL).

All appended statistic columns are stored as single-precision ``FLOAT`` (the
arithmetic is done in DOUBLE and only the stored value is narrowed). These
~7-per-score columns are otherwise high-entropy and compress poorly, so this
roughly halves the output size; the raw score columns are themselves already
``FLOAT`` (carried through from the wide score table).

Engine / memory strategy
------------------------
Same DuckDB out-of-core philosophy as the tables it consumes. Each
gene-annotation output is a single streaming ``COPY (...) TO`` over a 2-table
``read_parquet`` INNER JOIN; DuckDB builds one hash table (the small 5-column
linker) and probes it with the score table, spilling to ``--temp-dir`` under
``--memory-limit`` as needed.

The statistics step favors memory/disk safety: per-gene statistics are computed
**one score column at a time** into a tiny per-gene dimension table (~one row per
gene), reading only the ``ensg`` + that one score column from the (columnar)
Parquet each pass and spilling holistic aggregates (median/quantile) to
``--temp-dir``. The final ``{name}_ensg_stats.parquet`` is then a single
streaming ``COPY`` of the ``_ensg`` table LEFT-JOINed to that tiny dimension
table (computing the per-row z-scores on the fly) -- no wide ~80M-row
intermediate is ever materialized.

Existing outputs are skipped (resumable) unless ``--overwrite`` is set.

Run from the ``src/`` directory::

    python create_variant_scores_gene_aggregation.py
    python create_variant_scores_gene_aggregation.py --input variant_scores_all_outer.parquet
    python create_variant_scores_gene_aggregation.py --gene-agg-stats
    python create_variant_scores_gene_aggregation.py --gene-agg-stats --memory-limit 20GB
    python create_variant_scores_gene_aggregation.py --dry-run
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import duckdb

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: merge_v2/src/...)
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # merge_v2/
SCORES_DIR = PROJECT_DIR / "data" / "processed_data" / "scores"
LINKER_PARQUET = PROJECT_DIR / "data" / "raw_data" / "linker" / "linker_all.parquet"
DEFAULT_OUTPUT_DIR = SCORES_DIR / "gene_aggregated"

# The non-percentile variant score tables this step gene-annotates by default.
DEFAULT_INPUTS = [
    SCORES_DIR / "variant_scores_all_outer.parquet",
    SCORES_DIR / "variant_scores_all_inner.parquet",
]

# Variant key shared by the score tables and the linker; everything else in a
# score table is a score column carried through unchanged.
VARIANT_KEYS = ["chrom", "pos", "ref", "alt"]
# Gene column contributed by the linker.
ENSG_COL = "ensg"

OUTPUT_SUFFIX = "_ensg"
STATS_SUFFIX = "_ensg_stats"

# Percentile-score statistics, computed with the discrete max-tie convention of
# create_percentile_score_tables.py (smallest value v with cume(v) >= p).
STAT_PERCENTILES = {"p90": 0.90, "p95": 0.95}
# Consistency constant for the modified z-score: 0.6745 is the 0.75 quantile of
# the standard normal, making MAD a consistent estimator of sigma.
MOD_Z_CONST = 0.6745


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


def score_columns(con: duckdb.DuckDBPyConnection, reader: str, input_path: Path) -> list[str]:
    """Return the score columns of an input (all columns except the variant key).

    Order is preserved from the input schema. Errors if a key column is missing.
    """
    all_cols = column_names(con, reader)
    missing = [k for k in VARIANT_KEYS if k not in all_cols]
    if missing:
        sys.exit(
            f"ERROR: {input_path.name} is missing expected key column(s) {missing}. "
            f"Found columns: {all_cols}"
        )
    scores = [c for c in all_cols if c not in VARIANT_KEYS]
    if ENSG_COL in all_cols:
        sys.exit(
            f"ERROR: {input_path.name} already has an '{ENSG_COL}' column; "
            f"this step expects a variant-level (non-gene-annotated) table."
        )
    return scores


def join_sql(score_reader: str, linker_reader: str, scores: list[str]) -> str:
    """Build the INNER-JOIN projection that gene-annotates a score table.

    Joins the score table (``s``) to the linker (``l``) on the variant key and
    projects ``key, ensg, <scores>``. The join can duplicate a score row when a
    variant maps to several genes (one output row per gene) -- this is intended.
    """
    on = " AND ".join(f"s.{q(k)} = l.{q(k)}" for k in VARIANT_KEYS)
    select_cols = [f"s.{q(k)}" for k in VARIANT_KEYS]
    select_cols.append(f"l.{q(ENSG_COL)} AS {q(ENSG_COL)}")
    select_cols += [f"s.{q(c)}" for c in scores]
    return (
        f"SELECT {', '.join(select_cols)} "
        f"FROM {score_reader} AS s "
        f"JOIN {linker_reader} AS l ON {on}"
    )


def copy_sql(source_sql: str, out_path: Path, compression: str, row_group_size: int) -> str:
    """Build a ``COPY (...) TO '<path>' (FORMAT PARQUET, ...)`` statement."""
    return (
        f"COPY ({source_sql}) TO {sql_str(str(out_path))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(compression)}, "
        f"ROW_GROUP_SIZE {int(row_group_size)})"
    )


def output_path(out_dir: Path, input_path: Path) -> Path:
    """Derive the gene-annotated output path for an input score table."""
    return out_dir / f"{input_path.stem}{OUTPUT_SUFFIX}.parquet"


def stats_output_path(out_dir: Path, input_path: Path) -> Path:
    """Derive the stats-augmented output path for an input score table."""
    return out_dir / f"{input_path.stem}{STATS_SUFFIX}.parquet"


def notna(col: str, alias: str | None = None) -> str:
    """SQL predicate: a score value is present (not NULL and not NaN)."""
    ref = f"{alias}.{q(col)}" if alias else q(col)
    return f"({ref} IS NOT NULL AND NOT isnan({ref}))"


def build_gene_stats(
    con: duckdb.DuckDBPyConnection,
    reader: str,
    scores: list[str],
    stats_table: str = "gene_stats",
) -> None:
    """Build a tiny per-gene dimension table of statistics, one score at a time.

    The result has one row per ``ensg`` and, for each score ``c``, the columns
    ``c_mean, c_median, c_max, c_p90, c_p95, c_std, c_mad`` (``std``/``mad`` are
    kept for the z-score formulas; they are not emitted in the final output).

    Each score is processed independently, reading only ``ensg`` + that one score
    column from the (columnar) Parquet, so peak memory stays bounded to a single
    score's per-gene aggregation (which DuckDB spills to ``--temp-dir``).
    """
    con.execute(
        f"CREATE OR REPLACE TABLE {q(stats_table)} AS "
        f"SELECT DISTINCT {q(ENSG_COL)} AS ensg FROM {reader}"
    )

    p90 = STAT_PERCENTILES["p90"]
    p95 = STAT_PERCENTILES["p95"]
    n = len(scores)
    for i, c in enumerate(scores):
        qc = q(c)
        # Base aggregates (mean, continuous median, max, population std).
        con.execute(
            "CREATE OR REPLACE TEMP TABLE _base AS "
            f"SELECT {q(ENSG_COL)} AS ensg, avg({qc}) AS mean, median({qc}) AS median, "
            f"max({qc}) AS max, stddev_pop({qc}) AS std "
            f"FROM {reader} WHERE {notna(c)} GROUP BY {q(ENSG_COL)}"
        )
        # Discrete max-tie percentiles per gene: smallest value v whose
        # cumulative fraction (count of values <= v) / N reaches the threshold.
        con.execute(
            "CREATE OR REPLACE TEMP TABLE _pct AS "
            "WITH dv AS ("
            f"SELECT {q(ENSG_COL)} AS ensg, {qc} AS val, count(*) AS cnt "
            f"FROM {reader} WHERE {notna(c)} GROUP BY {q(ENSG_COL)}, {qc}"
            "), cf AS ("
            "SELECT ensg, val, "
            "CAST(SUM(cnt) OVER (PARTITION BY ensg ORDER BY val) AS DOUBLE) "
            "/ SUM(cnt) OVER (PARTITION BY ensg) AS pct FROM dv) "
            "SELECT ensg, "
            f"min(val) FILTER (WHERE pct >= {p90!r}) AS p90, "
            f"min(val) FILTER (WHERE pct >= {p95!r}) AS p95 "
            "FROM cf GROUP BY ensg"
        )
        # MAD = median(|x - median_g|), using the gene's continuous median.
        con.execute(
            "CREATE OR REPLACE TEMP TABLE _mad AS "
            f"SELECT r.{q(ENSG_COL)} AS ensg, median(abs(r.{qc} - b.median)) AS mad "
            f"FROM {reader} AS r JOIN _base AS b ON r.{q(ENSG_COL)} = b.ensg "
            f"WHERE {notna(c, 'r')} GROUP BY r.{q(ENSG_COL)}"
        )
        # Merge this score's stats onto the dimension table (tiny, ~one row/gene).
        new_table = f"{stats_table}__new"
        con.execute(
            f"CREATE TABLE {q(new_table)} AS SELECT g.*, "
            f"b.mean AS {q(c + '_mean')}, b.median AS {q(c + '_median')}, "
            f"b.max AS {q(c + '_max')}, p.p90 AS {q(c + '_p90')}, "
            f"p.p95 AS {q(c + '_p95')}, b.std AS {q(c + '_std')}, "
            f"m.mad AS {q(c + '_mad')} "
            f"FROM {q(stats_table)} AS g "
            "LEFT JOIN _base AS b USING (ensg) "
            "LEFT JOIN _pct AS p USING (ensg) "
            "LEFT JOIN _mad AS m USING (ensg)"
        )
        con.execute(f"DROP TABLE {q(stats_table)}")
        con.execute(f"ALTER TABLE {q(new_table)} RENAME TO {q(stats_table)}")
        print(f"    [{i + 1}/{n}] per-gene stats for {c} ...", flush=True)

    for tmp in ("_base", "_pct", "_mad"):
        con.execute(f"DROP TABLE IF EXISTS {tmp}")


def stats_projection(scores: list[str]) -> str:
    """Projection of the appended statistic columns for the final stats COPY.

    ``t`` aliases the ``_ensg`` table (per-row scores) and ``g`` the per-gene
    dimension table. mean/median/max/p90/p95 are broadcast from ``g``; the
    z-score and modified z-score are computed per row, NULL where the gene has no
    spread (std/MAD 0 or NULL) or where this row's own score is missing.

    Every appended statistic column is stored as ``FLOAT`` (single precision):
    the arithmetic is done in DOUBLE and only the stored value is narrowed, which
    roughly halves the (otherwise high-entropy, poorly-compressing) footprint of
    these ~7-per-score columns. The raw score columns (carried through via
    ``t.*``) are themselves already ``FLOAT`` in the source table.
    """
    parts: list[str] = []
    for c in scores:
        tc = f"t.{q(c)}"
        g_mean = f"g.{q(c + '_mean')}"
        g_median = f"g.{q(c + '_median')}"
        g_std = f"g.{q(c + '_std')}"
        g_mad = f"g.{q(c + '_mad')}"
        parts.append(f"CAST({g_mean} AS FLOAT) AS {q(c + '_mean')}")
        parts.append(f"CAST({g_median} AS FLOAT) AS {q(c + '_median')}")
        parts.append(f"CAST(g.{q(c + '_max')} AS FLOAT) AS {q(c + '_max')}")
        parts.append(f"CAST(g.{q(c + '_p90')} AS FLOAT) AS {q(c + '_p90')}")
        parts.append(f"CAST(g.{q(c + '_p95')} AS FLOAT) AS {q(c + '_p95')}")
        parts.append(
            f"CAST(CASE WHEN {g_std} IS NULL OR {g_std} = 0 OR NOT {notna(c, 't')} "
            f"THEN NULL ELSE ({tc} - {g_mean}) / {g_std} END AS FLOAT) "
            f"AS {q(c + '_zscore')}"
        )
        parts.append(
            f"CAST(CASE WHEN {g_mad} IS NULL OR {g_mad} = 0 OR NOT {notna(c, 't')} "
            f"THEN NULL ELSE {MOD_Z_CONST!r} * ({tc} - {g_median}) / {g_mad} END "
            f"AS FLOAT) AS {q(c + '_modified_zscore')}"
        )
    return ",\n  ".join(parts)


def write_gene_stats(
    con: duckdb.DuckDBPyConnection,
    ensg_path: Path,
    scores: list[str],
    stats_out: Path,
    args: argparse.Namespace,
) -> None:
    """Compute per-gene stats from a ``_ensg.parquet`` and write ``_ensg_stats.parquet``."""
    if stats_out.exists() and not args.overwrite:
        print(f"  stats skip ({stats_out.name} exists; use --overwrite to rebuild)")
        return

    reader = reader_sql(ensg_path.resolve())
    print(f"  computing per-gene statistics for {len(scores)} score(s) ...", flush=True)
    build_gene_stats(con, reader, scores)

    source_sql = (
        f"SELECT t.*,\n  {stats_projection(scores)}\n"
        f"FROM {reader} AS t LEFT JOIN gene_stats AS g ON t.{q(ENSG_COL)} = g.ensg"
    )
    print(f"  writing {stats_out.name} ...", flush=True)
    con.execute(copy_sql(source_sql, stats_out, args.compression, args.row_group_size))
    con.execute("DROP TABLE IF EXISTS gene_stats")

    n_rows = con.execute(
        f"SELECT count(*) FROM read_parquet({sql_str(str(stats_out))})"
    ).fetchone()[0]
    n_cols = len(column_names(con, reader_sql(stats_out.resolve())))
    print(f"  stats done: {n_rows:,} rows x {n_cols} columns -> {stats_out.name}")


def configure(
    con: duckdb.DuckDBPyConnection,
    memory_limit: str | None,
    threads: int | None,
    temp_dir: Path,
) -> None:
    """Apply memory / spill settings so the large join stays out-of-core."""
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET temp_directory = {sql_str(str(temp_dir))}")
    # Lower peak memory for the streaming COPY by not preserving row order.
    con.execute("SET preserve_insertion_order = false")
    if memory_limit:
        con.execute(f"SET memory_limit = {sql_str(memory_limit)}")
    if threads:
        con.execute(f"SET threads = {int(threads)}")


def annotate_table(
    con: duckdb.DuckDBPyConnection,
    input_path: Path,
    linker_reader: str,
    out_dir: Path,
    args: argparse.Namespace,
) -> None:
    """INNER-JOIN one score table with the linker and write ``{name}_ensg.parquet``."""
    if not input_path.exists():
        sys.exit(
            f"ERROR: input not found: {input_path}\n"
            f"       Have you built {input_path.name} "
            f"(create_variant_scores_all_table.py / create_percentile_score_tables.py)?"
        )

    score_reader = reader_sql(input_path.resolve())
    scores = score_columns(con, score_reader, input_path)
    out_path = output_path(out_dir, input_path)
    stats_out = stats_output_path(out_dir, input_path)
    source_sql = join_sql(score_reader, linker_reader, scores)

    print(f"\n=== {input_path.name} ===")
    print(f"Input:  {input_path}")
    print(f"Linker: {LINKER_PARQUET}")
    print(f"Join key: {VARIANT_KEYS}  (INNER JOIN; rows may duplicate per gene)")
    print(f"Score columns ({len(scores)}): {scores}")
    print(f"Output: {out_path}")
    if args.gene_agg_stats:
        print(f"Stats:  {stats_out}  (+7 columns per score: "
              "mean, median, max, p90, p95, zscore, modified_zscore)")

    if args.dry_run:
        print("  (dry run -- nothing written)")
        print(copy_sql(source_sql, out_path, args.compression, args.row_group_size) + ";")
        if args.gene_agg_stats:
            print("  [stats] per-gene statistics are computed at run time from "
                  f"{out_path.name}; output -> {stats_out.name}")
        return

    if out_path.exists() and not args.overwrite:
        print("  skip (output exists; use --overwrite to rebuild)")
    else:
        print("  joining + writing parquet ...", flush=True)
        con.execute(copy_sql(source_sql, out_path, args.compression, args.row_group_size))
        n_in = con.execute(f"SELECT count(*) FROM {score_reader}").fetchone()[0]
        n_out = con.execute(
            f"SELECT count(*) FROM read_parquet({sql_str(str(out_path))})"
        ).fetchone()[0]
        print(f"  done: {n_in:,} variant rows -> {n_out:,} (variant, gene) rows")

    if args.gene_agg_stats:
        write_gene_stats(con, out_path, scores, stats_out, args)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", nargs="+", type=Path, default=None, metavar="PARQUET",
        help="Score table(s) to gene-annotate (default: the non-percentile "
             f"variant tables {[p.name for p in DEFAULT_INPUTS]}). A bare "
             "filename is resolved against the scores directory.",
    )
    parser.add_argument(
        "--linker", type=Path, default=LINKER_PARQUET,
        help=f"variant->gene linker parquet (default: {LINKER_PARQUET}).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for the *_ensg.parquet outputs (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--gene-agg-stats", "--gene_agg_stats", dest="gene_agg_stats",
        action="store_true",
        help="After each {name}_ensg.parquet, also compute per-gene statistics "
             "(mean, median, max, p90, p95, z-score, modified z-score) for every "
             "score and write {name}_ensg_stats.parquet with those columns appended.",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Rebuild outputs that already exist (default: skip them).",
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
             "(default: <output dir>/.duckdb_build). Must live on a volume with "
             "enough free space.",
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
        help="Print the plan (inputs, columns, outputs, SQL) and exit without writing.",
    )
    args = parser.parse_args()

    # Resolve inputs: bare filenames resolve against the scores directory.
    if args.input is None:
        inputs = list(DEFAULT_INPUTS)
    else:
        inputs = [p if p.parent != Path(".") else SCORES_DIR / p for p in args.input]

    if not args.linker.exists():
        sys.exit(
            f"ERROR: linker not found: {args.linker}\n"
            f"       Run: python download_source_data.py --linker-table all"
        )
    linker_reader = reader_sql(args.linker.resolve())

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # On-disk build DB so the per-gene stats intermediates never have to fit in
    # RAM; DuckDB also spills operators here under --memory-limit.
    temp_dir = args.temp_dir or (out_dir / ".duckdb_build")
    temp_dir.mkdir(parents=True, exist_ok=True)
    db_path = temp_dir / "build.duckdb"

    con = duckdb.connect(str(db_path))
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        for input_path in inputs:
            annotate_table(con, input_path, linker_reader, out_dir, args)
    finally:
        con.close()
        shutil.rmtree(temp_dir, ignore_errors=True)

    print("\nAll done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
