#!/usr/bin/env python3
"""Build anchor-based pairwise score tables from the wide score parquet files.

A *pairwise* operation pairs one **anchor** score with every other
(*non-anchor*) score and restricts to the rows where **both** scores are
present -- the per-pair intersection. For an anchor ``polyphen`` and 16 other
scores this produces 16 small files (e.g. ``AM_polyphen_anchor.parquet``), each
carrying the key columns plus the anchor and non-anchor columns. The
intersection of every pair is different, so each file has its own row set.

Output files are named ``{non_anchor}_{anchor}_anchor.parquet`` for the raw
flavor and ``{non_anchor}_{anchor}_anchor_{pre|post}_percentile.parquet`` for
the percentile flavors.

Three flavors are produced, one per subdirectory inside the ``pairwise``
directory (default ``../data/processed_data/scores/pairwise``):

* ``pairwise_raw/``  -- the raw anchor/non-anchor score values on the pair's
  intersection. No percentile transform. Built from ``<base>_all.parquet``.

* ``pairwise_pre/``  -- percentiles computed **before** the pairwise subsetting:
  each column is already its *outer* (whole-distribution) percentile, then the
  pair is intersected. Built from ``<base>_outer_pre_percentile.parquet``
  (produced by ``create_percentile_score_tables.py``), so no percentile is
  recomputed here -- the precomputed outer percentiles are simply restricted to
  the pair's intersection.

* ``pairwise_post/`` -- percentiles computed **after** the pairwise subsetting:
  the pair is intersected first, then ``CUME_DIST`` is computed over that
  intersection for each column. These cannot reuse a precomputed table because
  each pair's intersection is a different subset, so they are computed from
  scratch from ``<base>_all.parquet``.

The ``pairwise_raw`` files keep the original score names (``polyphen``, ``AM``,
...). The percentile flavors rename each score column to record the transform:
``{score}_{anchor}_anchor_{pre|post}_percentile`` (e.g. ``AM_polyphen_anchor_
pre_percentile``), mirroring the column-renaming convention in
``create_percentile_score_tables.py``. The ``pre`` values are mathematically the
same as the corresponding non-pairwise outer percentiles; the rename is purely
for bookkeeping and downstream naming consistency.

Percentile definition (max tie-break)
-------------------------------------
``pairwise_post`` uses ``CUME_DIST`` over the pair's intersection: for ``N``
rows a value ``v`` maps to ``(count of values <= v) / N`` -- max tie-breaking,
so the maximum value is exactly ``1.0`` even when duplicated. ``pairwise_pre``
inherits the same convention from the precomputed outer-percentile table.

Datasets
--------
By default the script processes both ``variant_scores`` and ``ensg_scores``
(gene level). ``--variant`` / ``--ensg`` restrict it to one. Each requested
dataset's inputs must exist (the script errors otherwise), so while the
gene-level tables are not yet built, pass ``--variant``.

Anchor
------
``--anchor`` selects the anchor score; it defaults to ``polyphen`` when that
column is present in the dataset, otherwise it must be given explicitly.

Engine / memory strategy
------------------------
Same DuckDB out-of-core philosophy as the tables it consumes. Each pair is a
single streaming ``COPY (... ) TO`` over a ``read_parquet`` scan; ``pairwise_post``
additionally sorts the pair's intersection for ``CUME_DIST`` (DuckDB spills to
``--temp-dir`` under ``--memory-limit``). Only one pair is in flight at a time,
and existing output files are skipped (resumable) unless ``--overwrite`` is set.

Run from the ``src/`` directory::

    python create_pairwise_score_tables.py --variant
    python create_pairwise_score_tables.py --variant --anchor cadd
    python create_pairwise_score_tables.py --variant --modes raw post
    python create_pairwise_score_tables.py --variant --memory-limit 20GB
    python create_pairwise_score_tables.py --variant --dry-run
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

# Per-dataset key columns carried through unchanged; everything else is a score.
VARIANT_KEYS = ["chrom", "pos", "ref", "alt"]
ENSG_KEYS = ["ensg"]

DEFAULT_ANCHOR = "polyphen"
MODES = ("raw", "pre", "post")

# create_percentile_score_tables.py renames the score columns of the outer
# pre-percentile table to "<score>_pre_percentile". The pre flavor reads that
# table, so its input columns carry this suffix while the anchor/non-anchor are
# always referred to by their original score names.
PRE_INPUT_SUFFIX = "_pre_percentile"


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


def configure(
    con: duckdb.DuckDBPyConnection,
    memory_limit: str | None,
    threads: int | None,
    temp_dir: Path,
) -> None:
    """Apply memory / spill settings so large sorts stay out-of-core."""
    temp_dir.mkdir(parents=True, exist_ok=True)
    # temp_directory MUST be set for an in-memory connection to spill to disk.
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


def pair_filename(mode: str, anchor: str, na: str) -> str:
    """Output basename for one pair: ``{non_anchor}_{anchor}_anchor[...]``.

    ``raw``  -> ``{na}_{anchor}_anchor.parquet``
    ``pre``  -> ``{na}_{anchor}_anchor_pre_percentile.parquet``
    ``post`` -> ``{na}_{anchor}_anchor_post_percentile.parquet``
    """
    base = f"{na}_{anchor}_anchor"
    if mode == "raw":
        return f"{base}.parquet"
    return f"{base}_{mode}_percentile.parquet"


def pair_input_col(mode: str, score: str) -> str:
    """Input column name to read for a score in one mode.

    ``raw``/``post`` read the ``_all`` table, whose columns are the original
    score names. ``pre`` reads the outer pre-percentile table, whose columns are
    suffixed ``_pre_percentile``.
    """
    return f"{score}{PRE_INPUT_SUFFIX}" if mode == "pre" else score


def mode_score_names(
    con: duckdb.DuckDBPyConnection, reader: str, key_cols: list[str], mode: str
) -> list[str]:
    """Original score names available in a mode's input (suffix stripped for pre)."""
    cols = discover_score_columns(con, reader, key_cols)
    if mode != "pre":
        return cols
    originals = []
    for c in cols:
        if not c.endswith(PRE_INPUT_SUFFIX):
            sys.exit(
                f"ERROR: pre input column {c!r} does not end with "
                f"{PRE_INPUT_SUFFIX!r}; is this an outer_pre_percentile table "
                f"from create_percentile_score_tables.py?"
            )
        originals.append(c[: -len(PRE_INPUT_SUFFIX)])
    return originals


def pair_colname(mode: str, score: str, anchor: str) -> str:
    """Output column name for a score in one pair/mode.

    ``raw`` keeps the original score name; the percentile flavors rename each
    score column to ``{score}_{anchor}_anchor_{pre|post}_percentile`` so the
    column name records the pairwise transform (matching the file name and the
    column-renaming convention in ``create_percentile_score_tables.py``).
    """
    if mode == "raw":
        return score
    return f"{score}_{anchor}_anchor_{mode}_percentile"


def pair_query(
    mode: str, reader: str, key_cols: list[str], anchor: str, na: str
) -> str:
    """Build the SELECT for one anchor/non-anchor pair in the given mode.

    All modes carry the key columns plus the anchor and non-anchor columns; the
    values and the score-column names differ:

    * ``raw``  -- the raw values on the pair's intersection (original names).
    * ``pre``  -- the precomputed outer percentiles restricted to the pair's
      intersection (no recompute; the input columns already hold percentiles).
    * ``post`` -- ``CUME_DIST`` over the pair's intersection (recomputed because
      each pair's intersection is a different subset).

    For ``pre``/``post`` each score column is renamed via :func:`pair_colname`.
    The anchor/non-anchor are passed as original score names; the actual input
    column is resolved via :func:`pair_input_col` (the pre input is suffixed).
    """
    keys_sql = ", ".join(q(k) for k in key_cols)
    anchor_in = q(pair_input_col(mode, anchor))
    na_in = q(pair_input_col(mode, na))
    anchor_out = q(pair_colname(mode, anchor, anchor))
    na_out = q(pair_colname(mode, na, anchor))
    if mode in ("raw", "post"):
        # Raw inputs: a score is present iff not NULL and not NaN.
        where = (
            f"({anchor_in} IS NOT NULL AND NOT isnan({anchor_in})) AND "
            f"({na_in} IS NOT NULL AND NOT isnan({na_in}))"
        )
    else:  # pre: the input already holds percentiles (NULL iff the raw was).
        where = f"{anchor_in} IS NOT NULL AND {na_in} IS NOT NULL"

    if mode == "post":
        # CUME_DIST is computed in DOUBLE but stored as single-precision FLOAT to
        # roughly halve the on-disk footprint; percentiles in (0, 1] retain ample
        # precision at float32. (raw/pre carry their inputs through, which are
        # already FLOAT.)
        select = (
            f"{keys_sql}, "
            f"CAST(CUME_DIST() OVER (ORDER BY {anchor_in}) AS FLOAT) AS {anchor_out}, "
            f"CAST(CUME_DIST() OVER (ORDER BY {na_in}) AS FLOAT) AS {na_out}"
        )
    else:
        select = f"{keys_sql}, {anchor_in} AS {anchor_out}, {na_in} AS {na_out}"

    return f"SELECT {select} FROM {reader} WHERE {where}"


def resolve_anchor(explicit: str | None, score_cols: list[str], label: str) -> str:
    """Pick the anchor: the explicit one if given, else the default if present."""
    if explicit is not None:
        if explicit not in score_cols:
            sys.exit(
                f"ERROR: {label}: anchor {explicit!r} is not a score column. "
                f"Available scores: {score_cols}"
            )
        return explicit
    if DEFAULT_ANCHOR in score_cols:
        return DEFAULT_ANCHOR
    sys.exit(
        f"ERROR: {label}: no --anchor given and the default {DEFAULT_ANCHOR!r} "
        f"is not present. Pass --anchor with one of: {score_cols}"
    )


def dataset_inputs(label: str, args: argparse.Namespace) -> tuple[Path, Path]:
    """Return the (raw `_all`, `_outer_pre_percentile`) inputs for a dataset."""
    if label == "variant":
        return args.variant_score_path, args.variant_pre_path
    return args.ensg_score_path, args.ensg_pre_path


def process_dataset(
    con: duckdb.DuckDBPyConnection,
    label: str,
    key_cols: list[str],
    out_base_dir: Path,
    args: argparse.Namespace,
) -> None:
    """Build the requested pairwise flavors for one dataset (variant/ensg)."""
    raw_path, pre_path = dataset_inputs(label, args)
    # Each mode reads from its own input: raw/post from the `_all` table, pre
    # from the precomputed `_outer_pre_percentile` table.
    mode_input = {"raw": raw_path, "post": raw_path, "pre": pre_path}

    # Validate that every input a requested mode needs exists.
    for mode in args.modes:
        ip = mode_input[mode]
        if not ip.exists():
            hint = (
                "Have you built it with create_variant_scores_all_table.py?"
                if mode in ("raw", "post")
                else "Have you built it with create_percentile_score_tables.py?"
            )
            sys.exit(
                f"ERROR: {label} '{mode}' input not found: {ip}\n       {hint}"
            )

    # Resolve the anchor against a reference input (prefer raw, else pre). The
    # anchor is always an original score name (pre input columns are suffixed).
    ref_mode = "raw" if any(m in ("raw", "post") for m in args.modes) else "pre"
    ref_cols = mode_score_names(
        con, reader_sql(mode_input[ref_mode].resolve()), key_cols, ref_mode
    )
    anchor = resolve_anchor(args.anchor, ref_cols, label)

    print(f"\n=== {label} ===")
    print(f"Anchor: {anchor}")
    print(f"Keys carried through: {key_cols}")
    print(f"Output base dir: {out_base_dir}")

    for mode in args.modes:
        reader = reader_sql(mode_input[mode].resolve())
        score_cols = mode_score_names(con, reader, key_cols, mode)
        if anchor not in score_cols:
            sys.exit(
                f"ERROR: {label} '{mode}': anchor {anchor!r} not among the "
                f"scores in {mode_input[mode].name}: {score_cols}"
            )
        non_anchor = [c for c in score_cols if c != anchor]
        out_dir = out_base_dir / f"pairwise_{mode}"

        print(f"\n--- {label} / pairwise_{mode} ({len(non_anchor)} pairs) ---")
        print(f"Input: {mode_input[mode]}")
        print(f"Output dir: {out_dir}")
        for na in non_anchor:
            print(f"  - {pair_filename(mode, anchor, na)}")

        if args.dry_run:
            continue

        out_dir.mkdir(parents=True, exist_ok=True)
        for i, na in enumerate(non_anchor, start=1):
            out_path = out_dir / pair_filename(mode, anchor, na)
            if out_path.exists() and not args.overwrite:
                print(f"  [{i}/{len(non_anchor)}] {out_path.name} -- SKIPPED (exists)",
                      flush=True)
                continue
            query = pair_query(mode, reader, key_cols, anchor, na)
            print(f"  [{i}/{len(non_anchor)}] {anchor} x {na} ...", flush=True)
            con.execute(copy_sql(query, out_path, args.compression, args.row_group_size))
            n = con.execute(
                f"SELECT count(*) FROM read_parquet({sql_str(str(out_path))})"
            ).fetchone()[0]
            print(f"      wrote {out_path.name}: {n:,} rows", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--variant", action="store_true",
        help="Only produce pairwise tables for the variant scores "
             "(default: both variant and ensg).",
    )
    selection.add_argument(
        "--ensg", action="store_true",
        help="Only produce pairwise tables for the ensg (gene-level) scores "
             "(default: both).",
    )
    parser.add_argument(
        "--anchor", default=None,
        help=f"Anchor score column (default: {DEFAULT_ANCHOR!r} when present).",
    )
    parser.add_argument(
        "--modes", nargs="+", choices=list(MODES), default=list(MODES),
        metavar="MODE",
        help=f"Which flavors to produce (default: {' '.join(MODES)}).",
    )
    parser.add_argument(
        "--variant-score-path", "--variant_score_path", dest="variant_score_path",
        type=Path, default=SCORES_DIR / "variant_scores_all_outer.parquet",
        help="Variant `_all_outer` parquet (raw/post input). "
             f"Default: {SCORES_DIR / 'variant_scores_all_outer.parquet'}",
    )
    parser.add_argument(
        "--variant-pre-path", "--variant_pre_path", dest="variant_pre_path",
        type=Path, default=SCORES_DIR / "variant_scores_outer_pre_percentile.parquet",
        help="Variant outer-pre-percentile parquet (pre input). "
             f"Default: {SCORES_DIR / 'variant_scores_outer_pre_percentile.parquet'}",
    )
    parser.add_argument(
        "--ensg-score-path", "--ensg_score_path", dest="ensg_score_path",
        type=Path, default=SCORES_DIR / "ensg_scores_all.parquet",
        help="Ensg `_all` parquet (raw/post input). "
             f"Default: {SCORES_DIR / 'ensg_scores_all.parquet'}",
    )
    parser.add_argument(
        "--ensg-pre-path", "--ensg_pre_path", dest="ensg_pre_path",
        type=Path, default=SCORES_DIR / "ensg_scores_outer_pre_percentile.parquet",
        help="Ensg outer-pre-percentile parquet (pre input). "
             f"Default: {SCORES_DIR / 'ensg_scores_outer_pre_percentile.parquet'}",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=SCORES_DIR / "pairwise",
        help="Base directory holding the pairwise_raw/pre/post subdirs "
             f"(default: {SCORES_DIR / 'pairwise'}).",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Recompute and overwrite existing pair files (default: skip them).",
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
        help="Print the plan (datasets, anchor, pairs, outputs) and exit.",
    )
    args = parser.parse_args()

    # Selection: default is both; --variant / --ensg restrict to one.
    want_variant = args.variant or not args.ensg
    want_ensg = args.ensg or not args.variant

    datasets: list[tuple[str, list[str]]] = []
    if want_variant:
        datasets.append(("variant", VARIANT_KEYS))
    if want_ensg:
        datasets.append(("ensg", ENSG_KEYS))

    out_base_dir = args.output_dir
    out_base_dir.mkdir(parents=True, exist_ok=True)

    temp_dir = args.temp_dir or (out_base_dir / ".duckdb_spill")
    temp_dir.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        for label, key_cols in datasets:
            process_dataset(con, label, key_cols, out_base_dir, args)
    finally:
        con.close()
        shutil.rmtree(temp_dir, ignore_errors=True)

    print("\nAll done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
