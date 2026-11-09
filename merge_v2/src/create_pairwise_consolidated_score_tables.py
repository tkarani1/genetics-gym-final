#!/usr/bin/env python3
"""Consolidate the per-pair pairwise score tables into one table per flavor.

``create_pairwise_score_tables.py`` writes one small file per anchor/non-anchor
pair, in three flavors (``pairwise_raw`` / ``pairwise_pre`` / ``pairwise_post``).
For 16 non-anchor scores that is 16 files per flavor. This step collapses each
flavor's per-pair files into a **single** wide score table whose row universe is
the outer join of those per-pair tables -- i.e. every key where the anchor is
present together with at least one non-anchor score. The consolidated tables are
the score-table inputs to the ``pairwise_consolidated`` group of
``create_analysis_tables.py``, so that a single analysis table can be built per
flavor instead of one per pair.

Self-contained: the only inputs are the per-pair files under the ``pairwise``
directory (default ``../data/processed_data/scores/pairwise``). The non-anchor
scores and the anchor are discovered from the file names / schema, so this step
works for any anchor and any set of pairs that ``create_pairwise_score_tables.py``
produced.

Column / collision handling
----------------------------
Within a flavor every per-pair file shares the same key columns and the same
*anchor* column name, and carries one unique *non-anchor* column.

* Non-anchor columns are unique across the per-pair files (each non-anchor
  appears in exactly one file), so they are carried through unchanged.
* The anchor column collides across all the files. How it is resolved depends on
  the flavor, because that is exactly what distinguishes the flavors:

  - ``raw`` / ``pre`` -- the anchor value for a given key is **identical** in
    every pair (the raw value, resp. the precomputed *outer* percentile), so the
    16 copies are coalesced into a **single** anchor column (original name kept).

  - ``post`` -- the anchor's percentile is ``CUME_DIST`` recomputed over **each
    pair's own intersection**, so the value differs per pair. Every copy is kept
    and renamed ``{anchor_col}_{anchor}_x_{non_anchor}`` (e.g.
    ``polyphen_polyphen_anchor_post_percentile_polyphen_x_AM``) so all are
    preserved and uniquely named. The non-anchor post columns are already unique
    and keep their names.

Engine / memory strategy
------------------------
Same DuckDB out-of-core philosophy as the upstream scripts. The per-pair files
are combined with a sequential full-outer join on the key columns, materializing
one stage table at a time (each stage drops the previous one) so peak memory is
bounded; DuckDB spills to ``--temp-dir`` under ``--memory-limit``. Existing
consolidated outputs are skipped (resumable) unless ``--overwrite`` is set.

Run from the ``src/`` directory::

    python create_pairwise_consolidated_score_tables.py --variant
    python create_pairwise_consolidated_score_tables.py --variant --modes post
    python create_pairwise_consolidated_score_tables.py --variant --anchor cadd
    python create_pairwise_consolidated_score_tables.py --variant --memory-limit 20GB
    python create_pairwise_consolidated_score_tables.py --variant --dry-run
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

# Per-dataset key columns carried through unchanged; everything else is a score.
VARIANT_KEYS = ["chrom", "pos", "ref", "alt"]
ENSG_KEYS = ["ensg"]

DEFAULT_ANCHOR = "polyphen"
MODES = ("raw", "pre", "post")

# raw/pre share one anchor value across all pairs (coalesced); post recomputes
# the anchor percentile per pair (every copy kept and renamed).
COALESCE_MODES = {"raw", "pre"}


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


def copy_sql(source_sql: str, out_path: Path, compression: str, row_group_size: int) -> str:
    """Build a ``COPY (...) TO '<path>' (FORMAT PARQUET, ...)`` statement."""
    return (
        f"COPY ({source_sql}) TO {sql_str(str(out_path))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(compression)}, "
        f"ROW_GROUP_SIZE {int(row_group_size)})"
    )


# ---------------------------------------------------------------------------
# Naming (mirrors create_pairwise_score_tables.py so this step can locate and
# read the per-pair files it produced).
# ---------------------------------------------------------------------------
def flavor_dirname(mode: str) -> str:
    return f"pairwise_{mode}"


def pair_filename(mode: str, anchor: str, na: str) -> str:
    """Per-pair file basename written by create_pairwise_score_tables.py."""
    base = f"{na}_{anchor}_anchor"
    return f"{base}.parquet" if mode == "raw" else f"{base}_{mode}_percentile.parquet"


def file_suffix(mode: str, anchor: str) -> str:
    """The constant filename suffix following the non-anchor score token."""
    base = f"_{anchor}_anchor"
    return f"{base}.parquet" if mode == "raw" else f"{base}_{mode}_percentile.parquet"


def anchor_colname(mode: str, anchor: str) -> str:
    """Anchor column name inside a per-pair file (matches pair_colname)."""
    return anchor if mode == "raw" else f"{anchor}_{anchor}_anchor_{mode}_percentile"


def na_colname(mode: str, anchor: str, na: str) -> str:
    """Non-anchor column name inside a per-pair file (matches pair_colname)."""
    return na if mode == "raw" else f"{na}_{anchor}_anchor_{mode}_percentile"


def post_anchor_colname(mode: str, anchor: str, na: str) -> str:
    """Per-pair renamed anchor column for the post flavor (collision break)."""
    return f"{anchor_colname(mode, anchor)}_{anchor}_x_{na}"


def consolidated_filename(label: str, mode: str) -> str:
    prefix = "" if label == "variant" else f"{label}_"
    return f"{prefix}pairwise_{mode}_consolidated.parquet"


def discover_non_anchors(flavor_dir: Path, mode: str, anchor: str) -> list[str]:
    """Find the non-anchor score tokens from the per-pair file names (sorted)."""
    suffix = file_suffix(mode, anchor)
    nas: list[str] = []
    for p in sorted(flavor_dir.glob(f"*{suffix}")):
        if p.name.startswith("."):
            continue
        token = p.name[: -len(suffix)]
        if token:
            nas.append(token)
    return nas


def column_names(con: duckdb.DuckDBPyConnection, reader: str) -> list[str]:
    return [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {reader}").fetchall()]


def pair_projection(
    flavor_dir: Path, mode: str, key_cols: list[str], anchor: str, na: str
) -> str:
    """SELECT for one per-pair file, ready to feed the consolidation join.

    Coalesce flavors keep the anchor column under its original name (it is
    coalesced across files later); the post flavor renames the anchor per pair
    so all copies survive the join.
    """
    path = flavor_dir / pair_filename(mode, anchor, na)
    keys_sql = ", ".join(q(k) for k in key_cols)
    anchor_in = anchor_colname(mode, anchor)
    na_in = na_colname(mode, anchor, na)
    if mode in COALESCE_MODES:
        anchor_sel = q(anchor_in)
    else:
        anchor_sel = f"{q(anchor_in)} AS {q(post_anchor_colname(mode, anchor, na))}"
    return (
        f"SELECT {keys_sql}, {anchor_sel}, {q(na_in)} "
        f"FROM {reader_sql(path.resolve())}"
    )


def consolidate_flavor(
    con: duckdb.DuckDBPyConnection,
    label: str,
    key_cols: list[str],
    flavor_dir: Path,
    mode: str,
    anchor: str,
    out_path: Path,
    args: argparse.Namespace,
) -> None:
    """Full-outer-join the per-pair files of one flavor into ``out_path``."""
    nas = discover_non_anchors(flavor_dir, mode, anchor)
    if not nas:
        sys.exit(
            f"ERROR: no per-pair files found in {flavor_dir} matching "
            f"'*{file_suffix(mode, anchor)}'. Has create_pairwise_score_tables.py "
            f"run for anchor {anchor!r}, mode {mode!r}?"
        )

    anchor_in = anchor_colname(mode, anchor)
    # Validate the anchor column is actually present in a sample file.
    sample = flavor_dir / pair_filename(mode, anchor, nas[0])
    sample_cols = column_names(con, reader_sql(sample.resolve()))
    if anchor_in not in sample_cols:
        sys.exit(
            f"ERROR: expected anchor column {anchor_in!r} not found in "
            f"{sample.name} (columns: {sample_cols}). Wrong --anchor?"
        )

    print(f"\n--- {label} / {flavor_dirname(mode)} ({len(nas)} pairs) ---")
    print(f"Input dir: {flavor_dir}")
    print(f"Anchor column: {anchor_in} "
          f"({'coalesced' if mode in COALESCE_MODES else 'renamed per pair'})")
    print(f"Output: {out_path}")

    if args.dry_run:
        print("  (dry run -- nothing written)")
        return

    keys_excl = ", ".join(q(k) for k in key_cols)
    coalesce = mode in COALESCE_MODES

    stage = "_cons_stage_0"
    con.execute(f"DROP TABLE IF EXISTS {stage}")
    con.execute(
        f"CREATE TEMP TABLE {stage} AS "
        f"{pair_projection(flavor_dir, mode, key_cols, anchor, nas[0])}"
    )

    for i, na in enumerate(nas[1:], start=1):
        cur = f"_cons_stage_{i}"
        proj = pair_projection(flavor_dir, mode, key_cols, anchor, na)
        on = " AND ".join(f"a.{q(k)} = b.{q(k)}" for k in key_cols)
        key_coal = ", ".join(f"COALESCE(a.{q(k)}, b.{q(k)}) AS {q(k)}" for k in key_cols)
        if coalesce:
            select = (
                f"{key_coal}, "
                f"COALESCE(a.{q(anchor_in)}, b.{q(anchor_in)}) AS {q(anchor_in)}, "
                f"a.* EXCLUDE ({keys_excl}, {q(anchor_in)}), "
                f"b.* EXCLUDE ({keys_excl}, {q(anchor_in)})"
            )
        else:
            select = (
                f"{key_coal}, a.* EXCLUDE ({keys_excl}), b.* EXCLUDE ({keys_excl})"
            )
        con.execute(f"DROP TABLE IF EXISTS {cur}")
        con.execute(
            f"CREATE TEMP TABLE {cur} AS SELECT {select} "
            f"FROM {stage} a FULL OUTER JOIN ({proj}) b ON {on}"
        )
        con.execute(f"DROP TABLE {stage}")
        stage = cur
        print(f"  joined + {na} ({i + 1}/{len(nas)})", flush=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    con.execute(copy_sql(f"SELECT * FROM {stage}", out_path, args.compression, args.row_group_size))
    con.execute(f"DROP TABLE {stage}")

    n = con.execute(
        f"SELECT count(*) FROM read_parquet({sql_str(str(out_path))})"
    ).fetchone()[0]
    ncols = len(column_names(con, reader_sql(out_path.resolve())))
    print(f"  wrote {out_path.name}: {n:,} rows, {ncols} cols", flush=True)


def process_dataset(
    con: duckdb.DuckDBPyConnection,
    label: str,
    key_cols: list[str],
    args: argparse.Namespace,
) -> None:
    """Consolidate the requested flavors for one dataset (variant/ensg)."""
    pairwise_dir = args.pairwise_dir
    out_dir = args.output_dir

    print(f"\n=== {label} ===")
    print(f"Anchor: {args.anchor}")
    print(f"Pairwise input dir: {pairwise_dir}")
    print(f"Consolidated output dir: {out_dir}")

    for mode in args.modes:
        flavor_dir = pairwise_dir / flavor_dirname(mode)
        if not flavor_dir.is_dir():
            sys.exit(
                f"ERROR: {label} '{mode}' input dir not found: {flavor_dir}\n"
                f"       Run create_pairwise_score_tables.py (it writes "
                f"pairwise_{mode}/) first."
            )
        out_path = out_dir / consolidated_filename(label, mode)
        if out_path.exists() and not args.overwrite:
            print(f"\n--- {label} / {flavor_dirname(mode)} -- SKIPPED "
                  f"({out_path.name} exists; use --overwrite) ---")
            continue
        consolidate_flavor(
            con, label, key_cols, flavor_dir, mode, args.anchor, out_path, args
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--variant", action="store_true",
        help="Only consolidate the variant pairwise tables (default: both).",
    )
    selection.add_argument(
        "--ensg", action="store_true",
        help="Only consolidate the ensg (gene-level) pairwise tables (default: both).",
    )
    parser.add_argument(
        "--anchor", default=DEFAULT_ANCHOR,
        help=f"Anchor score used by create_pairwise_score_tables.py "
             f"(default: {DEFAULT_ANCHOR!r}).",
    )
    parser.add_argument(
        "--modes", nargs="+", choices=list(MODES), default=list(MODES),
        metavar="MODE",
        help=f"Which flavors to consolidate (default: {' '.join(MODES)}).",
    )
    parser.add_argument(
        "--pairwise-dir", "--pairwise_dir", dest="pairwise_dir",
        type=Path, default=SCORES_DIR / "pairwise",
        help="Directory holding the pairwise_raw/pre/post subdirs "
             f"(default: {SCORES_DIR / 'pairwise'}).",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=SCORES_DIR / "pairwise_consolidated",
        help="Output directory for the consolidated tables "
             f"(default: {SCORES_DIR / 'pairwise_consolidated'}).",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Rebuild consolidated outputs that already exist (default: skip them).",
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
        help="Print the plan (datasets, flavors, inputs, outputs) and exit.",
    )
    args = parser.parse_args()

    want_variant = args.variant or not args.ensg
    want_ensg = args.ensg or not args.variant
    datasets: list[tuple[str, list[str]]] = []
    if want_variant:
        datasets.append(("variant", VARIANT_KEYS))
    if want_ensg:
        datasets.append(("ensg", ENSG_KEYS))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or (args.output_dir / ".duckdb_spill")
    temp_dir.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        for label, key_cols in datasets:
            process_dataset(con, label, key_cols, args)
    finally:
        con.close()
        shutil.rmtree(temp_dir, ignore_errors=True)

    print("\nAll done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
