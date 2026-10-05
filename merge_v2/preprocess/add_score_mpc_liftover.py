#!/usr/bin/env python3
"""Fold ``mpc_liftover`` into the local ``variant_scores_all_outer.parquet``
without a full ``create_variant_scores_all_table.py`` rebuild.

``create_variant_scores_all_table.py`` (the real 21-source merge step) is
currently *disabled* in ``config.json``: the local 21-col
``variant_scores_all_outer.parquet`` (4 keys + 17 legacy scores, msa_pairformer
excluded) was instead produced by ``preprocess/shortcut_exclude_msa_pairformer.py``,
a bypass that reuses an already-built table. This script is the analogous
bypass for *adding* the new gnomAD v2.1.1 MPC-liftover score: rather than
re-running the full 21-source merge (tens of minutes, ~74GB of spill per the
merge script's own docstring), it reuses that merge script's own
``analyze_source``/``dedup_source_sql`` helpers (imported, not reimplemented)
to run just **one** additional ``FULL JOIN`` of the new source onto the
existing outer table -- the same per-step shape as one iteration of
``create_variant_scores_all_table.py``'s ``build_statements`` stage loop.

The new source's manifest entry is read directly from
``data_config/score_input_data.json`` (single source of truth for the
``{src_col: out_col}`` mapping), matched by output column name.

Safety model (mirrors ``shortcut_exclude_msa_pairformer.py``): writes a
sidecar first (``variant_scores_all_outer.22col.parquet``); only with
``--commit`` does it back up the current file
(``variant_scores_all_outer.21col.bak.parquet``) and atomically rename the
sidecar into place.

Usage::

    python add_score_mpc_liftover.py --dry-run
    python add_score_mpc_liftover.py                 # writes the sidecar only
    python add_score_mpc_liftover.py --commit         # swaps it into place
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import duckdb

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent  # merge_v2/
SRC_DIR = PROJECT_DIR / "src"
sys.path.insert(0, str(SRC_DIR))

import create_variant_scores_all_table as cvsat  # noqa: E402  (reused, not reimplemented)

SCORES_DIR_DEFAULT = PROJECT_DIR / "data" / "processed_data" / "scores"
MANIFEST = PROJECT_DIR / "data_config" / "score_input_data.json"

NEW_OUTPUT_COL = "mpc_liftover"
JOIN_KEYS = cvsat.JOIN_KEYS  # ["chrom", "pos", "ref", "alt"]

COMPRESSION = "zstd"
ROW_GROUP_SIZE = 512_000


def q(identifier: str) -> str:
    return cvsat.q(identifier)


def sql_str(value: str) -> str:
    return cvsat.sql_str(value)


def find_manifest_entry(output_col: str) -> dict:
    """Find the score_input_data.json entry whose score_fields maps to output_col."""
    import json

    with MANIFEST.open() as fh:
        data = json.load(fh)
    for entry in data.get("variant_level", []):
        for mapping in entry.get("score_fields", []):
            if output_col in mapping.values():
                return entry
    sys.exit(
        f"ERROR: no entry in {MANIFEST} has an output score name '{output_col}'. "
        f"Register it first (score_input_data.json)."
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--scores-dir", type=Path, default=SCORES_DIR_DEFAULT)
    parser.add_argument("--output-col", default=NEW_OUTPUT_COL)
    parser.add_argument("--memory-limit", default="16GB")
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--temp-dir", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--commit", action="store_true",
        help="After writing the sidecar, back up the current outer table and "
             "atomically swap the sidecar into its place.",
    )
    args = parser.parse_args()

    if args.dry_run and args.commit:
        sys.exit("ERROR: --dry-run and --commit are mutually exclusive.")

    scores_dir = args.scores_dir
    current = scores_dir / "variant_scores_all_outer.parquet"
    sidecar = scores_dir / "variant_scores_all_outer.22col.parquet"
    backup = scores_dir / "variant_scores_all_outer.21col.bak.parquet"

    if not current.exists():
        sys.exit(f"ERROR: current outer table not found: {current}")

    entry = find_manifest_entry(args.output_col)
    print(f"Manifest entry: {entry['file_path']} -> score_fields {entry['score_fields']}")

    con = duckdb.connect()
    temp_dir = Path(args.temp_dir) if args.temp_dir else (scores_dir / ".duckdb_add_score")
    cvsat.configure(con, args.memory_limit, args.threads, temp_dir)

    reader, key_casts, score_casts, snv_pred = cvsat.analyze_source(con, entry, compact=False)
    if list(score_casts.keys()) != [args.output_col]:
        sys.exit(
            f"ERROR: expected exactly one score column '{args.output_col}' from "
            f"this source, got {list(score_casts.keys())}."
        )
    new_source_sql = cvsat.dedup_source_sql(reader, key_casts, score_casts, snv_pred)

    current_reader = f"read_parquet({sql_str(str(current))})"
    using = "(" + ", ".join(q(k) for k in JOIN_KEYS) + ")"
    join_sql = (
        f"SELECT * FROM {current_reader}\n"
        f"FULL JOIN (\n{new_source_sql}\n) USING {using}"
    )
    copy_sql = (
        f"COPY ({join_sql}) TO {sql_str(str(sidecar))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(COMPRESSION)}, ROW_GROUP_SIZE {ROW_GROUP_SIZE})"
    )

    print(f"Current outer table: {current}")
    print(f"New source:          {entry['file_path']} (output col: {args.output_col})")
    print(f"Sidecar output:      {sidecar}")

    if args.dry_run:
        print("\n(dry run -- nothing written)")
        print(copy_sql + ";")
        return 0

    n_before = con.execute(
        f"SELECT count(*) FROM read_parquet({sql_str(str(current))})"
    ).fetchone()[0]
    cols_before = [
        r[0] for r in con.execute(
            f"DESCRIBE SELECT * FROM read_parquet({sql_str(str(current))})"
        ).fetchall()
    ]
    print(f"\nBefore: {n_before:,} rows x {len(cols_before)} cols")

    if sidecar.exists():
        print(f"sidecar already exists, skipping write: {sidecar}")
    else:
        t0 = time.time()
        print("Joining + writing sidecar ...", flush=True)
        con.execute(copy_sql)
        print(f"  done in {time.time() - t0:.1f}s")

    n_after = con.execute(
        f"SELECT count(*) FROM read_parquet({sql_str(str(sidecar))})"
    ).fetchone()[0]
    cols_after = [
        r[0] for r in con.execute(
            f"DESCRIBE SELECT * FROM read_parquet({sql_str(str(sidecar))})"
        ).fetchall()
    ]
    n_null_new = con.execute(
        f"SELECT count(*) FROM read_parquet({sql_str(str(sidecar))}) "
        f"WHERE {q(args.output_col)} IS NULL"
    ).fetchone()[0]
    n_new_only = n_after - n_before  # rows whose key wasn't already in `current`

    print(f"After:  {n_after:,} rows x {len(cols_after)} cols")
    print(f"  new rows contributed only by the new source: {n_new_only:,}")
    print(f"  '{args.output_col}' non-null: {n_after - n_null_new:,} ({(n_after - n_null_new) / n_after * 100:.2f}%)")
    print(f"  columns added: {[c for c in cols_after if c not in cols_before]}")

    con.close()
    import shutil
    shutil.rmtree(temp_dir, ignore_errors=True)

    if not args.commit:
        print("\nSidecar written. Inspect it, then rerun with --commit to swap it into place.")
        return 0

    if backup.exists():
        sys.exit(f"ERROR: refusing to commit; backup path already exists: {backup}")
    current.rename(backup)
    sidecar.rename(current)
    print(f"\nCommitted: {current.name} -> {backup.name}; sidecar promoted to {current.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
