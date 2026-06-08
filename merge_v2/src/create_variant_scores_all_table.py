#!/usr/bin/env python3
"""Merge all variant-level score fields into a single ``variant_scores_all_outer.parquet``.

Reads the manifest at
``../data/processed_data/scores/score_input_data.json`` and, for every entry
under ``variant_level``, projects out the requested score column(s) and merges
everything on the variant key ``(chrom, pos, ref, alt)`` into one wide Parquet
file, keeping every variant from every source (full-outer semantics).

Engine / memory strategy
------------------------
This is a **DuckDB out-of-core** pipeline. Each score table is ~78M rows and
the full set is tens of GB, so a fully in-memory merge is not viable.

The merge is built **one source at a time into an on-disk DuckDB table**:
each source is first deduped to one row per variant key (``max()`` per score,
ignoring NULLs), then ``FULL JOIN``-ed onto the growing wide table, which is
re-materialized each step (``CREATE TABLE stageN AS ...; DROP TABLE stageN-1``).

This avoids the two traps that exhausted the disk (~74 GB spill) at this scale:

* A single 16-way join pipelines every hash-table build simultaneously.
* A ``UNION ALL`` + ``GROUP BY`` inflates the input ~16x into a ~1.25B-row
  long-format intermediate before aggregating.

Here only one 2-table join runs at a time and the row count stays
~one-per-variant throughout, so peak memory/disk is bounded by a single join
plus the (one) growing wide table. ``max()`` dedup also collapses duplicate
keys within a source deterministically.

DuckDB is used because it:

* **Materializes intermediates on disk** (persistent build DB) and **spills**
  operators to ``--temp-dir`` (bounded by ``--memory-limit``), so nothing has
  to fit in RAM. The build DB + spill scratch are removed on completion.
* **Reads ``.tsv.bgz`` directly** via streaming gzip decompression -- no
  multi-GB decompressed copies are written to disk (an earlier pure-Polars
  attempt had to inflate the bgz TSVs to ~54GB of scratch before scanning).
* Streams the result straight to Parquet with ``COPY ... TO``.

Runtime is intentionally traded for a bounded memory/disk footprint.

Source formats
--------------
* Partitioned Parquet directories (Spark/Hail ``part-*.parquet`` + ``_SUCCESS``)
  -> read via a ``**/*.parquet`` glob (the ``_SUCCESS`` marker is ignored).
* Bgzipped TSVs (``*.tsv.bgz``) -> ``read_csv(..., compression='gzip')`` with
  quoting disabled, since columns such as ``alleles``/``values`` contain
  unescaped JSON; only the key + score columns are projected anyway.

The chromosome key is normalized across formats: Parquet sources expose
``chrom`` while the TSV sources expose ``chr`` -- both map to canonical
``chrom``.

Run from the ``src/`` directory::

    python create_variant_scores_all_table.py
    python create_variant_scores_all_table.py --memory-limit 8GB --threads 4
    python create_variant_scores_all_table.py --dry-run

Source selection (mutually exclusive; match by file_path, basename, or
substring)::

    python create_variant_scores_all_table.py --exclude coalesced_snp_misfit.tsv.bgz
    python create_variant_scores_all_table.py --include AM_AM_score.parquet sift.parquet
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import duckdb

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: merge_v2/src/...)
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # merge_v2/
INPUT_JSON = PROJECT_DIR / "data" / "processed_data" / "scores" / "score_input_data.json"
DEFAULT_OUTPUT = PROJECT_DIR / "data" / "processed_data" / "scores" / "variant_scores_all_outer.parquet"

# Canonical join key -> accepted source column aliases (first match wins).
KEY_SPEC: dict[str, list[str]] = {
    "chrom": ["chrom", "chr"],
    "pos": ["pos"],
    "ref": ["ref"],
    "alt": ["alt"],
}
JOIN_KEYS = list(KEY_SPEC.keys())

# Null sentinels for the TSV sources. (Score columns are additionally made
# null-safe by TRY_CAST below, which nulls any non-numeric value; this list
# mainly protects the key columns and matches the house convention in
# merge/table_io.py.)
CSV_NULL_VALUES = ["NA", "N/A", "N/a", "n/a", "na", "Na", "NaN", "nan", ""]


def q(identifier: str) -> str:
    """Quote a SQL identifier (DuckDB double-quotes; e.g. 'ref' is reserved)."""
    return '"' + identifier.replace('"', '""') + '"'


def sql_str(value: str) -> str:
    """Quote a SQL string literal."""
    return "'" + value.replace("'", "''") + "'"


def load_config() -> list[dict]:
    """Load the variant_level entries from the manifest."""
    if not INPUT_JSON.exists():
        sys.exit(f"ERROR: manifest not found: {INPUT_JSON}")
    with INPUT_JSON.open() as fh:
        data = json.load(fh)
    entries = data.get("variant_level")
    if not entries:
        sys.exit("ERROR: manifest has no 'variant_level' entries to merge.")
    return entries


def _entry_matches(entry: dict, token: str) -> bool:
    """Whether a manifest entry matches an --include/--exclude token.

    Matches on the full ``file_path``, on its basename, or as a substring of the
    ``file_path`` (trailing slashes ignored).
    """
    fp = entry["file_path"].rstrip("/")
    base = Path(fp).name
    tok = token.rstrip("/")
    return tok == fp or tok == base or tok in fp


def select_entries(
    entries: list[dict], include: list[str] | None, exclude: list[str] | None
) -> list[dict]:
    """Filter manifest entries by --include / --exclude (mutually exclusive)."""
    if include and exclude:
        sys.exit("ERROR: --include and --exclude are mutually exclusive.")
    tokens = include or exclude
    if not tokens:
        return entries

    unmatched = [t for t in tokens if not any(_entry_matches(e, t) for e in entries)]
    if unmatched:
        available = "\n  ".join(e["file_path"] for e in entries)
        sys.exit(
            "ERROR: these "
            f"{'--include' if include else '--exclude'} token(s) matched no "
            f"manifest entry: {unmatched}\nAvailable file_path values:\n  {available}"
        )

    if include:
        selected = [e for e in entries if any(_entry_matches(e, t) for t in tokens)]
    else:
        selected = [e for e in entries if not any(_entry_matches(e, t) for t in tokens)]
    if not selected:
        sys.exit("ERROR: the include/exclude selection left no sources to merge.")
    return selected


def assert_unique_outputs(entries: list[dict]) -> None:
    """Reject duplicate output score names -- they would collide in the join."""
    seen: dict[str, str] = {}
    for entry in entries:
        for mapping in entry.get("score_fields", []):
            for _src_col, out_col in mapping.items():
                if out_col in seen:
                    sys.exit(
                        f"ERROR: duplicate output score name '{out_col}' "
                        f"(from {entry['file_path']} and {seen[out_col]})."
                    )
                seen[out_col] = entry["file_path"]


def resolve_path(file_path: str) -> Path:
    """Resolve a manifest file_path (relative to the project root)."""
    p = Path(file_path)
    return p if p.is_absolute() else (PROJECT_DIR / p)


def is_parquet(file_path: str) -> bool:
    return ".parquet" in file_path


def reader_sql(file_path: str, resolved: Path) -> str:
    """Build the DuckDB table-function expression that reads a source file."""
    if is_parquet(file_path):
        pattern = str(resolved / "**" / "*.parquet") if resolved.is_dir() else str(resolved)
        return f"read_parquet({sql_str(pattern)})"

    compression = "gzip" if resolved.suffix in (".bgz", ".gz") else "auto"
    nullstr = "[" + ", ".join(sql_str(v) for v in CSV_NULL_VALUES) + "]"
    return (
        "read_csv("
        f"{sql_str(str(resolved))}, "
        "delim='\t', header=true, quote='', "
        f"nullstr={nullstr}, "
        f"compression={sql_str(compression)}, "
        "all_varchar=true)"
    )


def source_columns(con: duckdb.DuckDBPyConnection, reader: str) -> list[str]:
    """Return the column names exposed by a reader expression (no full scan)."""
    return [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {reader}").fetchall()]


def analyze_source(con: duckdb.DuckDBPyConnection, entry: dict) -> tuple[str, dict[str, str], dict[str, str]]:
    """Introspect one source and return its reader, key casts, and score casts.

    Returns
    -------
    reader : str
        The DuckDB table-function expression reading the file.
    key_casts : dict[canonical_key -> sql_expr]
        Cast expression that yields each canonical join key (resolving the
        ``chr``/``chrom`` alias and casting ``pos`` to BIGINT).
    score_casts : dict[output_name -> sql_expr]
        Cast expression (``TRY_CAST(... AS FLOAT)``) for each score this
        source contributes, keyed by the output column name. Scores are stored
        as single-precision ``FLOAT`` to roughly halve the on-disk footprint.
    """
    file_path = entry["file_path"]
    resolved = resolve_path(file_path)
    if not resolved.exists():
        sys.exit(
            f"ERROR: source not found: {resolved}\n"
            f"       (declared in manifest as '{file_path}'). "
            f"Have you run download_source_data.py?"
        )

    reader = reader_sql(file_path, resolved)
    available = set(source_columns(con, reader))

    key_casts: dict[str, str] = {}
    for canon, aliases in KEY_SPEC.items():
        src_col = next((a for a in aliases if a in available), None)
        if src_col is None:
            sys.exit(
                f"ERROR: {file_path}: no column for join key '{canon}' "
                f"(looked for {aliases}). Available: {sorted(available)}"
            )
        cast_type = "BIGINT" if canon == "pos" else "VARCHAR"
        key_casts[canon] = f"TRY_CAST({q(src_col)} AS {cast_type})"

    score_casts: dict[str, str] = {}
    for mapping in entry.get("score_fields", []):
        for src_col, out_col in mapping.items():
            if src_col not in available:
                sys.exit(
                    f"ERROR: {file_path}: score column '{src_col}' not found. "
                    f"Available: {sorted(available)}"
                )
            score_casts[out_col] = f"TRY_CAST({q(src_col)} AS FLOAT)"

    return reader, key_casts, score_casts


def dedup_source_sql(reader: str, key_casts: dict[str, str], score_casts: dict[str, str]) -> str:
    """A per-source subquery: keys + this source's scores, deduped on the key.

    Casts keys/scores, then collapses to one row per ``(chrom, pos, ref, alt)``
    via ``max()`` (ignores NULLs, deterministic). Deduping each source *before*
    joining is what keeps the row count bounded at ~one-per-variant and prevents
    a duplicate key from multiplying rows across the join chain.
    """
    cast_cols = [f"{key_casts[k]} AS {q(k)}" for k in JOIN_KEYS]
    cast_cols += [f"{expr} AS {q(name)}" for name, expr in score_casts.items()]
    inner = "SELECT " + ", ".join(cast_cols) + f" FROM {reader}"

    keys_sql = ", ".join(q(k) for k in JOIN_KEYS)
    agg_sql = ", ".join(f"max({q(name)}) AS {q(name)}" for name in score_casts)
    select_list = f"{keys_sql}, {agg_sql}" if agg_sql else keys_sql
    # Drop rows with an incomplete variant key (e.g. blank/malformed source
    # rows): a NULL key component cannot join and is not a usable variant.
    nonnull = " AND ".join(f"{q(k)} IS NOT NULL" for k in JOIN_KEYS)
    return f"SELECT {select_list} FROM (\n  {inner}\n) WHERE {nonnull} GROUP BY {keys_sql}"


def build_statements(
    sources: list[tuple[str, dict[str, str], dict[str, str]]],
    final_table: str,
) -> tuple[list[str], list[str]]:
    """Build the sequential, materialized full-outer-join statements.

    Instead of one giant 16-way join (which pipelines all hash-table builds at
    once and spills tens of GB) or a ``UNION ALL`` + ``GROUP BY`` (which inflates
    the input ~Nx into a long-format intermediate), this materializes the wide
    table one source at a time::

        CREATE TABLE stage0 AS <deduped source 0>;
        CREATE TABLE stage1 AS
            SELECT * FROM stage0 FULL JOIN <deduped source 1> USING (keys);
        DROP TABLE stage0;
        ... (repeat) ...
        ALTER TABLE stageN RENAME TO <final_table>;

    Only one 2-table join runs at a time and the row count stays ~one-per-variant
    throughout, so peak memory/disk is bounded by a single join plus the (one)
    growing wide table -- not by all sources at once.
    """
    score_names: list[str] = []
    for _reader, _keys, score_casts in sources:
        score_names.extend(score_casts.keys())

    using = "(" + ", ".join(q(k) for k in JOIN_KEYS) + ")"
    statements: list[str] = []

    reader0, keys0, scores0 = sources[0]
    statements.append(f"CREATE TABLE stage0 AS\n{dedup_source_sql(reader0, keys0, scores0)}")

    prev = "stage0"
    for i in range(1, len(sources)):
        reader, keys, scores = sources[i]
        cur = f"stage{i}"
        statements.append(
            f"CREATE TABLE {cur} AS\nSELECT * FROM {prev}\n"
            f"FULL JOIN (\n{dedup_source_sql(reader, keys, scores)}\n) USING {using}"
        )
        statements.append(f"DROP TABLE {prev}")
        prev = cur

    statements.append(f"ALTER TABLE {prev} RENAME TO {q(final_table)}")
    return statements, score_names


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
    # Lower peak memory for the streaming COPY by not preserving row order.
    con.execute("SET preserve_insertion_order = false")
    if memory_limit:
        con.execute(f"SET memory_limit = {sql_str(memory_limit)}")
    if threads:
        con.execute(f"SET threads = {int(threads)}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT,
        help=f"Output parquet path (default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--memory-limit", default=None,
        help="DuckDB memory limit, e.g. '8GB' (default: DuckDB's ~80%% of RAM). "
             "Lower it to force earlier spilling on constrained machines.",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="DuckDB worker threads (default: DuckDB auto).",
    )
    parser.add_argument(
        "--temp-dir", type=Path, default=None,
        help="Directory for DuckDB spill files (default: <output dir>/.duckdb_spill). "
             "Must live on a volume with enough free space.",
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
        help="Print the assembled SQL and exit without writing.",
    )
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--include", nargs="+", metavar="FILE",
        help="Only merge these manifest sources (match by file_path, basename, "
             "or substring). Mutually exclusive with --exclude.",
    )
    selection.add_argument(
        "--exclude", nargs="+", metavar="FILE",
        help="Merge every manifest source except these (same matching as "
             "--include). Mutually exclusive with --include.",
    )
    args = parser.parse_args()

    entries = load_config()
    entries = select_entries(entries, args.include, args.exclude)
    assert_unique_outputs(entries)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Scratch dir holds the persistent build DB + DuckDB spill files. The build
    # DB lives on disk so intermediates never have to fit in RAM.
    temp_dir = args.temp_dir or (args.output.parent / ".duckdb_build")
    temp_dir.mkdir(parents=True, exist_ok=True)
    db_path = temp_dir / "build.duckdb"

    final_table = "variant_scores_all"
    copy_sql = (
        f"COPY (SELECT * FROM {q(final_table)}) TO {sql_str(str(args.output))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(args.compression)}, "
        f"ROW_GROUP_SIZE {int(args.row_group_size)})"
    )

    con = duckdb.connect(str(db_path))
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        sources = [analyze_source(con, entry) for entry in entries]
        statements, score_names = build_statements(sources, final_table)

        print(f"Merging {len(entries)} source(s) -> {len(score_names)} score field(s):")
        for name in score_names:
            print(f"  - {name}")
        print(f"Merge: sequential materialized FULL OUTER JOIN on {JOIN_KEYS}")
        print(f"Build DB / spill dir: {temp_dir}")
        print(f"Output: {args.output}\n")

        if args.dry_run:
            for stmt in statements:
                print(stmt + ";\n")
            print(copy_sql + ";")
            return 0

        print("Running out-of-core merge (this may take a while)...")
        for i, stmt in enumerate(statements, start=1):
            print(f"  [step {i}/{len(statements)}] {stmt.splitlines()[0]} ...", flush=True)
            con.execute(stmt)

        print("  writing parquet ...", flush=True)
        con.execute(copy_sql)

        n_rows = con.execute(
            f"SELECT count(*) FROM read_parquet({sql_str(str(args.output))})"
        ).fetchone()[0]
        cols = [
            r[0]
            for r in con.execute(
                f"DESCRIBE SELECT * FROM read_parquet({sql_str(str(args.output))})"
            ).fetchall()
        ]
    finally:
        con.close()
        # Remove the on-disk build DB + spill scratch.
        shutil.rmtree(temp_dir, ignore_errors=True)

    print(f"\nDone. Wrote {n_rows:,} rows x {len(cols)} columns to {args.output}")
    print(f"Columns: {cols}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
