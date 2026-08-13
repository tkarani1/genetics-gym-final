#!/usr/bin/env python3
"""Merge all gene-level eval fields into a single ``ensg_evals_all.parquet``.

Gene-level counterpart to ``create_variant_eval_all_table.py``. Reads the
manifest at ``../data/processed_data/evals/evals_input_data.json`` and, for
every entry under ``gene_level``, projects out the requested eval field(s) and
merges everything on the gene key ``ensg`` into one wide Parquet file, keeping
every gene from every source (full-outer semantics).

Field typing
------------
Eval fields are a mix of two kinds, typed per field:

* ``is_pos``-style labels -> normalized to ``BOOLEAN`` (handles ``true``/``false``
  strings, native booleans, and ``0``/``1``). A field is treated as boolean if
  its output name starts with ``is_pos`` or its source column is BOOLEAN.
  (The clinvar / genebass burden gene labels.)
* ``n_case_*`` / ``n_ctrl_*`` (and any other numeric) -> single-precision
  ``FLOAT``. A real (rather than BIGINT) type stays parallel with the variant
  merge and avoids silently nulling any future non-integer gene-level metric;
  ``FLOAT`` rather than DOUBLE roughly halves the on-disk footprint.

Engine / memory strategy
------------------------
Identical to ``create_variant_eval_all_table.py``: the wide table is built
**one source at a time into an on-disk DuckDB table**. Each source is first
deduped to one row per gene key, then ``FULL JOIN``-ed onto the growing table,
which is re-materialized each step. Only one 2-table join runs at a time and
the row count stays ~one-per-gene, so peak memory/disk is bounded. (Gene-level
sources are small -- a handful of MB each, one row per gene -- so in practice
this is fast; the staged design is kept purely for parity with the variant
script.)

Source formats
--------------
* Partitioned Parquet directories -> read via a ``**/*.parquet`` glob.
* TSV / bgzipped TSV -> ``read_csv`` with quoting disabled and ``all_varchar``;
  values are cast explicitly. (All current gene-level sources are plain TSVs.)

Run from the ``src/`` directory::

    python create_ensg_eval_all_table.py
    python create_ensg_eval_all_table.py --memory-limit 20GB
    python create_ensg_eval_all_table.py --dry-run

Source selection (mutually exclusive; match by file_path, basename, or
substring)::

    # merge only the two case/control LoF count sources for dd and asd
    python create_ensg_eval_all_table.py --include dd-lof-vep asd-lof-vep
    # merge everything except the schema source
    python create_ensg_eval_all_table.py --exclude schema_ensg_eval
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
INPUT_JSON = PROJECT_DIR / "data" / "processed_data" / "evals" / "evals_input_data.json"
DEFAULT_OUTPUT = PROJECT_DIR / "data" / "processed_data" / "evals" / "ensg_evals_all.parquet"

# Canonical join key -> accepted source column aliases (first match wins). All
# current gene-level sources expose the Ensembl gene id directly as ``ensg``;
# the extra aliases are defensive for future sources.
KEY_SPEC: dict[str, list[str]] = {
    "ensg": ["ensg", "gene_id", "ensembl_gene_id", "gene"],
}
JOIN_KEYS = list(KEY_SPEC.keys())

# Null sentinels for the TSV sources. (Numeric fields are additionally made
# null-safe by TRY_CAST; this list mainly protects the key/boolean columns and
# matches the house convention in merge/table_io.py.)
CSV_NULL_VALUES = ["NA", "N/A", "N/a", "n/a", "na", "Na", "NaN", "nan", ""]

# Truthy / falsy spellings accepted when normalizing a label to BOOLEAN.
BOOL_TRUE = ["true", "t", "1", "yes", "y"]
BOOL_FALSE = ["false", "f", "0", "no", "n"]


def q(identifier: str) -> str:
    """Quote a SQL identifier (DuckDB double-quotes; e.g. 'ref' is reserved)."""
    return '"' + identifier.replace('"', '""') + '"'


def sql_str(value: str) -> str:
    """Quote a SQL string literal."""
    return "'" + value.replace("'", "''") + "'"


def load_config() -> list[dict]:
    """Load the gene_level entries from the manifest."""
    if not INPUT_JSON.exists():
        sys.exit(f"ERROR: manifest not found: {INPUT_JSON}")
    with INPUT_JSON.open() as fh:
        data = json.load(fh)
    entries = data.get("gene_level")
    if not entries:
        sys.exit("ERROR: manifest has no 'gene_level' entries to merge.")
    return entries


def _entry_matches(entry: dict, token: str) -> bool:
    """Whether a manifest entry matches an --include/--exclude token.

    Matches on the full ``file_path``, on its basename, or as a substring of the
    ``file_path`` (trailing slashes ignored). So ``dd-lof-vep-case-ctrl-counts-ensg``,
    ``dd-lof-vep-case-ctrl-counts-ensg.tsv``, ``dd-lof-vep`` and the full
    relative path all select the dd LoF count source.
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
    """Reject duplicate output eval names -- they would collide in the join."""
    seen: dict[str, str] = {}
    for entry in entries:
        for mapping in entry.get("eval_fields", []):
            for _src_col, out_col in mapping.items():
                if out_col in seen:
                    sys.exit(
                        f"ERROR: duplicate output eval name '{out_col}' "
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


def source_column_types(con: duckdb.DuckDBPyConnection, reader: str) -> dict[str, str]:
    """Return {column_name: column_type} for a reader expression (no full scan)."""
    return {
        r[0]: r[1]
        for r in con.execute(f"DESCRIBE SELECT * FROM {reader}").fetchall()
    }


def _is_bool_field(out_col: str, src_type: str) -> bool:
    return out_col.lower().startswith("is_pos") or src_type.upper().startswith("BOOL")


def _bool_cast(src_col: str) -> str:
    """Robustly normalize a label column to BOOLEAN across encodings."""
    as_text = f"lower(TRY_CAST({q(src_col)} AS VARCHAR))"
    true_in = ", ".join(sql_str(v) for v in BOOL_TRUE)
    false_in = ", ".join(sql_str(v) for v in BOOL_FALSE)
    return (
        f"CASE WHEN {as_text} IN ({true_in}) THEN TRUE "
        f"WHEN {as_text} IN ({false_in}) THEN FALSE END"
    )


def analyze_source(
    con: duckdb.DuckDBPyConnection, entry: dict
) -> tuple[str, dict[str, str], dict[str, tuple[str, str]]]:
    """Introspect one source: reader, key casts, and per-field (cast, agg).

    Returns
    -------
    reader : str
        DuckDB table-function expression reading the file.
    key_casts : dict[canonical_key -> sql_expr]
        Cast yielding each canonical key (resolves the ``ensg`` alias to VARCHAR).
    field_specs : dict[output_name -> (cast_expr, agg_func)]
        Per eval field: the cast expression and the dedup aggregate to use
        (``bool_or`` for boolean labels, ``max`` for numeric fields).
    """
    file_path = entry["file_path"]
    resolved = resolve_path(file_path)
    if not resolved.exists():
        sys.exit(
            f"ERROR: source not found: {resolved}\n"
            f"       (declared in manifest as '{file_path}'). "
            f"Have you run download_source_data.py --eval-tables gene_level?"
        )

    reader = reader_sql(file_path, resolved)
    types = source_column_types(con, reader)
    available = set(types)

    key_casts: dict[str, str] = {}
    for canon, aliases in KEY_SPEC.items():
        src_col = next((a for a in aliases if a in available), None)
        if src_col is None:
            sys.exit(
                f"ERROR: {file_path}: no column for join key '{canon}' "
                f"(looked for {aliases}). Available: {sorted(available)}"
            )
        key_casts[canon] = f"TRY_CAST({q(src_col)} AS VARCHAR)"

    field_specs: dict[str, tuple[str, str]] = {}
    for mapping in entry.get("eval_fields", []):
        for src_col, out_col in mapping.items():
            if src_col not in available:
                sys.exit(
                    f"ERROR: {file_path}: eval column '{src_col}' not found. "
                    f"Available: {sorted(available)}"
                )
            if _is_bool_field(out_col, types[src_col]):
                field_specs[out_col] = (_bool_cast(src_col), "bool_or")
            else:
                field_specs[out_col] = (f"TRY_CAST({q(src_col)} AS FLOAT)", "max")

    return reader, key_casts, field_specs


def dedup_source_sql(
    reader: str, key_casts: dict[str, str], field_specs: dict[str, tuple[str, str]]
) -> str:
    """Per-source subquery: keys + this source's fields, deduped to one row/key.

    Casts the key/fields, drops rows with a null gene key, then collapses to one
    row per ``ensg`` using each field's aggregate. Deduping each source *before*
    joining keeps the row count ~one-per-gene and stops a duplicated key from
    multiplying rows across the join chain.
    """
    cast_cols = [f"{key_casts[k]} AS {q(k)}" for k in JOIN_KEYS]
    cast_cols += [f"{cast} AS {q(name)}" for name, (cast, _agg) in field_specs.items()]
    inner = "SELECT " + ", ".join(cast_cols) + f" FROM {reader}"

    keys_sql = ", ".join(q(k) for k in JOIN_KEYS)
    agg_sql = ", ".join(f"{agg}({q(name)}) AS {q(name)}" for name, (_c, agg) in field_specs.items())
    select_list = f"{keys_sql}, {agg_sql}" if agg_sql else keys_sql
    nonnull = " AND ".join(f"{q(k)} IS NOT NULL" for k in JOIN_KEYS)
    return f"SELECT {select_list} FROM (\n  {inner}\n) WHERE {nonnull} GROUP BY {keys_sql}"


def build_statements(
    sources: list[tuple[str, dict[str, str], dict[str, tuple[str, str]]]],
    final_table: str,
) -> tuple[list[str], list[str]]:
    """Build the sequential, materialized full-outer-join statements.

    Materializes the wide table one source at a time::

        CREATE TABLE stage0 AS <deduped source 0>;
        CREATE TABLE stage1 AS
            SELECT * FROM stage0 FULL JOIN <deduped source 1> USING (ensg);
        DROP TABLE stage0;
        ... (repeat) ...
        ALTER TABLE stageN RENAME TO <final_table>;

    Only one 2-table join runs at a time, so peak memory/disk is bounded.
    """
    field_names: list[str] = []
    for _reader, _keys, field_specs in sources:
        field_names.extend(field_specs.keys())

    using = "(" + ", ".join(q(k) for k in JOIN_KEYS) + ")"
    statements: list[str] = []

    reader0, keys0, fields0 = sources[0]
    statements.append(f"CREATE TABLE stage0 AS\n{dedup_source_sql(reader0, keys0, fields0)}")

    prev = "stage0"
    for i in range(1, len(sources)):
        reader, keys, fields = sources[i]
        cur = f"stage{i}"
        statements.append(
            f"CREATE TABLE {cur} AS\nSELECT * FROM {prev}\n"
            f"FULL JOIN (\n{dedup_source_sql(reader, keys, fields)}\n) USING {using}"
        )
        statements.append(f"DROP TABLE {prev}")
        prev = cur

    statements.append(f"ALTER TABLE {prev} RENAME TO {q(final_table)}")
    return statements, field_names


def configure(
    con: duckdb.DuckDBPyConnection,
    memory_limit: str | None,
    threads: int | None,
    temp_dir: Path,
) -> None:
    """Apply memory / spill settings so large joins stay out-of-core."""
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"SET temp_directory = {sql_str(str(temp_dir))}")
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
        help="DuckDB memory limit, e.g. '20GB' (default: DuckDB's ~80%% of RAM).",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="DuckDB worker threads (default: DuckDB auto).",
    )
    parser.add_argument(
        "--temp-dir", type=Path, default=None,
        help="Directory for the on-disk build DB + spill (default: "
             "<output dir>/.duckdb_build). Must have enough free space.",
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
    temp_dir = args.temp_dir or (args.output.parent / ".duckdb_build")
    temp_dir.mkdir(parents=True, exist_ok=True)
    db_path = temp_dir / "build.duckdb"

    final_table = "ensg_evals_all"
    copy_sql = (
        f"COPY (SELECT * FROM {q(final_table)}) TO {sql_str(str(args.output))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(args.compression)}, "
        f"ROW_GROUP_SIZE {int(args.row_group_size)})"
    )

    con = duckdb.connect(str(db_path))
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        sources = [analyze_source(con, entry) for entry in entries]
        statements, field_names = build_statements(sources, final_table)

        print(f"Merging {len(entries)} source(s) -> {len(field_names)} eval field(s):")
        for _reader, _keys, field_specs in sources:
            for name, (_cast, agg) in field_specs.items():
                kind = "bool" if agg == "bool_or" else "float"
                print(f"  - {name} ({kind})")
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
        shutil.rmtree(temp_dir, ignore_errors=True)

    print(f"\nDone. Wrote {n_rows:,} rows x {len(cols)} columns to {args.output}")
    print(f"Columns: {cols}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
