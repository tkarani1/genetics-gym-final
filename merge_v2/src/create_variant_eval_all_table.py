#!/usr/bin/env python3
"""Merge all variant-level eval fields into a single ``variant_evals_all.parquet``.

Reads the manifest at
``../data_config/evals_input_data.json`` and, for every entry
under ``variant_level``, projects out the requested eval field(s) and merges
everything on the variant key ``(chrom, pos, ref, alt)`` into one wide Parquet
file, keeping every variant from every source (full-outer semantics).

Field typing
------------
Eval fields are a mix of two kinds, typed per field:

* ``is_pos``-style labels -> normalized to ``BOOLEAN`` (handles ``true``/``false``
  strings, native booleans, and ``0``/``1``). A field is treated as boolean if
  its output name starts with ``is_pos`` or its source column is BOOLEAN.
* ``observed_*`` / ``expected_*`` (and any other numeric) -> single-precision
  ``FLOAT`` (halves the on-disk footprint versus DOUBLE).

Engine / memory strategy
------------------------
Identical to ``create_variant_scores_all_table.py``: the wide table is built
**one source at a time into an on-disk DuckDB table**. Each source is first
deduped to one row per variant key, then ``FULL JOIN``-ed onto the growing
table, which is re-materialized each step. Only one 2-table join runs at a
time and the row count stays ~one-per-variant, so peak memory/disk is bounded.

The per-source dedup is what makes the observed/expected tables safe: those are
keyed per *variant* but stored per *transcript*, so a raw outer join would fan
out up to ~23x per key. Collapsing each source to one row per key first
(``bool_or`` for booleans, ``max`` for numerics -- any aggregate is safe since
there are no within-key conflicts) prevents the blow-up while preserving every
distinct variant.

Source formats
--------------
* Partitioned Parquet directories -> read via a ``**/*.parquet`` glob.
* TSV / bgzipped TSV -> ``read_csv`` with quoting disabled and ``all_varchar``;
  values are cast explicitly. The chromosome key is normalized (``chr``->``chrom``).

Run from the ``src/`` directory::

    python create_variant_eval_all_table.py
    python create_variant_eval_all_table.py --memory-limit 20GB
    python create_variant_eval_all_table.py --dry-run

Source selection (mutually exclusive; match by file_path, basename, or
substring)::

    # merge everything except the three obs/exp sources
    python create_variant_eval_all_table.py --exclude _obs_exp_mis
    # merge only the two named sources
    python create_variant_eval_all_table.py --include clinvar.tsv.bgz gnomad_independent_mis.tsv.bgz
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import duckdb

import variant_dtypes as vdt

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: merge_v2/src/...)
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # merge_v2/
INPUT_JSON = PROJECT_DIR / "data_config" / "evals_input_data.json"
DEFAULT_OUTPUT = PROJECT_DIR / "data" / "processed_data" / "evals" / "variant_evals_all.parquet"

# Canonical join key -> accepted source column aliases (first match wins).
KEY_SPEC: dict[str, list[str]] = {
    "chrom": ["chrom", "chr"],
    "pos": ["pos"],
    "ref": ["ref"],
    "alt": ["alt"],
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
    ``file_path`` (trailing slashes ignored). So ``asd_obs_exp_mis``,
    ``asd_obs_exp_mis.parquet``, ``_obs_exp_mis`` and the full relative path all
    select the asd obs/exp source.
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


def _alleles_clean(alleles_col: str) -> str:
    """Strip Hail's ``["A","G"]`` array text down to a bare ``A,G`` string."""
    expr = f"TRY_CAST({q(alleles_col)} AS VARCHAR)"
    for ch in ("[", "]", '"', " "):
        expr = f"replace({expr}, {sql_str(ch)}, '')"
    return expr


def derive_key_from_locus_alleles(canon: str, locus_col: str, alleles_col: str) -> str:
    """Derive a canonical key from Hail-style ``locus``/``alleles`` columns.

    ``locus`` looks like ``chr1:925942`` (split on ``:`` -> chrom, pos) and
    ``alleles`` like ``["A","G"]`` (cleaned to ``A,G``, split on ``,`` -> ref, alt).
    Used only for sources (the hail_evaluation tables) that lack a direct
    ``chrom/pos/ref/alt`` quartet. The chromosome keeps its ``chr`` prefix to
    match the sources that carry ``chrom`` directly.
    """
    locus = f"TRY_CAST({q(locus_col)} AS VARCHAR)"
    if canon == "chrom":
        return f"split_part({locus}, ':', 1)"
    if canon == "pos":
        return f"TRY_CAST(split_part({locus}, ':', 2) AS BIGINT)"
    clean = _alleles_clean(alleles_col)
    if canon == "ref":
        return f"split_part({clean}, ',', 1)"
    if canon == "alt":
        return f"split_part({clean}, ',', 2)"
    raise ValueError(f"unknown canonical key '{canon}'")


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
    con: duckdb.DuckDBPyConnection,
    entry: dict,
    compact: bool = False,
) -> tuple[str, dict[str, str], dict[str, tuple[str, str]], str | None]:
    """Introspect one source: reader, key casts, per-field (cast, agg), SNV filter.

    Parameters
    ----------
    compact : bool
        When true, keys are emitted in the compact numeric profile (``chrom``
        UTINYINT, ``pos`` UINTEGER; ``ref``/``alt`` the single ASCII byte
        ``UTINYINT``, with non-SNV rows dropped at build time -- see
        ``dedup_source_sql``); ``is_pos`` labels stay ``BOOLEAN`` and numeric
        eval fields stay single-precision ``FLOAT``. When false (default), the
        historical types are kept (VARCHAR/BIGINT keys) and all rows retained.

    Returns
    -------
    reader : str
        DuckDB table-function expression reading the file.
    key_casts : dict[canonical_key -> sql_expr]
        Cast yielding each canonical key in the active type profile.
    field_specs : dict[output_name -> (cast_expr, agg_func)]
        Per eval field: the cast expression and the dedup aggregate to use.
    snv_pred : str | None
        A predicate over the **raw VARCHAR** alleles keeping only SNVs (compact
        mode), or ``None`` -- applied before ``ref``/``alt`` are encoded.
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
    types = source_column_types(con, reader)
    available = set(types)

    # Prefer a direct chrom/pos/ref/alt quartet; otherwise fall back to deriving
    # the key from Hail-style locus/alleles columns (the hail_evaluation tables).
    has_locus_alleles = "locus" in available and "alleles" in available
    key_casts: dict[str, str] = {}
    raw_keys: dict[str, str] = {}
    for canon, aliases in KEY_SPEC.items():
        src_col = next((a for a in aliases if a in available), None)
        if src_col is not None:
            cast_type = "BIGINT" if canon == "pos" else "VARCHAR"
            raw = f"TRY_CAST({q(src_col)} AS {cast_type})"
        elif has_locus_alleles:
            raw = derive_key_from_locus_alleles(canon, "locus", "alleles")
        else:
            sys.exit(
                f"ERROR: {file_path}: no column for join key '{canon}' "
                f"(looked for {aliases}, or a locus+alleles pair). "
                f"Available: {sorted(available)}"
            )
        raw_keys[canon] = raw
        key_casts[canon] = vdt.key_expr(canon, raw, compact)

    snv_pred = vdt.snv_only_predicate(raw_keys["ref"], raw_keys["alt"], compact)

    field_specs: dict[str, tuple[str, str]] = {}
    for mapping in entry.get("eval_fields", []):
        for src_col, out_col in mapping.items():
            if src_col not in available:
                sys.exit(
                    f"ERROR: {file_path}: eval column '{src_col}' not found. "
                    f"Available: {sorted(available)}"
                )
            if _is_bool_field(out_col, types[src_col]):
                field_specs[out_col] = vdt.bool_field_spec(_bool_cast(src_col), compact)
            else:
                field_specs[out_col] = (f"TRY_CAST({q(src_col)} AS FLOAT)", "max")

    return reader, key_casts, field_specs, snv_pred


def dedup_source_sql(
    reader: str, key_casts: dict[str, str], field_specs: dict[str, tuple[str, str]],
    snv_pred: str | None = None,
) -> str:
    """Per-source subquery: keys + this source's fields, deduped to one row/key.

    Casts keys/fields, drops rows with an incomplete key, then collapses to one
    row per ``(chrom, pos, ref, alt)`` using each field's aggregate. Deduping
    each source *before* joining keeps the row count ~one-per-variant and stops
    a duplicated key (e.g. per-transcript observed/expected rows) from
    multiplying rows across the join chain.

    ``snv_pred`` (compact mode) drops non-SNV rows on the **raw** alleles in the
    inner scan -- before ``ref``/``alt`` are encoded to a single ASCII byte.
    """
    cast_cols = [f"{key_casts[k]} AS {q(k)}" for k in JOIN_KEYS]
    cast_cols += [f"{cast} AS {q(name)}" for name, (cast, _agg) in field_specs.items()]
    inner = "SELECT " + ", ".join(cast_cols) + f" FROM {reader}"
    if snv_pred:
        inner += f" WHERE {snv_pred}"

    keys_sql = ", ".join(q(k) for k in JOIN_KEYS)
    agg_sql = ", ".join(f"{agg}({q(name)}) AS {q(name)}" for name, (_c, agg) in field_specs.items())
    select_list = f"{keys_sql}, {agg_sql}" if agg_sql else keys_sql
    nonnull = " AND ".join(f"{q(k)} IS NOT NULL" for k in JOIN_KEYS)
    return f"SELECT {select_list} FROM (\n  {inner}\n) WHERE {nonnull} GROUP BY {keys_sql}"


def build_statements(
    sources: list[tuple[str, dict[str, str], dict[str, tuple[str, str]], str | None]],
    final_table: str,
) -> tuple[list[str], list[str]]:
    """Build the sequential, materialized full-outer-join statements.

    Materializes the wide table one source at a time::

        CREATE TABLE stage0 AS <deduped source 0>;
        CREATE TABLE stage1 AS
            SELECT * FROM stage0 FULL JOIN <deduped source 1> USING (keys);
        DROP TABLE stage0;
        ... (repeat) ...
        ALTER TABLE stageN RENAME TO <final_table>;

    Only one 2-table join runs at a time, so peak memory/disk is bounded.
    """
    field_names: list[str] = []
    for _reader, _keys, field_specs, _snv in sources:
        field_names.extend(field_specs.keys())

    using = "(" + ", ".join(q(k) for k in JOIN_KEYS) + ")"
    statements: list[str] = []

    reader0, keys0, fields0, snv0 = sources[0]
    statements.append(
        f"CREATE TABLE stage0 AS\n{dedup_source_sql(reader0, keys0, fields0, snv0)}"
    )

    prev = "stage0"
    for i in range(1, len(sources)):
        reader, keys, fields, snv = sources[i]
        cur = f"stage{i}"
        statements.append(
            f"CREATE TABLE {cur} AS\nSELECT * FROM {prev}\n"
            f"FULL JOIN (\n{dedup_source_sql(reader, keys, fields, snv)}\n) USING {using}"
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
        "--compact-dtypes", action=argparse.BooleanOptionalAction, default=False,
        help=vdt.FLAG_HELP,
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

    final_table = "variant_evals_all"
    copy_sql = (
        f"COPY (SELECT * FROM {q(final_table)}) TO {sql_str(str(args.output))} "
        f"(FORMAT PARQUET, COMPRESSION {sql_str(args.compression)}, "
        f"ROW_GROUP_SIZE {int(args.row_group_size)})"
    )

    con = duckdb.connect(str(db_path))
    try:
        configure(con, args.memory_limit, args.threads, temp_dir)
        sources = [
            analyze_source(con, entry, args.compact_dtypes)
            for entry in entries
        ]
        statements, field_names = build_statements(sources, final_table)

        profile = "compact numeric" if args.compact_dtypes else "default"
        print(f"Merging {len(entries)} source(s) -> {len(field_names)} eval field(s):")
        for _reader, _keys, field_specs, _snv_pred in sources:
            for name in field_specs:
                kind = "label" if name.lower().startswith("is_pos") else "float"
                print(f"  - {name} ({kind})")
        print(f"Dtype profile: {profile}")
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
