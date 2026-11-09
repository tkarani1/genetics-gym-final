#!/usr/bin/env python3
"""Append per-chromosome mechanism annotations to a gene-aggregated analysis table.

Reads a gene-aggregated parquet keyed on ``(chrom, pos, ref, alt, ensg)`` and
LEFT-JOINs the per-chromosome mechanism shards under
``../data/raw_data/mechanisms/`` (``linker_mech_chrom_{1..22,X}.tsv.gz``) onto
it, appending the mech table's non-key columns to the output. Every non-key
mech column is emitted with a ``mech_`` prefix so it cannot collide with
existing filter / eval columns (e.g. ``aromatic`` -> ``mech_aromatic``,
``foldx_ddg`` -> ``mech_foldx_ddg``).

Iteration is **per chromosome**: one chromosome slice + one mech shard at a
time, each written to a per-chrom Parquet checkpoint under
``.ckpt_<output stem>/``. After all shards are processed, the checkpoints are
concatenated into the final output. This keeps peak memory bounded to a single
chromosome's slice + its mech shard, and lets a crash lose at most one
chromosome's worth of work. On restart, chromosomes whose checkpoint already
exists are skipped.

Rows in the input whose ``chrom`` is not one of the 23 mech shards (e.g.
``chrY``, ``chrM``) are preserved with NULL mech columns via a final
**residual step** that writes ``_residual.parquet`` alongside the per-chrom
checkpoints; the consolidation includes it after ``chrX.parquet``. This
keeps the LEFT JOIN contract (no input row is dropped) end-to-end even
though the per-chromosome loop can only visit chromosomes that have a mech
shard.

All 5 join key columns must be present in the input; a missing key column is
a hard error.

The mech schema is introspected **once** from the first shard (``chr1``) with
a full-file type inference sample, then that ``{name: sql_type}`` mapping is
pinned on every subsequent ``read_csv`` call. This guarantees every per-chrom
checkpoint has the same output columns / types, so the final consolidation
concatenation is schema-stable.

Deduplication
-------------
The mech shards are *not* 1-to-1 on ``(chrom, pos, ref, alt, ensg)``: each
5-key can appear twice (typically once with an all-NULL payload and once with
the real values -- a residue of the source table having an additional
distinguishing column, e.g. transcript/isoform, that was dropped when the
mech linker was materialised). To keep the LEFT JOIN 1-to-1 (no fan-out) the
mech shard is deduped in-place inside the join subquery: ``GROUP BY`` the
5-key and aggregate every non-key column with ``max()`` (NULL-skipping;
picks the filled value when the group is one-NULL + one-filled, which is the
observed pattern for every duplicate group). For BOOLEAN columns ``max`` acts
as ``bool_or`` (``FALSE < TRUE``); for numerics it picks the non-null value.

Run from ``merge_v2/src/`` (requires the input analysis table and the mech
shards to be present locally first)::

    python download_source_data.py --mechanism-tables all   # one-time
    python create_gene_aggregated_mechanism_appended_table.py
    python create_gene_aggregated_mechanism_appended_table.py --input <path>.parquet
    python create_gene_aggregated_mechanism_appended_table.py --dry-run
    python create_gene_aggregated_mechanism_appended_table.py --overwrite

Prerequisites:
    * A gene-aggregated parquet keyed on ``(chrom, pos, ref, alt, ensg)`` --
      typically the output of ``create_gene_aggregated_filtered_table.py``.
    * The mechanism shards downloaded to ``../data/raw_data/mechanisms/``.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import duckdb

# Reuse the small SQL helpers (identifier/string quoting, parquet reader,
# DuckDB configuration, COPY builder, CSV null sentinels, key constants) from
# the variant-level filter script so this module speaks the same dialect.
import create_variant_scores_filtered_tables as vf  # noqa: E402
import variant_dtypes as vdt  # noqa: E402

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: merge_v2/src/...).
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # merge_v2/
DATA_DIR = PROJECT_DIR / "data"
GENE_AGG_DIR = (
    DATA_DIR / "processed_data" / "full_analysis_tables" / "gene_aggregated"
)
MECHANISMS_DIR = DATA_DIR / "raw_data" / "mechanisms"

DEFAULT_INPUT = (
    GENE_AGG_DIR / "variant_scores_all_outer_ensg_stats_eval_filtered.parquet"
)

# Mechanism file naming and mech column output convention.
MECH_FILE_PREFIX = "linker_mech_chrom_"
MECH_FILE_SUFFIX = ".tsv.gz"
MECH_COL_PREFIX = "mech_"

# Suffix added to the input stem to build the default output path.
MECH_OUTPUT_SUFFIX = "_mech"

# The 5-column join key: variant key + ensg.
JOIN_KEYS = vf.VARIANT_KEYS + [vf.ENSG_COL]

# Canonical chromosome order for both the loop and the consolidated output.
CHROMOSOME_ORDER = [f"chr{n}" for n in range(1, 23)] + ["chrX"]

# Filenames like ``linker_mech_chrom_22.tsv.gz`` / ``linker_mech_chrom_X.tsv.gz``.
_MECH_NAME_RE = re.compile(
    rf"^{re.escape(MECH_FILE_PREFIX)}([0-9]+|X){re.escape(MECH_FILE_SUFFIX)}$"
)

PLAN_FILE = "plan.json"
CKPT_PREFIX = ".ckpt_"

# Checkpoint filename for input rows whose ``chrom`` is not covered by any
# mech shard (e.g. chrY, chrM). Kept alongside the per-chromosome checkpoints
# so the final consolidation sees a schema-identical piece for every input
# row and the LEFT JOIN contract (no input row is dropped) holds end-to-end.
RESIDUAL_CHECKPOINT = "_residual.parquet"


# ---------------------------------------------------------------------------
# Discovery of the mech shards on disk.
# ---------------------------------------------------------------------------
def discover_mech_files(directory: Path) -> list[tuple[str, Path]]:
    """Return ``[(chr_label, path), ...]`` for the mech shards, in chrom order.

    Files that do not match the ``linker_mech_chrom_{N}.tsv.gz`` pattern are
    ignored (extra unrelated files in the directory don't break the loop).
    """
    found: dict[str, Path] = {}
    for p in sorted(directory.glob(f"{MECH_FILE_PREFIX}*{MECH_FILE_SUFFIX}")):
        m = _MECH_NAME_RE.match(p.name)
        if m is None:
            continue
        found["chr" + m.group(1)] = p
    return [(c, found[c]) for c in CHROMOSOME_ORDER if c in found]


# ---------------------------------------------------------------------------
# Mech-shard readers. The schema is inferred once from a reference shard
# (full-file sample) and pinned on every subsequent read via ``columns={...}``
# so all 23 checkpoints share identical column types.
# ---------------------------------------------------------------------------
def _csv_common_kwargs() -> str:
    """The reader kwargs shared by every mech ``read_csv`` call."""
    nullstr = "[" + ", ".join(vf.sql_str(v) for v in vf.CSV_NULL_VALUES) + "]"
    return (
        "delim='\\t', header=true, quote='', "
        f"nullstr={nullstr}, "
        "compression='gzip'"
    )


def _mech_reader_auto(path: Path) -> str:
    """Basic ``read_csv`` with full-file type inference; used only at planning."""
    return (
        "read_csv("
        f"{vf.sql_str(str(path))}, "
        f"{_csv_common_kwargs()}, "
        "sample_size=-1)"
    )


def _mech_reader_pinned(path: Path, schema: list[tuple[str, str]]) -> str:
    """``read_csv`` with an explicit ``columns={...}`` map pinning every dtype."""
    columns_map = (
        "{"
        + ", ".join(
            f"{vf.sql_str(name)}: {vf.sql_str(sql_type)}"
            for name, sql_type in schema
        )
        + "}"
    )
    return (
        "read_csv("
        f"{vf.sql_str(str(path))}, "
        f"{_csv_common_kwargs()}, "
        f"columns={columns_map})"
    )


def introspect_mech_schema(
    con: duckdb.DuckDBPyConnection, reference_path: Path
) -> list[tuple[str, str]]:
    """``[(col_name, sql_type), ...]`` for the reference shard, in file order.

    The 5 key columns (chrom, pos, ref, alt, ensg) come first and are followed
    by the mech's non-key annotation columns.
    """
    rows = con.execute(
        f"DESCRIBE SELECT * FROM {_mech_reader_auto(reference_path)}"
    ).fetchall()
    return [(r[0], r[1]) for r in rows]


def mech_value_cols(schema: list[tuple[str, str]]) -> list[str]:
    """Non-key mech columns in file order (the appended payload)."""
    keys = set(JOIN_KEYS)
    return [name for name, _ in schema if name not in keys]


def prefixed_mech_cols(schema: list[tuple[str, str]]) -> list[str]:
    """The mech value columns after applying the ``mech_`` output prefix."""
    return [MECH_COL_PREFIX + c for c in mech_value_cols(schema)]


# ---------------------------------------------------------------------------
# Per-chromosome SQL builder.
# ---------------------------------------------------------------------------
def _mech_dedup_sql(mech_reader: str, schema: list[tuple[str, str]], compact: bool) -> str:
    """Dedup a mech shard on the 5-key, aggregating non-key columns with max().

    ``max`` skips NULLs, so on the observed "one-NULL + one-filled" duplicate
    groups it deterministically picks the filled value. For BOOLEAN columns it
    behaves as ``bool_or`` (``FALSE < TRUE``); for DOUBLE/BIGINT it picks the
    non-null value. Output columns are the 5 keys followed by the non-key mech
    columns, each already renamed with the ``mech_`` prefix.

    In the compact profile the SNV drop runs on the **raw** VARCHAR ref/alt
    inside this subquery, before the outer key encoding turns them into a
    single UTINYINT byte.
    """
    key_sel = ", ".join(vf.q(k) for k in JOIN_KEYS)
    aggs = [
        f"max({vf.q(name)}) AS {vf.q(MECH_COL_PREFIX + name)}"
        for name in mech_value_cols(schema)
    ]
    snv_pred = vdt.snv_only_predicate(vf.q("ref"), vf.q("alt"), compact)
    where = f" WHERE {snv_pred}" if snv_pred else ""
    return (
        f"SELECT {key_sel}"
        + (f", {', '.join(aggs)}" if aggs else "")
        + f" FROM {mech_reader}{where} "
        f"GROUP BY {key_sel}"
    )


def _mech_encoded_sql(
    mech_reader: str, schema: list[tuple[str, str]], compact: bool
) -> str:
    """The mech-shard subquery ready to be joined: deduped, then key-encoded.

    In the default profile the key encoding is a no-op (raw VARCHAR/BIGINT
    types already match a default-profile input). In the compact profile the
    keys are re-encoded to the input's UTINYINT/UINTEGER shape.
    """
    dedup = _mech_dedup_sql(mech_reader, schema, compact)
    prefixed_cols = [MECH_COL_PREFIX + n for n in mech_value_cols(schema)]
    key_projection = [
        f"{vdt.key_expr(k, vf.q(k), compact)} AS {vf.q(k)}"
        for k in vf.VARIANT_KEYS
    ] + [f"{vf.q(vf.ENSG_COL)} AS {vf.q(vf.ENSG_COL)}"]
    value_passthrough = [vf.q(c) for c in prefixed_cols]
    return (
        f"SELECT {', '.join(key_projection)}"
        + (f", {', '.join(value_passthrough)}" if value_passthrough else "")
        + f" FROM (\n{dedup}\n)"
    )


def _chrom_literal(chrom_label: str, compact: bool) -> str:
    """SQL literal matching the input's ``chrom`` column for the given label.

    Default profile: ``'chr22'`` (VARCHAR). Compact profile: the UTINYINT
    encoding used by ``variant_dtypes`` (``chr1``..``chr22`` -> 1..22,
    ``chrX`` -> 23).
    """
    if not compact:
        return vf.sql_str(chrom_label)
    if chrom_label == "chrX":
        return str(vdt.CHR_X)
    return chrom_label.removeprefix("chr")


def _keys_using() -> str:
    return "(" + ", ".join(vf.q(k) for k in JOIN_KEYS) + ")"


def chromosome_join_sql(
    input_reader: str,
    mech_reader: str,
    schema: list[tuple[str, str]],
    chrom_label: str,
    compact: bool,
) -> str:
    """The full SELECT that produces one chromosome's slice-and-LEFT-JOIN.

    Shape::

        SELECT * FROM (SELECT * FROM <input> WHERE chrom = <lit>)
        LEFT JOIN (
            SELECT <encoded keys>, mech_<col1>, mech_<col2>, ...
            FROM (SELECT chrom,pos,ref,alt,ensg, max(col1) AS mech_col1, ...
                  FROM <mech> [WHERE <snv>] GROUP BY 5-key)
        ) USING (chrom, pos, ref, alt, ensg)

    The inner GROUP BY dedupes the mech shard so the LEFT JOIN cannot fan out;
    the outer SELECT re-encodes the mech's raw keys to match the input's
    profile (a no-op in default mode, UTINYINT/UINTEGER in compact mode). The
    USING clause emits each join key exactly once, so output columns are the
    input's columns followed by the ``mech_``-prefixed non-key mech columns
    in mech-shard file order.
    """
    chrom_lit = _chrom_literal(chrom_label, compact)
    input_slice = (
        f"SELECT * FROM {input_reader} "
        f"WHERE {vf.q('chrom')} = {chrom_lit}"
    )
    mech_subq = _mech_encoded_sql(mech_reader, schema, compact)
    return (
        f"SELECT * FROM (\n{input_slice}\n) "
        f"LEFT JOIN (\n{mech_subq}\n) USING {_keys_using()}"
    )


# ---------------------------------------------------------------------------
# Checkpoint bookkeeping.
# ---------------------------------------------------------------------------
def ckpt_dir_for(output: Path) -> Path:
    return output.parent / f"{CKPT_PREFIX}{output.stem}"


def ckpt_path_for(ckpt_dir: Path, chrom_label: str) -> Path:
    return ckpt_dir / f"{chrom_label}.parquet"


def build_plan_dict(
    args: argparse.Namespace,
    mech_files: list[tuple[str, Path]],
    schema: list[tuple[str, str]],
) -> dict:
    """Serialisable resume-validation plan; the checkpoint dir is rejected if
    this dictionary no longer matches on a subsequent run."""
    return {
        "input": str(args.input),
        "output": str(args.output),
        "compact_dtypes": bool(args.compact_dtypes),
        "chromosomes": [c for c, _ in mech_files],
        "mech_files": {c: p.name for c, p in mech_files},
        "mech_schema": [(name, sql_type) for name, sql_type in schema],
        "mech_col_prefix": MECH_COL_PREFIX,
    }


def write_plan(ckpt_dir: Path, plan: dict) -> None:
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    (ckpt_dir / PLAN_FILE).write_text(json.dumps(plan, indent=2))


def load_plan(ckpt_dir: Path) -> dict | None:
    p = ckpt_dir / PLAN_FILE
    if not p.exists():
        return None
    return json.loads(p.read_text())


def validate_resume(ckpt_dir: Path, plan: dict) -> None:
    """Refuse to resume against a plan that no longer matches the inputs."""
    existing = load_plan(ckpt_dir)
    if existing is None:
        return
    # JSON round-trips tuples to lists; normalise before comparing.
    expected = json.loads(json.dumps(plan))
    if existing != expected:
        sys.exit(
            "ERROR: existing checkpoints in\n"
            f"  {ckpt_dir}\n"
            "were built with a different plan (input path, output path, mech "
            "schema, chromosome list, or --compact-dtypes flag changed). "
            "Re-run with --overwrite to start clean, or restore the original "
            "inputs."
        )


# ---------------------------------------------------------------------------
# DuckDB execution helpers (single-purpose wrappers around the sibling's
# ``configure`` / ``copy_sql``).
# ---------------------------------------------------------------------------
def _residual_select_sql(
    input_reader: str,
    schema: list[tuple[str, str]],
    mech_chromosomes: list[str],
    compact: bool,
) -> str:
    """SELECT for the residual: every input row whose chrom is not covered by
    any mech shard, projected with NULL mech columns cast to their pinned
    types (so the schema matches the per-chrom checkpoints exactly)."""
    lits = ", ".join(_chrom_literal(c, compact) for c in mech_chromosomes)
    null_cols = [
        f"CAST(NULL AS {sql_type}) AS {vf.q(MECH_COL_PREFIX + name)}"
        for name, sql_type in schema
        if name not in set(JOIN_KEYS)
    ]
    return (
        "SELECT *"
        + (f", {', '.join(null_cols)}" if null_cols else "")
        + f" FROM {input_reader}"
        + f" WHERE {vf.q('chrom')} NOT IN ({lits})"
    )


def run_residual_batch(
    input_path: Path,
    schema: list[tuple[str, str]],
    mech_chromosomes: list[str],
    out_path: Path,
    temp_dir: Path,
    args: argparse.Namespace,
) -> int:
    """Write the residual checkpoint (input rows with chrom outside the mech
    shard set). May legitimately be zero rows if the input has no such rows;
    we still emit an empty Parquet so the consolidation's UNION ALL sees a
    schema-consistent piece for every position in the parts list."""
    shutil.rmtree(temp_dir, ignore_errors=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    input_reader = vf.parquet_reader(input_path)
    select_sql = _residual_select_sql(
        input_reader, schema, mech_chromosomes, args.compact_dtypes
    )
    con = duckdb.connect()
    try:
        vf.configure(con, args.memory_limit, args.threads, temp_dir)
        print(
            f"    streaming residual (chrom NOT IN mech shards) -> "
            f"{out_path.name} ...",
            flush=True,
        )
        con.execute(
            vf.copy_sql(select_sql, out_path, args.compression, args.row_group_size)
        )
        n_rows = con.execute(
            f"SELECT count(*) FROM {vf.parquet_reader(out_path)}"
        ).fetchone()[0]
    finally:
        con.close()
    shutil.rmtree(temp_dir, ignore_errors=True)
    return n_rows


def run_chromosome_batch(
    input_path: Path,
    mech_path: Path,
    schema: list[tuple[str, str]],
    chrom_label: str,
    out_path: Path,
    temp_dir: Path,
    args: argparse.Namespace,
) -> int:
    """Stream one chromosome slice + LEFT JOIN mech shard as one ``COPY`` and
    return the row count of the resulting Parquet checkpoint."""
    shutil.rmtree(temp_dir, ignore_errors=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    input_reader = vf.parquet_reader(input_path)
    mech_reader = _mech_reader_pinned(mech_path, schema)
    select_sql = chromosome_join_sql(
        input_reader, mech_reader, schema, chrom_label, args.compact_dtypes
    )

    con = duckdb.connect()
    try:
        vf.configure(con, args.memory_limit, args.threads, temp_dir)
        print(
            f"    streaming {chrom_label} -> {out_path.name} "
            f"(from {mech_path.name}) ...",
            flush=True,
        )
        con.execute(
            vf.copy_sql(select_sql, out_path, args.compression, args.row_group_size)
        )
        n_rows = con.execute(
            f"SELECT count(*) FROM {vf.parquet_reader(out_path)}"
        ).fetchone()[0]
    finally:
        con.close()
    shutil.rmtree(temp_dir, ignore_errors=True)
    return n_rows


def consolidate(
    ckpt_dir: Path,
    chrom_labels: list[str],
    out_path: Path,
    temp_dir: Path,
    args: argparse.Namespace,
) -> int:
    """Concat all per-chrom checkpoints plus the residual into the final
    Parquet.

    Reads each ``<ckpt_dir>/<chrom_label>.parquet`` in canonical chromosome
    order followed by ``_residual.parquet`` via a UNION ALL
    ``read_parquet(list)``, so DuckDB streams them sequentially rather than
    treating them as a partitioned dataset (avoids any glob-order surprises).
    """
    shutil.rmtree(temp_dir, ignore_errors=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    parts = [str(ckpt_path_for(ckpt_dir, c)) for c in chrom_labels]
    parts.append(str(ckpt_dir / RESIDUAL_CHECKPOINT))
    parts_sql = "[" + ", ".join(vf.sql_str(p) for p in parts) + "]"
    select_sql = f"SELECT * FROM read_parquet({parts_sql})"

    con = duckdb.connect()
    try:
        vf.configure(con, args.memory_limit, args.threads, temp_dir)
        print(f"  streaming consolidation -> {out_path.name} ...", flush=True)
        con.execute(
            vf.copy_sql(select_sql, out_path, args.compression, args.row_group_size)
        )
        n_rows = con.execute(
            f"SELECT count(*) FROM {vf.parquet_reader(out_path)}"
        ).fetchone()[0]
    finally:
        con.close()
    shutil.rmtree(temp_dir, ignore_errors=True)
    return n_rows


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------
def default_output(input_path: Path) -> Path:
    """``<stem>.parquet`` -> ``<stem>_mech.parquet`` beside the input."""
    return input_path.with_name(f"{input_path.stem}{MECH_OUTPUT_SUFFIX}.parquet")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT,
        help=f"Gene-aggregated analysis parquet to append mech columns onto "
        f"(default: {DEFAULT_INPUT}).",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Output parquet path (default: <input stem>_mech.parquet beside "
        "the input).",
    )
    parser.add_argument(
        "--mechanisms-dir", "--mechanisms_dir", dest="mechanisms_dir",
        type=Path, default=MECHANISMS_DIR,
        help=f"Directory of linker_mech_chrom_*.tsv.gz shards "
        f"(default: {MECHANISMS_DIR}).",
    )
    parser.add_argument(
        "--compact-dtypes", action=argparse.BooleanOptionalAction, default=False,
        help="Match the compact numeric key profile of the input: encode the "
        "mech shard's chrom/pos/ref/alt the same way (so the LEFT JOIN keys "
        "match) and drop non-SNV rows up front. Must match the flag the input "
        "table was built with (the default gene-aggregated tables use the "
        "non-compact profile).",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Discard any existing output and checkpoints and start clean.",
    )
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Keep the per-chromosome checkpoint directory after a successful "
        "run (for debugging).",
    )
    parser.add_argument(
        "--memory-limit", default=None,
        help="DuckDB memory limit, e.g. '12GB' (default: DuckDB's ~80%% of RAM).",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="DuckDB worker threads (default: DuckDB auto).",
    )
    parser.add_argument(
        "--temp-dir", type=Path, default=None,
        help="Directory for DuckDB spill files during each per-chromosome "
        "batch (default: <output dir>/.duckdb_spill).",
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
        help="Print the per-chromosome plan (and a representative SQL "
        "statement) and exit.",
    )
    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(
            f"ERROR: input not found: {args.input}\n"
            f"       Build it with create_gene_aggregated_filtered_table.py "
            f"first, or pass --input <path>."
        )
    if not args.mechanisms_dir.exists():
        sys.exit(
            f"ERROR: mechanisms directory not found: {args.mechanisms_dir}\n"
            f"       Run: python download_source_data.py --mechanism-tables all"
        )

    mech_files = discover_mech_files(args.mechanisms_dir)
    if not mech_files:
        sys.exit(
            f"ERROR: no linker_mech_chrom_*.tsv.gz files found under\n"
            f"  {args.mechanisms_dir}\n"
            f"       Run: python download_source_data.py --mechanism-tables all"
        )

    args.output = args.output or default_output(args.input)
    out_dir = args.output.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.temp_dir or (out_dir / ".duckdb_spill")
    ckpt_dir = ckpt_dir_for(args.output)

    # Plan: introspect input & mech schemas, validate collisions, build plan
    # dict. All done through a short-lived DuckDB connection.
    plan_con = duckdb.connect()
    try:
        input_cols = vf.column_names(plan_con, vf.parquet_reader(args.input))
        reference_path = mech_files[0][1]  # chr1 = largest & most-populated shard
        mech_schema = introspect_mech_schema(plan_con, reference_path)
    finally:
        plan_con.close()

    missing_keys = [k for k in JOIN_KEYS if k not in input_cols]
    if missing_keys:
        sys.exit(
            f"ERROR: input is missing required join key column(s) "
            f"{missing_keys}: {args.input}\n"
            f"       Expected all of {JOIN_KEYS}."
        )

    prefixed = prefixed_mech_cols(mech_schema)
    input_col_lower = {c.lower(): c for c in input_cols}
    collisions = sorted(
        input_col_lower[c.lower()] for c in prefixed if c.lower() in input_col_lower
    )
    if collisions:
        sys.exit(
            f"ERROR: input already has column(s) that would collide with "
            f"prefixed mech columns: {collisions}. The mech append would "
            f"produce duplicate column names."
        )

    plan = build_plan_dict(args, mech_files, mech_schema)

    print(f"Input:  {args.input}")
    print(f"Output: {args.output}")
    print(f"Mechanisms dir: {args.mechanisms_dir}")
    print(f"Mech shards: {len(mech_files)} "
          f"({', '.join(c for c, _ in mech_files)})")
    print(f"Mech columns to append: {len(prefixed)} "
          f"(prefixed with '{MECH_COL_PREFIX}')")
    print(f"Join type: LEFT OUTER on {JOIN_KEYS}")
    print(f"Dtype profile: {'compact numeric' if args.compact_dtypes else 'default'}")
    print(f"Checkpoints: {ckpt_dir}")
    print(f"Spill dir: {temp_dir}\n")

    if args.dry_run:
        print("Per-chromosome plan:")
        for c, p in mech_files:
            print(f"  {c:6s} -> {p.name}  ({ckpt_path_for(ckpt_dir, c).name})")
        print(f"\nMech columns (all prefixed with '{MECH_COL_PREFIX}'):")
        for c in prefixed:
            print(f"  {c}")
        print("\nRepresentative SQL (first chromosome):")
        first_chrom, first_path = mech_files[0]
        sql = chromosome_join_sql(
            vf.parquet_reader(args.input),
            _mech_reader_pinned(first_path, mech_schema),
            mech_schema,
            first_chrom,
            args.compact_dtypes,
        )
        print(vf.copy_sql(sql, ckpt_path_for(ckpt_dir, first_chrom),
                          args.compression, args.row_group_size) + ";")
        print("\nResidual SQL (input rows with chrom outside mech shards):")
        residual_sql = _residual_select_sql(
            vf.parquet_reader(args.input),
            mech_schema,
            [c for c, _ in mech_files],
            args.compact_dtypes,
        )
        print(vf.copy_sql(residual_sql, ckpt_dir / RESIDUAL_CHECKPOINT,
                          args.compression, args.row_group_size) + ";")
        return 0

    if args.overwrite:
        shutil.rmtree(ckpt_dir, ignore_errors=True)
        args.output.unlink(missing_ok=True)

    if args.output.exists():
        print("Output already exists; nothing to do (use --overwrite to rebuild).")
        return 0

    validate_resume(ckpt_dir, plan)
    write_plan(ckpt_dir, plan)

    # Per-chromosome loop. Each iteration writes exactly one checkpoint; if a
    # checkpoint already exists we skip it (that's how we resume). The
    # per-chromosome memory footprint = one chrom slice + one mech shard,
    # nothing carries between iterations except the checkpoint on disk.
    total_rows = 0
    for i, (chrom_label, mech_path) in enumerate(mech_files, start=1):
        ckpt_path = ckpt_path_for(ckpt_dir, chrom_label)
        header = f"[{i}/{len(mech_files)}] {chrom_label}"
        if ckpt_path.exists():
            existing = duckdb.connect()
            try:
                n = existing.execute(
                    f"SELECT count(*) FROM {vf.parquet_reader(ckpt_path)}"
                ).fetchone()[0]
            finally:
                existing.close()
            total_rows += n
            print(f"  {header}: checkpoint exists ({n:,} rows) -- skipping.")
            continue

        print(f"  {header}: joining {mech_path.name} ...")
        n = run_chromosome_batch(
            args.input, mech_path, mech_schema, chrom_label,
            ckpt_path, temp_dir, args,
        )
        total_rows += n
        print(f"     checkpoint written: {n:,} rows -> {ckpt_path.name}")

    # Residual step: preserve any input rows whose chrom is outside the mech
    # shard set (e.g. chrY, chrM). Emitted as its own checkpoint with NULL
    # mech columns so the LEFT JOIN contract holds for every input row.
    chrom_labels = [c for c, _ in mech_files]
    residual_path = ckpt_dir / RESIDUAL_CHECKPOINT
    if residual_path.exists():
        existing = duckdb.connect()
        try:
            n = existing.execute(
                f"SELECT count(*) FROM {vf.parquet_reader(residual_path)}"
            ).fetchone()[0]
        finally:
            existing.close()
        total_rows += n
        print(f"  residual: checkpoint exists ({n:,} rows) -- skipping.")
    else:
        print("  residual: writing input rows with chrom outside mech shards ...")
        n = run_residual_batch(
            args.input, mech_schema, chrom_labels,
            residual_path, temp_dir, args,
        )
        total_rows += n
        print(f"     checkpoint written: {n:,} rows -> {residual_path.name}")

    # Consolidation: single COPY reading all per-chrom checkpoints plus the
    # residual in canonical chromosome order.
    final_rows = consolidate(ckpt_dir, chrom_labels, args.output, temp_dir, args)

    # Sanity: the consolidated file should have the same row count as the sum
    # of the per-chromosome checkpoints.
    if final_rows != total_rows:
        print(
            f"WARNING: consolidated row count {final_rows:,} != "
            f"sum of per-chrom checkpoints {total_rows:,}",
            file=sys.stderr,
        )

    if not args.keep_intermediates:
        shutil.rmtree(ckpt_dir, ignore_errors=True)
    shutil.rmtree(temp_dir, ignore_errors=True)

    final_con = duckdb.connect()
    try:
        cols = vf.column_names(final_con, vf.parquet_reader(args.output))
    finally:
        final_con.close()
    print(
        f"\nDone. Wrote {final_rows:,} rows x {len(cols)} columns to "
        f"{args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
