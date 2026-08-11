#!/usr/bin/env python3
"""Download source data from GCP into the local data tree.

This script reads the GCS locations declared in
``data_config/input_data_locations.json`` and downloads them into
``data/raw_data/scores`` (for ``score_tables``),
``data/raw_data/evals`` (for ``eval_tables``),
``data/raw_data/linker`` (for the standalone ``linker_table``),
``data/raw_data/filters`` (for ``filter_tables``), or
``data/raw_data/mechanisms`` (for ``mechanism_tables``).

The ``filter_tables`` group is special: its two subcategories land in
*different* destinations -- ``variant_level`` -> ``data/raw_data/filters/variant``
and ``gene_level`` -> ``data/raw_data/filters/ensg`` -- so the variant-level and
gene/``ensg``-level filters stay cleanly separated on disk.

The ``mechanism_tables`` group is opt-in: unlike the other groups it is
**not** included in the "no flags = download everything" default. It only
downloads when ``--mechanism-tables all`` is passed explicitly. This mirrors
the downstream pipeline, where the mechanism-append step is also opt-in.

Paths are resolved relative to this file, so the script can be invoked
from anywhere, but it is intended to be run from the parallel ``src``
directory::

    python download_source_data.py                       # download everything (except mechanisms)
    python download_source_data.py --score-tables variant_level
    python download_source_data.py --eval-tables variant_level,gene_level
    python download_source_data.py --score-tables all --eval-tables none
    python download_source_data.py --linker-table all          # linker only
    python download_source_data.py --filter-tables all         # filters only
    python download_source_data.py --filter-tables variant_level
    python download_source_data.py --mechanism-tables all      # mechanisms only (opt-in)

Selection semantics:
    * With no selection flags, everything except ``mechanism_tables`` is
      downloaded (scores, evals, filters, linker). Mechanisms are opt-in.
    * If any flag is supplied, only the explicitly requested categories
      are downloaded; the unspecified groups default to ``none``. Use
      ``all`` / ``none`` (or a comma-separated list of subcategory keys such
      as ``variant_level,gene_level``) to control each group precisely.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths (resolved relative to this file: duckdb_v2/src/download_source_data.py)
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent  # duckdb_v2/
DATA_DIR = PROJECT_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw_data"
LOCATIONS_JSON = PROJECT_DIR / "data_config" / "input_data_locations.json"

# Top-level group in the JSON -> local destination directory. Every
# subcategory of these groups lands in the same group directory.
GROUP_DESTINATIONS = {
    "score_tables": RAW_DATA_DIR / "scores",
    "eval_tables": RAW_DATA_DIR / "evals",
}

# The filter group is handled specially: each subcategory has its *own*
# destination subfolder (variant-level vs gene/ensg-level), under raw_data
# like the other downloaded source groups.
FILTER_KEY = "filter_tables"
FILTER_DEST_BASE = RAW_DATA_DIR / "filters"
FILTER_SUBCAT_DESTINATIONS = {
    "variant_level": FILTER_DEST_BASE / "variant",
    "gene_level": FILTER_DEST_BASE / "ensg",
}

# Standalone single-URI entry in the JSON -> local destination directory.
LINKER_KEY = "linker_table"
LINKER_DEST = RAW_DATA_DIR / "linker"

# Standalone flat-list entry in the JSON -> local destination directory.
# Mechanisms are opt-in: they are NOT downloaded by the default "no flags"
# path; the caller must pass ``--mechanism-tables all`` explicitly.
MECHANISM_KEY = "mechanism_tables"
MECHANISM_DEST = RAW_DATA_DIR / "mechanisms"

# Sentinel selection keywords.
SELECT_ALL = "all"
SELECT_NONE = "none"


def load_locations() -> dict:
    """Load and minimally validate the input locations JSON."""
    if not LOCATIONS_JSON.exists():
        sys.exit(f"ERROR: locations file not found: {LOCATIONS_JSON}")
    with LOCATIONS_JSON.open() as fh:
        data = json.load(fh)
    for group in GROUP_DESTINATIONS:
        if group not in data:
            sys.exit(f"ERROR: expected top-level key '{group}' in {LOCATIONS_JSON.name}")
    return data


def parse_selection(raw: str | None, available: list[str]) -> list[str]:
    """Resolve a ``--score-tables``/``--eval-tables`` argument to subcategory keys.

    ``raw`` is the comma-separated CLI value (or ``None`` if the flag was not
    passed). ``available`` is the list of subcategory keys present in the JSON
    for this group. Returns the (validated) list of keys to download.
    """
    if raw is None or raw.strip().lower() == SELECT_ALL:
        return list(available)
    if raw.strip().lower() == SELECT_NONE:
        return []

    requested = [item.strip() for item in raw.split(",") if item.strip()]
    unknown = [item for item in requested if item not in available]
    if unknown:
        sys.exit(
            f"ERROR: unknown subcategory(ies) {unknown}. "
            f"Available options: {available + [SELECT_ALL, SELECT_NONE]}"
        )
    return requested


def collect_sources(
    data: dict,
    score_selection: list[str],
    eval_selection: list[str],
    filter_selection: list[str],
    include_linker: bool,
    include_mechanisms: bool,
) -> list[tuple[str, Path]]:
    """Build a flat list of (gcs_uri, destination_dir) pairs to download."""
    selections = {
        "score_tables": score_selection,
        "eval_tables": eval_selection,
    }
    jobs: list[tuple[str, Path]] = []
    for group, subcats in selections.items():
        dest_dir = GROUP_DESTINATIONS[group]
        for subcat in subcats:
            for uri in data[group].get(subcat, []):
                jobs.append((uri, dest_dir))
    # Filters route each subcategory to its own destination subfolder.
    for subcat in filter_selection:
        dest_dir = FILTER_SUBCAT_DESTINATIONS[subcat]
        for uri in data.get(FILTER_KEY, {}).get(subcat, []):
            jobs.append((uri, dest_dir))
    if include_linker and data.get(LINKER_KEY):
        jobs.append((data[LINKER_KEY], LINKER_DEST))
    if include_mechanisms:
        for uri in data.get(MECHANISM_KEY, []):
            jobs.append((uri, MECHANISM_DEST))
    return jobs


def download(uri: str, dest_dir: Path, dry_run: bool) -> bool:
    """Download a single GCS object/prefix into ``dest_dir``.

    Uses ``gcloud storage cp -r`` so that both single files (e.g. ``*.tsv``,
    ``*.tsv.bgz``) and partitioned directory prefixes (e.g. ``*.parquet/``)
    are handled uniformly. gcloud creates the appropriately named entry inside
    ``dest_dir``. Returns True on success.
    """
    # Normalize trailing slash so directory prefixes land in a named subfolder.
    normalized = uri.rstrip("/")
    cmd = ["gcloud", "storage", "cp", "-r", normalized, f"{dest_dir}/"]
    print(f"  -> {normalized}\n     into {dest_dir}/")
    if dry_run:
        print(f"     [dry-run] {' '.join(cmd)}")
        return True
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"     FAILED (exit {result.returncode}): {uri}", file=sys.stderr)
        return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--score-tables",
        "--score_tables",
        dest="score_tables",
        metavar="CSV",
        default=None,
        help="Comma-separated score subcategories to download "
        "(e.g. 'variant_level,gene_level'), or 'all' / 'none'.",
    )
    parser.add_argument(
        "--eval-tables",
        "--eval_tables",
        dest="eval_tables",
        metavar="CSV",
        default=None,
        help="Comma-separated eval subcategories to download "
        "(e.g. 'variant_level,gene_level'), or 'all' / 'none'.",
    )
    parser.add_argument(
        "--filter-tables",
        "--filter_tables",
        dest="filter_tables",
        metavar="CSV",
        default=None,
        help="Comma-separated filter subcategories to download "
        "(e.g. 'variant_level,gene_level'), or 'all' / 'none'. "
        "variant_level -> data/raw_data/filters/variant, "
        "gene_level -> data/raw_data/filters/ensg.",
    )
    parser.add_argument(
        "--linker-table",
        "--linker_table",
        dest="linker_table",
        metavar="all|none",
        default=None,
        help="Whether to download the linker table: 'all' or 'none'.",
    )
    parser.add_argument(
        "--mechanism-tables",
        "--mechanism_tables",
        dest="mechanism_tables",
        metavar="all|none",
        default=None,
        help="Whether to download the mechanism tables: 'all' or 'none'. "
        "Opt-in: NOT included in the default 'no flags = download everything' "
        "path; must be requested explicitly. Downloads to "
        "data/raw_data/mechanisms.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the gcloud commands without downloading anything.",
    )
    args = parser.parse_args()

    if shutil.which("gcloud") is None and not args.dry_run:
        sys.exit(
            "ERROR: 'gcloud' not found on PATH. Install the Google Cloud SDK "
            "and run 'gcloud auth login' before downloading."
        )

    data = load_locations()

    score_available = list(data["score_tables"].keys())
    eval_available = list(data["eval_tables"].keys())
    # Only filter subcategories we know how to route are downloadable here.
    filter_available = [
        k for k in data.get(FILTER_KEY, {}) if k in FILTER_SUBCAT_DESTINATIONS
    ]

    # Default (no flags) downloads scores/evals/filters/linker. Mechanisms are
    # opt-in and must be requested explicitly via --mechanism-tables (see
    # mechanism_raw below). --mechanism-tables IS included in any_flag so
    # requesting only mechanisms suppresses the default-all of the other
    # groups, matching the behavior of the other selection flags.
    any_flag = (
        args.score_tables is not None
        or args.eval_tables is not None
        or args.filter_tables is not None
        or args.linker_table is not None
        or args.mechanism_tables is not None
    )
    score_raw = args.score_tables if args.score_tables is not None else (SELECT_NONE if any_flag else SELECT_ALL)
    eval_raw = args.eval_tables if args.eval_tables is not None else (SELECT_NONE if any_flag else SELECT_ALL)
    filter_raw = args.filter_tables if args.filter_tables is not None else (SELECT_NONE if any_flag else SELECT_ALL)
    linker_raw = args.linker_table if args.linker_table is not None else (SELECT_NONE if any_flag else SELECT_ALL)
    # Mechanisms are opt-in regardless of any_flag: they only download when
    # explicitly requested.
    mechanism_raw = args.mechanism_tables if args.mechanism_tables is not None else SELECT_NONE

    score_selection = parse_selection(score_raw, score_available)
    eval_selection = parse_selection(eval_raw, eval_available)
    filter_selection = parse_selection(filter_raw, filter_available)

    linker_choice = linker_raw.strip().lower()
    if linker_choice not in (SELECT_ALL, SELECT_NONE):
        sys.exit(
            f"ERROR: invalid --linker-table value '{linker_raw}'. "
            f"Available options: {[SELECT_ALL, SELECT_NONE]}"
        )
    include_linker = linker_choice == SELECT_ALL

    mechanism_choice = mechanism_raw.strip().lower()
    if mechanism_choice not in (SELECT_ALL, SELECT_NONE):
        sys.exit(
            f"ERROR: invalid --mechanism-tables value '{mechanism_raw}'. "
            f"Available options: {[SELECT_ALL, SELECT_NONE]}"
        )
    include_mechanisms = mechanism_choice == SELECT_ALL

    # Ensure destination directories exist.
    for dest_dir in GROUP_DESTINATIONS.values():
        dest_dir.mkdir(parents=True, exist_ok=True)
    for subcat in filter_selection:
        FILTER_SUBCAT_DESTINATIONS[subcat].mkdir(parents=True, exist_ok=True)
    if include_linker:
        LINKER_DEST.mkdir(parents=True, exist_ok=True)
    if include_mechanisms:
        MECHANISM_DEST.mkdir(parents=True, exist_ok=True)

    jobs = collect_sources(
        data,
        score_selection,
        eval_selection,
        filter_selection,
        include_linker,
        include_mechanisms,
    )
    if not jobs:
        print("Nothing selected to download. Use --score-tables / --eval-tables / "
              "--filter-tables / --linker-table / --mechanism-tables with "
              "'all', 'none', or specific subcategory keys.")
        return 0

    print(f"Selected {len(jobs)} item(s) to download:")
    print(f"  score_tables:     {score_selection or 'none'}")
    print(f"  eval_tables:      {eval_selection or 'none'}")
    print(f"  filter_tables:    {filter_selection or 'none'}")
    print(f"  linker_table:     {'all' if include_linker else 'none'}")
    print(f"  mechanism_tables: {'all' if include_mechanisms else 'none'}\n")

    failures: list[str] = []
    for idx, (uri, dest_dir) in enumerate(jobs, start=1):
        print(f"[{idx}/{len(jobs)}]")
        if not download(uri, dest_dir, args.dry_run):
            failures.append(uri)

    print("\nDone.")
    if failures:
        print(f"{len(failures)} download(s) failed:", file=sys.stderr)
        for uri in failures:
            print(f"  {uri}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
