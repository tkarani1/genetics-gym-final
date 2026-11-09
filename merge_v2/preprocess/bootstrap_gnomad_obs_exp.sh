#!/usr/bin/env bash
# One-time bootstrap for the 4 new gnomad obs/exp eval sources.
#
# End-to-end: raw hail exports on GCS -> preprocessed *_obs_exp_mis.parquet on
# GCS, ready to be consumed by the standard eval merge pipeline (from that
# point on, download_source_data.py + create_variant_eval_all_table.py handle
# them like any other eval source; the URIs already exist in
# data_config/input_data_locations.json under eval_tables.variant_level, and
# the column projections already exist in data_config/evals_input_data.json).
#
# Steps performed per cohort:
#
#   [1/3] Download raw hail export
#         gs://.../parquet_files/scores/gnomad.<COHORT>.genetics_gym.per_variant.expected.parquet
#         -> <RAW_LOCAL_DIR>/gnomad.<COHORT>.genetics_gym.per_variant.expected.parquet/
#         Skipped if already present locally (idempotent re-run).
#
#   [2/3] Preprocess into the *_obs_exp_mis.parquet shape via
#         merge_v2/preprocess/preprocess_gnomad_obs_exp.py:
#         - filter to canonical/MANE Ensembl missense (one row per variant)
#         - split alleles[] into ref/alt, derive chrom/pos from locus.*
#         - rename observed/expected per cohort
#         Skipped if the output already exists locally (respects --overwrite).
#
#   [3/3] Upload the preprocessed file to the eval_obs_exp/ GCS prefix declared
#         in input_data_locations.json. Skipped if the GCS object already
#         exists (idempotent re-run).
#
# Local disk after full run: ~4.7 GiB of raw scratch under RAW_LOCAL_DIR plus
# the (smaller) preprocessed outputs under OUT_LOCAL_DIR. Set KEEP_RAW=0 to
# delete each raw directory after its cohort's preprocessed file is uploaded.
#
# Usage:
#   bash merge_v2/preprocess/bootstrap_gnomad_obs_exp.sh                  # do everything
#   DRY_RUN=1 bash merge_v2/preprocess/bootstrap_gnomad_obs_exp.sh        # print, don't run
#   KEEP_RAW=0 bash merge_v2/preprocess/bootstrap_gnomad_obs_exp.sh       # clean raw after upload
#   OVERWRITE=1 bash merge_v2/preprocess/bootstrap_gnomad_obs_exp.sh      # rebuild + reupload
#   COHORTS="gnomad_new gnomad_v2" bash ...                                # subset

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths (resolved relative to this script: merge_v2/preprocess/bootstrap_gnomad_obs_exp.sh)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_DIR="$( cd -- "$SCRIPT_DIR/.." &> /dev/null && pwd )"   # merge_v2/
RAW_LOCAL_DIR="$PROJECT_DIR/data/raw_data/scores"        # matches download_source_data.py conventions
OUT_LOCAL_DIR="$PROJECT_DIR/data/raw_data/evals"         # matches evals_input_data.json file_path
VENV_PYTHON="$PROJECT_DIR/.venv/bin/python"
PREPROC_SCRIPT="$SCRIPT_DIR/preprocess_gnomad_obs_exp.py"

# GCS prefixes. RAW_GCS_PREFIX is where the hail exports live; OUT_GCS_PREFIX
# is the eval_obs_exp/ prefix already referenced in input_data_locations.json.
RAW_GCS_PREFIX="gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/scores"
OUT_GCS_PREFIX="gs://grohlicek/genetics_gym_vsm_all_content/eval_obs_exp"

# Cohort <TAG> -> raw <BASENAME> mapping. Must match the COHORTS list baked
# into preprocess_gnomad_obs_exp.py so both agents see the same set. The
# preprocessed output is always <TAG>_obs_exp_mis.parquet, matching the eval
# JSON entries. A case statement is used (not an associative array) because
# macOS still ships bash 3.2, which pre-dates declare -A.
declare -a COHORT_TAGS=(gnomad_new gnomad_ukb gnomad_v2 gnomad_v2_ukb)

raw_basename_for() {
  case "$1" in
    gnomad_new)    echo "gnomad.new.genetics_gym.per_variant.expected.parquet" ;;
    gnomad_ukb)    echo "gnomad.ukb.genetics_gym.per_variant.expected.parquet" ;;
    gnomad_v2)     echo "gnomad.v2.genetics_gym.per_variant.expected.parquet" ;;
    gnomad_v2_ukb) echo "gnomad.v2.ukb.genetics_gym.per_variant.expected.parquet" ;;
    *) return 1 ;;
  esac
}

# ---------------------------------------------------------------------------
# Options (env vars)
# ---------------------------------------------------------------------------
DRY_RUN="${DRY_RUN:-0}"     # 1 = print commands without executing (still uses gcloud ls for existence checks)
KEEP_RAW="${KEEP_RAW:-1}"   # 0 = rm the raw hail export after that cohort is uploaded
OVERWRITE="${OVERWRITE:-0}" # 1 = rebuild preprocessed output and force-reupload (delete GCS obj first)

# COHORTS env var can restrict the run to a subset (space-separated tags).
if [[ -n "${COHORTS:-}" ]]; then
  IFS=' ' read -r -a REQUESTED <<< "$COHORTS"
  for tag in "${REQUESTED[@]}"; do
    if ! raw_basename_for "$tag" >/dev/null 2>&1; then
      echo "ERROR: unknown cohort tag '$tag'." >&2
      echo "       Known tags: ${COHORT_TAGS[*]}" >&2
      exit 2
    fi
  done
  COHORT_TAGS=("${REQUESTED[@]}")
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
run() {
  # Echo + execute (or just echo if DRY_RUN=1). Preserves argv boundaries.
  if [[ "$DRY_RUN" == "1" ]]; then
    printf '    [dry-run]'
    printf ' %q' "$@"
    printf '\n'
    return 0
  fi
  "$@"
}

# gcloud storage cp -r prints one "Copying gs://..." line per file, which
# saturates any modestly-sized log buffer when a cohort has thousands of
# part-*.parquet files (~4k in this dataset). NB: gcloud writes that per-file
# progress to *stderr*, not stdout, so a plain `>/dev/null` does not silence
# it. --verbosity=error is the targeted knob: suppresses INFO-level "Copying"
# and "Average throughput" messages but leaves real errors on stderr.
quiet_gcloud_cp() {
  if [[ "$DRY_RUN" == "1" ]]; then
    run gcloud --verbosity=error storage cp "$@"
    return
  fi
  printf '    running (per-file progress silenced via --verbosity=error):'
  printf ' %q' gcloud --verbosity=error storage cp "$@"
  printf '\n'
  gcloud --verbosity=error storage cp "$@"
}

require() {
  # Check that a required command is on PATH before we start.
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: required command '$1' not found on PATH." >&2
    exit 3
  fi
}

exists_local() {
  # A downloaded partitioned parquet dir is considered "present" if it exists
  # AND contains at least one *.parquet part-file (guards against a half-copy).
  local dir="$1"
  [[ -d "$dir" ]] && compgen -G "$dir/*.parquet" >/dev/null
}

exists_gcs() {
  # Test whether a GCS object/prefix exists without depending on the return
  # code of `gcloud storage ls` printing on stderr.
  local uri="$1"
  gcloud storage ls "$uri" >/dev/null 2>&1
}

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
require gcloud
require compgen
if [[ ! -x "$VENV_PYTHON" ]]; then
  echo "ERROR: venv python not found at $VENV_PYTHON" >&2
  echo "       Activate/create it first, or edit VENV_PYTHON at the top of this script." >&2
  exit 4
fi
if [[ ! -f "$PREPROC_SCRIPT" ]]; then
  echo "ERROR: preprocessor script not found at $PREPROC_SCRIPT" >&2
  exit 5
fi

mkdir -p "$RAW_LOCAL_DIR" "$OUT_LOCAL_DIR"

echo "############################################################"
echo "### bootstrap_gnomad_obs_exp.sh"
echo "### cohorts:     ${COHORT_TAGS[*]}"
echo "### raw src:     $RAW_GCS_PREFIX"
echo "### raw dst:     $RAW_LOCAL_DIR"
echo "### out src:     $OUT_LOCAL_DIR"
echo "### out dst:     $OUT_GCS_PREFIX"
echo "### DRY_RUN=$DRY_RUN  KEEP_RAW=$KEEP_RAW  OVERWRITE=$OVERWRITE"
echo "############################################################"

# ---------------------------------------------------------------------------
# Per-cohort loop
# ---------------------------------------------------------------------------
for tag in "${COHORT_TAGS[@]}"; do
  raw_basename="$(raw_basename_for "$tag")"
  raw_gcs="$RAW_GCS_PREFIX/$raw_basename"
  raw_local="$RAW_LOCAL_DIR/$raw_basename"

  out_basename="${tag}_obs_exp_mis.parquet"
  out_local="$OUT_LOCAL_DIR/$out_basename"
  out_gcs="$OUT_GCS_PREFIX/$out_basename"

  echo
  echo "============================================================"
  echo "== cohort: $tag"
  echo "============================================================"

  # ---- [1/3] download raw hail export --------------------------------------
  echo ">>> [1/3] download raw hail export"
  if exists_local "$raw_local"; then
    echo "    raw already present at $raw_local; skipping download"
  else
    quiet_gcloud_cp -r "$raw_gcs" "$RAW_LOCAL_DIR/"
  fi
  # Belt-and-suspenders cleanup: these Hail exports carry over stale Spark
  # task-attempt files under _temporary/ that are zero-byte or partial (DuckDB
  # rejects them). The canonical part-*.parquet at top level supersedes them.
  # The _SUCCESS marker is likewise unnecessary once the download is complete.
  if [[ "$DRY_RUN" != "1" && -d "$raw_local/_temporary" ]]; then
    n_stale=$(find "$raw_local/_temporary" -type f 2>/dev/null | wc -l | tr -d ' ')
    echo "    cleaning $n_stale stale Spark _temporary file(s) from $raw_local/_temporary"
    rm -rf "$raw_local/_temporary"
  fi
  [[ "$DRY_RUN" != "1" && -f "$raw_local/_SUCCESS" ]] && rm -f "$raw_local/_SUCCESS"

  # ---- [2/3] preprocess ----------------------------------------------------
  echo ">>> [2/3] preprocess -> $out_basename"
  if [[ -f "$out_local" && "$OVERWRITE" != "1" ]]; then
    echo "    $out_local already exists; skipping preprocess (OVERWRITE=1 to rebuild)"
  else
    # Delegate to the per-cohort mode of the preprocessor so this loop stays
    # the source of truth for which cohorts run (avoids double-work if the
    # user restricts COHORTS).
    preproc_args=(
      "$PREPROC_SCRIPT"
      --input "$raw_local"
      --cohort "$tag"
      --output-dir "$OUT_LOCAL_DIR"
    )
    [[ "$OVERWRITE" == "1" ]] && preproc_args+=(--overwrite)
    run "$VENV_PYTHON" "${preproc_args[@]}"
  fi

  # ---- [3/3] upload preprocessed to eval_obs_exp/ --------------------------
  echo ">>> [3/3] upload -> $out_gcs"
  if [[ "$OVERWRITE" == "1" ]] && exists_gcs "$out_gcs"; then
    echo "    OVERWRITE=1: removing existing GCS object before reupload"
    run gcloud storage rm -r "$out_gcs"
  fi
  if exists_gcs "$out_gcs" && [[ "$OVERWRITE" != "1" ]]; then
    echo "    $out_gcs already in GCS; skipping upload (OVERWRITE=1 to force)"
  else
    quiet_gcloud_cp -r "$out_local" "$OUT_GCS_PREFIX/"
  fi

  # ---- optional cleanup ----------------------------------------------------
  if [[ "$KEEP_RAW" == "0" ]]; then
    echo ">>> KEEP_RAW=0: deleting raw scratch $raw_local"
    run rm -rf "$raw_local"
  fi
done

echo
echo "BOOTSTRAP_ALL_DONE"
echo
echo "Next: run the standard eval pipeline (download_source_data.py picks up"
echo "the new URIs from input_data_locations.json, then"
echo "create_variant_eval_all_table.py merges the 8 new obs/exp columns from"
echo "evals_input_data.json)."
