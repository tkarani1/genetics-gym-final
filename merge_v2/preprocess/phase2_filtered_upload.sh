#!/usr/bin/env bash
# Phase 2: build the 4 gene-aggregated *_filtered tables one at a time (smallest
# first), uploading each _eval + _filtered to GCS and deleting the local copies
# as we go so peak local disk stays well under the free-space budget.
set -euo pipefail

SRC="/Users/grohlice/Source/genetics-gym-final/merge_v2/src"
GENE_DIR="/Users/grohlice/Source/genetics-gym-final/merge_v2/data/processed_data/full_analysis_tables/gene_aggregated"
VENV="/Users/grohlice/Source/genetics-gym-final/merge_v2/.venv/bin/python"
DEST="gs://grohlicek/genetics_gym_vsm_all_content/full_analysis_tables/2026_07_01_updated_tables"

# Smallest -> largest so early deletions free room for the big outer builds.
STEMS=(
  "variant_scores_all_inner_ensg_stats_eval"
  "variant_scores_all_inner_ensg_stats_eval_deduped"
  "variant_scores_all_outer_ensg_stats_eval_deduped"
  "variant_scores_all_outer_ensg_stats_eval"
)

cd "$SRC"

for stem in "${STEMS[@]}"; do
  eval_path="$GENE_DIR/$stem.parquet"
  filt_path="$GENE_DIR/${stem}_filtered.parquet"

  echo "############################################################"
  echo "### TABLE: $stem"
  echo "### free before: $(df -h "$GENE_DIR" | tail -1 | awk '{print $4}')"
  echo "############################################################"

  echo ">>> [1/4] upload _eval -> $DEST/$stem.parquet"
  if gcloud storage ls "$DEST/$stem.parquet" >/dev/null 2>&1; then
    echo "    already in GCS; skipping _eval upload"
  else
    gcloud storage cp "$eval_path" "$DEST/$stem.parquet"
  fi

  echo ">>> [2/4] build _filtered"
  # Fewer threads + smaller row groups keep the (non-spillable) Parquet write
  # buffer small: the final ~340-col batch OOM'd at 512k row groups x 6 threads.
  "$VENV" create_gene_aggregated_filtered_table.py \
    --input "$eval_path" --memory-limit 20GB --threads 3 \
    --row-group-size 100000 --overwrite

  echo ">>> [2.5] verify new score col + a filter col actually populated"
  "$VENV" - "$filt_path" <<'PY'
import sys, duckdb
p = sys.argv[1]
con = duckdb.connect()
cols = {r[0] for r in con.execute(
    f"describe select * from read_parquet('{p}')").fetchall()}
assert "MutScore" in cols, "MutScore score column missing from filtered table"
assert "buried" in cols, "buried filter column missing from filtered table"
n_score = con.execute(
    f"select count(MutScore) from read_parquet('{p}')").fetchone()[0]
n_true = con.execute(
    f"select count(*) from read_parquet('{p}') where buried").fetchone()[0]
print(f"    verify: non-null MutScore={n_score:,}  buried=TRUE rows={n_true:,}")
assert n_score > 0, "MutScore is entirely NULL -- score merge/carry-through broke"
assert n_true > 0, "buried filter is entirely FALSE -- filter key join mismatch"
PY

  echo ">>> [3/4] upload _filtered -> $DEST/${stem}_filtered.parquet"
  gcloud storage cp "$filt_path" "$DEST/${stem}_filtered.parquet"

  echo ">>> [4/4] delete local _eval + _filtered"
  rm -f "$eval_path" "$filt_path"

  echo "### free after: $(df -h "$GENE_DIR" | tail -1 | awk '{print $4}')"
  echo
done

echo "PHASE2_ALL_DONE"
