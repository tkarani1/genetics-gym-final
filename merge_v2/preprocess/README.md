# `preprocess/`

**Out-of-pipeline** one-time source-normalization scripts and dated upload drivers.

These are **not** part of the recurring `merge_v2` build. `full_data_aggregation_pipeline.py`
does not invoke anything in this directory. Everything here is either:

- a **one-time bootstrap** that reshapes a raw upstream source into a form that
  the standard pipeline can consume (after which the finished asset is uploaded
  to GCS and referenced from `merge_v2/data_config/input_data_locations.json`
  like any other source), or
- a **dated upload / migration driver** for a specific rebuild.

Nothing in `merge_v2/src/` imports from this directory. If a preprocess step
here starts being needed on every build, promote it into `src/` and add it to
the orchestrator; until then it lives here so `src/` stays focused on the
recurring merge pipeline itself.

## Contents

### `preprocess_gnomad_obs_exp.py` + `bootstrap_gnomad_obs_exp.sh`

One-time bootstrap for the four gnomAD `*_obs_exp_mis.parquet` eval sources.
The bootstrap: downloads each raw Hail export → runs the preprocessor to
normalize it to the shape of the existing `*_obs_exp_mis.parquet` eval files
(canonical/MANE Ensembl missense, one row per variant, split `ref`/`alt`) →
uploads the result to `gs://grohlicek/.../eval_obs_exp/`.

After this bootstrap ran once, the URIs in `input_data_locations.json` under
`eval_tables.variant_level` already point at the uploaded GCS assets, so the
normal `download_source_data.py` + `create_variant_eval_all_table.py` flow
handles them like any other eval source. Kept here only for reproducibility.

```bash
bash merge_v2/preprocess/bootstrap_gnomad_obs_exp.sh                  # do everything
DRY_RUN=1  bash merge_v2/preprocess/bootstrap_gnomad_obs_exp.sh       # print only
OVERWRITE=1 bash merge_v2/preprocess/bootstrap_gnomad_obs_exp.sh      # rebuild + reupload
COHORTS="gnomad_new gnomad_v2" bash merge_v2/preprocess/bootstrap_gnomad_obs_exp.sh
```

### `preprocess_gnomad_v4_obs_exp.py`

Preprocesses the raw gnomAD v4 Hail export at

    gs://grohlicek/genetics_gym_vsm_all_content/parquet_files/obs_exp/gnomad_v4_obs_exp_per_variant.parquet

into the same `*_obs_exp_mis.parquet` shape consumed by
`create_variant_eval_all_table.py`. The raw source differs from the four
gnomAD cohorts handled by `preprocess_gnomad_obs_exp.py` in several ways:

- **No transcript-level columns** (`canonical`, `mane_select`, `annotation`,
  `transcript`), so the `WHERE` filters from `preprocess_gnomad_obs_exp.py`
  don't apply.
- **Dual `gene_id` encoding** -- every gene appears twice (ENSG + numeric
  Entrez/NCBI ID), a Hail VEP artefact.
- **Exact duplicate rows** within the ENSG subset.
- **77% non-missense variants** -- without filtering, the downstream tables
  would bloat from ~120M to ~500M+ rows.

The script normalizes via:

1. **INNER JOIN with `linker_all.parquet`** to restrict to missense variants.
2. `gene_id LIKE 'ENSG%'` filter (drops the Entrez duplicates).
3. `SELECT DISTINCT` (collapses exact-duplicate rows).
4. Derives `chrom`, `pos`, `ref`, `alt` from the Hail `locus.*` / `alleles[]`
   columns; renames `gene_id` -> `ensg`, `observed` -> `observed_gnomad_v4`,
   `expected` -> `expected_gnomad_v4`.

```bash
python preprocess_gnomad_v4_obs_exp.py                  # do everything
python preprocess_gnomad_v4_obs_exp.py --dry-run        # print SQL only
python preprocess_gnomad_v4_obs_exp.py --overwrite      # rebuild
python preprocess_gnomad_v4_obs_exp.py --memory-limit 8GB
```

Prerequisites:

* `merge_v2/data/raw_data/evals/gnomad_v4_obs_exp_per_variant.parquet/`
  (download the raw Hail export from GCS first)
* `merge_v2/data/raw_data/linker/linker_all.parquet` (fetched by the standard
  downloader)

Output: `merge_v2/data/raw_data/evals/gnomad_v4_obs_exp_mis.parquet`
(schema: `locus.contig, locus.position, alleles, enst, ensg,
observed_gnomad_v4, expected_gnomad_v4, chrom, pos, ref, alt`).

The processed file is already uploaded to
`gs://grohlicek/genetics_gym_vsm_all_content/eval_obs_exp/gnomad_v4_obs_exp_mis.parquet`
and registered in `input_data_locations.json`, so `download_source_data.py`
fetches the processed version directly. This script is kept for
reproducibility.

### `phase2_filtered_upload.sh` (+ `.log`)

Dated (2026-07-01) upload driver for a specific batch of gene-aggregated
`*_filtered` tables. Builds each `*_filtered.parquet` one at a time (smallest
first), uploads it to `gs://.../full_analysis_tables/2026_07_01_updated_tables/`,
and deletes the local copy so peak disk stays bounded. The `.log` is the
transcript of the actual run.

This is a one-off and will not be re-run in its current form. It's kept for
provenance and as a template for future incremental publishes.

### `msa_pairformer_chr22_sanity_check.py`

One-off diagnostic (chr22 only) that quantifies the MSA-Pairformer
UniProt-grain score against the missense/UniProt linker: coverage of chr22
UniProt accessions, per-`(uniprot_id, aa_pos, aa_ref, aa_alt)` fan-out, ~11%
isoform-disagreement rate on `aa_ref`, and per-policy match rate for any /
canonical-isoform / MANE-select / Ensembl-canonical reductions. The output is
a human-readable report only; no files are produced.

Run once, share the numbers with the score's author, decide the final isoform
policy, then hand off to the actual coalescer (see next entry).

### `preprocess_msa_pairformer_llr.py`

Coalesces the raw MSA-Pairformer LLR (uniprot grain,
`(uniprot_id, position, aa_ref, aa_alt)`, ~193.5M rows) down to variant grain
`(chrom, pos, ref, alt)` for consumption by
`create_variant_scores_all_table.py`. Uses the **canonical UniProt isoform**
policy: linker rows are kept only where `uniprot_isoform IS NULL OR
uniprot_isoform = uniprot_id` (i.e. the linker's designation of "canonical by
default"). Emits a single score column:

```
msa_pairformer_llr = max(llr)  -- canonical UniProt isoform only; max over
                                   any residual (uniprot_id, chrom, pos, ref, alt)
                                   duplicates surviving the filter
```

Earlier iterations of this preprocessor also emitted MANE-select and
any-isoform columns; those were retired because canonical UniProt is the
endorsed reduction. The three-column shape is in the git history if it is
ever needed for a comparison.

Engine: per-chromosome loop over
`../data/raw_data/linker/linker_missense_enst_transcript_aa_uniprot_by_chrom_tsv/
linker_missense_enst_transcript_aa_uniprot_{1..22,X,Y}.tsv.bgz`, one shard at
a time, each streamed through an INNER JOIN + canonical-isoform WHERE filter
+ GROUP BY into a per-chromosome Parquet checkpoint under `.ckpt_<output
stem>/`. A single UNION ALL consolidation then concatenates the 24
checkpoints into the final output. Peak memory is bounded to one shard's
hash-build side (~500K-3M rows); a crash loses at most one chromosome's
worth of work, and finished checkpoints are picked up on restart. Consistent
with the `.bgz` reader recipe in
`merge_v2/.cursor/rules/duckdb-file-reading.mdc` (`compression='gzip'`,
`nullstr=['NA','']`, explicit `types={...}` map to skip the ~50 s/shard
inference scan).

```bash
python preprocess_msa_pairformer_llr.py                  # do everything
python preprocess_msa_pairformer_llr.py --dry-run        # print the plan only
python preprocess_msa_pairformer_llr.py --overwrite      # rebuild from clean
python preprocess_msa_pairformer_llr.py --memory-limit 20GB --threads 4
```

Prerequisites (both fetched by the standard downloader):

* `merge_v2/data/raw_data/scores/msa_pairformer_llr_chunks0-9.parquet`
* `merge_v2/data/raw_data/linker/linker_missense_enst_transcript_aa_uniprot_by_chrom_tsv/`

Output: `merge_v2/data/raw_data/scores/msa_pairformer_variant.parquet`
(schema `chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR,
msa_pairformer_llr FLOAT`).

After the coalesced parquet is produced, upload it to GCS and register the
URI in `merge_v2/data_config/input_data_locations.json` (under
`score_tables.variant_level`) and the field projection in
`merge_v2/data_config/score_input_data.json` -- from that point on the
standard pipeline picks it up like any other variant-level score source.

Because the canonical filter covers only ~37% of the missense variants that
map to a UniProt accession (measured on the join output; the rest map only to
non-canonical isoforms), folding `msa_pairformer_llr` into the outer score
merge will substantially shrink the drop-any-NA `variant_scores_all_inner`
table. If that shrinkage is unacceptable for a given build, add
`msa_pairformer_variant.parquet` to the `--exclude` list of the
`create_variant_scores_all_table` step in `merge_v2/src/config.json` -- the
file stays on disk and registered in the manifest, but is skipped for the
merge (analogous to the four sparse-ish 2026-07 score sources already
excluded there). See the git history of `input_data_locations.json` /
`score_input_data.json` for the initial registration.
