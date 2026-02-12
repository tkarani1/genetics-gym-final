"""Resources for processing missense scores."""

########################################################################################
# Missense score resource paths
########################################################################################

ESM1B_PATH = "gs://missense-scoring/esm1b-650m-brandes/proc/*.txt.bgz"
MISFIT_PATH = "gs://missense-scoring/misfit/raw/MisFit_by_Uniprot/*.txt.gz"
MISFIT_MAPPING_PATH = "gs://missense-scoring/misfit/raw/geneset_s_gene.txt"
POPEVE_PATH = "gs://missense-scoring/popEVE_ukbb_20250312/*.csv"
MPC_PATH = (
    "gs://asc-v17/mpc_gnomad_2.1.1/mpc_grch38_deduped_with_outliers_2024-04-30.ht"
)
RASP_PATH = (
    "gs://nnfc-fdp-konrad-public/RaSP/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.ht"
)
REVEL_PATH = "gs://missense-scoring/revel_with_transcript_ids"
CPT_PATH = "gs://missense-scoring/cpt_all/*.csv.gz"
PROTEINMPNN_PATH = "gs://missense-scoring/mpnn-outputs/AF_total_variants.ht"
PRIMATEAI3D_PATH = "gs://missense-scoring/primate_ai3d/PrimateAI-3D_scores.csv.bgz"
CONTEXT_VEP_ANNOTATED_PATH = "gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht"
CADD_PATH = "gs://genetics-gym-not-public/Trisha/whole_genome_SNVs.tsv.gz"
GPN_MSA_PATH = "gs://missense-scoring/GPN-MSA/scores.tsv.bgz"
AM_ISOFORMS_PATH = (
    "gs://dm_alphamissense/AlphaMissense_isoforms_aa_substitutions.tsv.gz"
)
AM_CANONICAL_PATH = "gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz"
AF2_UNIPROT_ISOFORM_MAPPING_PATH = (
    "gs://genetics-gym-not-public/Emily/af2db_uniprot_final.ht"
)

########################################################################################
# Missense score resource paths
########################################################################################

EVERYTHING_RAW_HT_PATH = "gs://missense-scoring/mutation/everything_raw.ht"

BASE_HT_PATH = "gs://gnomad-julia/genetics_gym/linkers/vsm_base_em_gene.ht"
KEYED_HT_PATH = "gs://gnomad-tmp-4day/genetics_gym/linkers/vsm_base_em_gene_keyed.ht"
PARTITIONS_HE_PATH = (
    "gs://gnomad-julia/genetics_gym/linkers/vsm_base_em_gene_partitions.he"
)