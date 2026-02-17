"""Resources for processing missense scores."""

########################################################################################
# Missense score resource paths
########################################################################################
MISSENSE_SCORE_RESOURCE_PATHS = {
    "ESM1B": "gs://missense-scoring/esm1b-650m-brandes/proc/*.txt.bgz",
    "MISFIT": "gs://missense-scoring/misfit/raw/MisFit_by_Uniprot/*.txt.gz",
    "MISFIT_MAPPING": "gs://missense-scoring/misfit/raw/geneset_s_gene.txt",
    "POPEVE": "gs://missense-scoring/popEVE_ukbb_20250312/*.csv",
    "MPC": "gs://asc-v17/mpc_gnomad_2.1.1/mpc_grch38_deduped_with_outliers_2024-04-30.ht",
    "RASP": "gs://nnfc-fdp-konrad-public/RaSP/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.ht",
    "REVEL": "gs://missense-scoring/revel_with_transcript_ids",
    "CPT": "gs://missense-scoring/cpt_all/*.csv.gz",
    "PROTEINMPNN": "gs://missense-scoring/mpnn-outputs/AF_total_variants.ht",
    "PRIMATEAI3D": "gs://missense-scoring/primate_ai3d/PrimateAI-3D_scores.csv.bgz",
    "CONTEXT_VEP_ANNOTATED_PATH": "gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht",
    "CADD": "gs://genetics-gym-not-public/Trisha/whole_genome_SNVs.tsv.gz",
    "GPN_MSA": "gs://missense-scoring/GPN-MSA/scores.tsv.bgz",
    "AM_ISOFORMS": "gs://dm_alphamissense/AlphaMissense_isoforms_aa_substitutions.tsv.gz",
    "AM_CANONICAL": "gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz"
    }

AF2_UNIPROT_ISOFORM_MAPPING_PATH = (
    "gs://genetics-gym-not-public/Emily/af2db_uniprot_final.ht"
)
CPT_UNIPROT_ENTRY_TO_ID_MAPPING_PATH = "gs://genetics-gym/linkers/CPT_uniprot_entry_to_accession.tsv"
ENST_TO_UNIPROT_FOLDER = "gs://genetics-gym/linkers/ensembl_tid_to_uniprot_11_2025/"
ENST_TO_UNIPROT_HT = "gs://genetics-gym/linkers/enst_to_uniprot_mapping.ht"

########################################################################################
# Missense score resource paths
########################################################################################

EVERYTHING_RAW_HT_PATH = "gs://missense-scoring/mutation/everything_raw.ht"

# BASE_HT_PATH = "gs://gnomad-julia/genetics_gym/linkers/vsm_base_em_gene.ht"
# KEYED_HT_PATH = "gs://gnomad-tmp-4day/genetics_gym/linkers/vsm_base_em_gene_keyed.ht"
# PARTITIONS_HE_PATH = (
#     "gs://gnomad-julia/genetics_gym/linkers/vsm_base_em_gene_partitions.he"
# )

LINKER_PATHS = {
    "MISSENSE_ONLY_SNP": "gs://genetics-gym/linkers/linker_missense_only_snp.ht",
    "MISSENSE_ENST_TRANSCRIPT": "gs://genetics-gym/linkers/linker_missense_enst_transcript.ht",
    "MISSENSE_ENST_TRANSCRIPT_AA": "gs://genetics-gym/linkers/linker_missense_enst_transcript_aa.ht",
    "MISSENSE_ENST_TRANSCRIPT_AA_UNIPROT": "gs://genetics-gym/linkers/linker_missense_enst_transcript_aa_uniprot.ht",
    "MISSENSE_REFSEQ_TRNSCRIPT_AA": "gs://genetics-gym/linkers/linker_missense_refseq_transcript_aa.ht",
}

########################################################################################
# PATHS CREATED
########################################################################################

WRITE_VSM_TABLES_PATH = "gs://genetics-gym/vsm-tables/ht"

FORMATTED_VSM_HT_PATHS = {
    "AM_ISOFORMS": f"{WRITE_VSM_TABLES_PATH}/AM_isoforms.ht",
    "AM_CANONICAL": f"{WRITE_VSM_TABLES_PATH}/AM_canonical.ht",
    "CADD": f"{WRITE_VSM_TABLES_PATH}/cadd.ht",
    "CPT": f"{WRITE_VSM_TABLES_PATH}/cpt.ht",
    "ESM1B": f"{WRITE_VSM_TABLES_PATH}/esm1b.ht",
    "GPN_MSA": f"{WRITE_VSM_TABLES_PATH}/gpn_msa.ht",
    "MISFIT": f"{WRITE_VSM_TABLES_PATH}/misfit.ht",
    "MPC_PATH": f"{WRITE_VSM_TABLES_PATH}/mpc.ht",
    "PRIMATEAI3D": f"{WRITE_VSM_TABLES_PATH}/PAI3D.ht",
    "POPEVE": f"{WRITE_VSM_TABLES_PATH}/popeve.ht",
    "POLYPHEN": f"{WRITE_VSM_TABLES_PATH}/polyphen.ht",
    "PROTEINMPNN": f"{WRITE_VSM_TABLES_PATH}/protein_mpnn.ht",
    "RASP": f"{WRITE_VSM_TABLES_PATH}/rasp.ht",
    "REVEL": f"{WRITE_VSM_TABLES_PATH}/revel_enst.ht",
}
"""Formatted VSM HT output paths from data_preprocessing/VSMs/1_import_each_VSM_to_ht."""

WRITE_VSM_LINKER_TABLES_PATH = "gs://genetics-gym/vsm-tables/linker_and_vsm/linker_"
"""The path to write the VSM tables."""
VSM_LINKER_TABLES = {
    "AM": f"{WRITE_VSM_LINKER_TABLES_PATH}AM.ht",
    "CADD": f"{WRITE_VSM_LINKER_TABLES_PATH}cadd.ht",
    "CPT": f"{WRITE_VSM_LINKER_TABLES_PATH}cpt.ht",
    "ESM1B": f"{WRITE_VSM_LINKER_TABLES_PATH}esm1b.ht",
    "GPN_MSA": f"{WRITE_VSM_LINKER_TABLES_PATH}gpn_msa.ht",
    "MISFIT": f"{WRITE_VSM_LINKER_TABLES_PATH}misfit.ht",
    "MPC": f"{WRITE_VSM_LINKER_TABLES_PATH}mpc.ht",
    "PRIMATEAI3D": f"{WRITE_VSM_LINKER_TABLES_PATH}PAI3D.ht",
    "POPEVE": f"{WRITE_VSM_LINKER_TABLES_PATH}popeve.ht",
    "POLYPHEN": f"{WRITE_VSM_LINKER_TABLES_PATH}polyphen.ht",
    "PROTEINMPNN": f"{WRITE_VSM_LINKER_TABLES_PATH}protein_mpnn.ht",
    "RASP": f"{WRITE_VSM_LINKER_TABLES_PATH}rasp.ht",
    "REVEL": f"{WRITE_VSM_LINKER_TABLES_PATH}revel.ht",
}

VSM_COUNTS_FILE = f'{WRITE_VSM_LINKER_TABLES_PATH}/VSM_counts.json'

WRITE_SNP_COALESCED_VSM_TABLES_PATH = "gs://genetics-gym/vsm-tables/coalesced/coalesced_snp_"
COALESCED_VSM_TABLES = {
    "AM": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}AM.ht",
    "CADD": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}cadd.ht",
    "CPT": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}cpt.ht",
    "ESM1B": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}esm1b.ht",
    "GPN_MSA": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}gpn_msa.ht",
    "MISFIT": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}misfit.ht",
    "MPC": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}mpc.ht",
    "PRIMATEAI3D": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}PAI3D.ht",
    "POPEVE": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}popeve.ht",
    "POLYPHEN": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}polyphen.ht",
    "PROTEINMPNN": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}protein_mpnn.ht",
    "RASP": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}rasp.ht",
    "REVEL": f"{WRITE_SNP_COALESCED_VSM_TABLES_PATH}revel.ht",
}