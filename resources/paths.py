WRITE_VSM_TABLES_PATH = "gs://genetics-gym/vsm-tables"
"""The path to write the VSM tables."""

FORMATTED_VSM_HT_PATHS = {
    "AM_ISOFORMS_PATH": f"{WRITE_VSM_TABLES_PATH}/AM_isoforms.ht",
    "AM_CANONICAL_PATH": f"{WRITE_VSM_TABLES_PATH}/AM_canonical.ht",
    "CADD_PATH": f"{WRITE_VSM_TABLES_PATH}/cadd.ht",
    "ESM1B_PATH": f"{WRITE_VSM_TABLES_PATH}/esm1b_scores.ht",
    "GPN_MSA_PATH": f"{WRITE_VSM_TABLES_PATH}/gpn_msa.ht",
    "MISFIT_PATH": f"{WRITE_VSM_TABLES_PATH}/misfit_scores.ht",
    "MPC_PATH": f"{WRITE_VSM_TABLES_PATH}/mpc_grch38_deduped_with_outliers_2024-04-30.ht",
    "PRIMATEAI3D_PATH": f"{WRITE_VSM_TABLES_PATH}/PAI3D.ht",
    "POPEVE_PATH": f"{WRITE_VSM_TABLES_PATH}/popeve_2025.ht",
    "POLYPHEN_PATH": f"{WRITE_VSM_TABLES_PATH}/polyphen.ht",
    "PROTEINMPNN_PATH": f"{WRITE_VSM_TABLES_PATH}/protein_mpnn.ht",
    "RASP_PATH": f"{WRITE_VSM_TABLES_PATH}/rasp_scores.ht",
    "REVEL_PATH": f"{WRITE_VSM_TABLES_PATH}/revel_enst.ht",
}
"""Formatted VSM HT output paths from data_preprocessing/VSMs/1_import_each_VSM_to_ht."""