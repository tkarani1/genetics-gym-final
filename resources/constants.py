## TODO: update naming conventions 

ESM1B_SCORE_FIELD = "esm1b"
MISFIT_D_SCORE_FIELD = "MisFit_D"
MISFIT_S_SCORE_FIELD = "MisFit_S"
POPEVE_SCORE_FIELD = "popEVE"
EVE_SCORE_FIELD = "EVE"
ESM_1V_SCORE_FIELD = "ESM_1v"
MPC_SCORE_FIELD = "mpc"
RASP_SCORE_FIELD = "rasp_score"
REVEL_SCORE_FIELD = "revel"
CPT1_SCORE_FIELD = "cpt1_score"
PROTEINMPNN_LLR_SCORE_FIELD = "proteinmpnn_llr"
PAI3D_SCORE_FIELD = "score_PAI3D"
POLYPHEN_SCORE_FIELD = "polyphen_score"
CADD_SCORE_FIELD = "cadd_score"
GPN_MSA_SCORE_FIELD = "gpn_msa_score"
AM_SCORE_FIELD = "AM_score"
AM_CANONICAL_SCORE_FIELD = "AM_score"

SCORE_FIELDS = {
    "esm1b": [ESM1B_SCORE_FIELD],
    "misfit": [MISFIT_D_SCORE_FIELD, MISFIT_S_SCORE_FIELD],
    "popeve": [POPEVE_SCORE_FIELD, EVE_SCORE_FIELD, ESM_1V_SCORE_FIELD],
    "mpc": [MPC_SCORE_FIELD],
    "rasp": [RASP_SCORE_FIELD],
    "revel": [REVEL_SCORE_FIELD],
    "cpt": [CPT1_SCORE_FIELD],
    "proteinmpnn": [PROTEINMPNN_LLR_SCORE_FIELD],
    "pai3d": [PAI3D_SCORE_FIELD],
    "polyphen": [POLYPHEN_SCORE_FIELD],
    "cadd": [CADD_SCORE_FIELD],
    "gpn_msa": [GPN_MSA_SCORE_FIELD],
    "am_isos": [AM_SCORE_FIELD],
    "am_canonical": [AM_SCORE_FIELD],
}
"""The score fields in each score table."""

HIGHER_IS_LESS_DELETERIOUS = {
    "esm1b": {ESM1B_SCORE_FIELD: True},
    "misfit": {MISFIT_D_SCORE_FIELD: False, MISFIT_S_SCORE_FIELD: False},
    "popeve": {
        POPEVE_SCORE_FIELD: True,
        EVE_SCORE_FIELD: False,
        ESM_1V_SCORE_FIELD: True,
    },
    "mpc": {MPC_SCORE_FIELD: False},
    "rasp": {RASP_SCORE_FIELD: False},
    "revel": {REVEL_SCORE_FIELD: False},
    "cpt": {CPT1_SCORE_FIELD: False},
    "proteinmpnn": {PROTEINMPNN_LLR_SCORE_FIELD: True},
    "pai3d": {PAI3D_SCORE_FIELD: False},
    "polyphen": {POLYPHEN_SCORE_FIELD: False},
    "cadd": {CADD_SCORE_FIELD: False},
    "gpn_msa": {GPN_MSA_SCORE_FIELD: True},
    "am_isos": {AM_SCORE_FIELD: False},
    "am_canonical": {AM_SCORE_FIELD: False},
}

