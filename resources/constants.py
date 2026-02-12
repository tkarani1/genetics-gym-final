ESM1B_SCORE_FIELD = "esm1b"
"""The ESM1B score field."""
MISFIT_D_SCORE_FIELD = "MisFit_D"
"""The MisFit D score field."""
MISFIT_S_SCORE_FIELD = "MisFit_S"
"""The MisFit S score field."""
POPEVE_SCORE_FIELD = "popEVE"
"""The popEVE score field."""
EVE_SCORE_FIELD = "EVE"
"""The EVE score field."""
ESM_1V_SCORE_FIELD = "ESM_1v"
"""The ESM-1v score field."""
MPC_SCORE_FIELD = "mpc"
"""The MPC score field."""
RASP_SCORE_FIELD = "rasp_score"
"""The RaSP score field."""
REVEL_SCORE_FIELD = "revel"
"""The REVEL score field."""
CPT1_SCORE_FIELD = "cpt1_score"
"""The CPT1 score field."""
PROTEINMPNN_LLR_SCORE_FIELD = "proteinmpnn_llr"
"""The ProteinMPNN LLR score field."""
PAI3D_SCORE_FIELD = "score_PAI3D"
"""The PAI3D score field."""
POLYPHEN_SCORE_FIELD = "polyphen_score"
"""The PolyPhen score field."""
CADD_SCORE_FIELD = "cadd_score"
"""The CADD score field."""
GPN_MSA_SCORE_FIELD = "gpn_msa_score"
"""The GPN-MSA score field."""
AM_SCORE_FIELD = "AM_score"
"""The AlphaMissense score field."""
AM_CANONICAL_SCORE_FIELD = "AM_score"
"""The AlphaMissense canonical score field."""

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
"""The score fields in each score resource."""

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

