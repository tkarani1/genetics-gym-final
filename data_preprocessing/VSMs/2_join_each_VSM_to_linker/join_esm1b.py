import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import FORMATTED_VSM_HT_PATHS, WRITE_VSM_LINKER_TABLES_PATH
from resources.files import LINKER_PATHS

method = 'ESM1B'
linker  = LINKER_PATHS['MISSENSE_ENST_TRANSCRIPT_AA_UNIPROT']

linker_ht = hl.read_table(linker)
linker_ht = linker_ht.key_by('uniprot_isoform', 'aa_pos', 'aa_ref', 'aa_alt')

vsm_ht = hl.read_table(FORMATTED_VSM_HT_PATHS[method])
vsm_ht = vsm_ht.key_by('uniprot_isoform', 'aa_pos', 'aa_ref', 'aa_alt')

ht = linker_ht.join(vsm_ht, how='right')

ht = ht.key_by('uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt')
ht = ht.join(vsm_ht, how='right')

ht = ht.annotate(
    esm1b = hl.coalesce(ht.esm_score, ht.esm_score_1)
)
ht = ht.checkpoint(WRITE_VSM_LINKER_TABLES_PATH[method], overwrite=True)