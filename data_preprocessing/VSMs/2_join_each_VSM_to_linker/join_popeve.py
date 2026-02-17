import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import FORMATTED_VSM_HT_PATHS, WRITE_VSM_LINKER_TABLES_PATH
from resources.files import LINKER_PATHS

method = 'POPEVE'
linker = LINKER_PATHS['MISSENSE_REFSEQ_TRNSCRIPT_AA']


linker_ht = hl.read_table(linker)
linker_ht = linker_ht.key_by('NP', 'aa_pos', 'aa_ref', 'aa_alt')

vsm_ht = hl.read_table(FORMATTED_VSM_HT_PATHS[method])
vsm_ht = vsm_ht.key_by('NP', 'aa_pos', 'aa_ref', 'aa_alt')

ht = linker_ht.join(vsm_ht, how='right')
ht = ht.checkpoint(f'gs://trisha-tmp/new_VSM_temp/vsm_{method}.ht', overwrite=True)
