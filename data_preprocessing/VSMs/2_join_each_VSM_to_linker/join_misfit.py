import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import FORMATTED_VSM_HT_PATHS, WRITE_VSM_LINKER_TABLES_PATH
from resources.files import LINKER_PATHS

METHOD = 'MISFIT'
linker = LINKER_PATHS['MISSENSE_ENST_TRANSCRIPT_AA']


# load linker
linker_ht = hl.read_table(linker)
linker_ht = linker_ht.key_by('enst', 'aa_pos', 'aa_alt')

# load method ht
vsm_ht = hl.read_table(FORMATTED_VSM_HT_PATHS[METHOD])
vsm_ht = vsm_ht.key_by('enst', 'aa_pos', 'aa_alt')

# join
ht = linker_ht.join(vsm_ht, how='right')

ht = ht.checkpoint(WRITE_VSM_LINKER_TABLES_PATH[METHOD], overwrite=True)