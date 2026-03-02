from gc import collect
import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import FORMATTED_VSM_HT_PATHS, WRITE_VSM_LINKER_TABLES_PATH, LINKER_PATHS, WRITE_VSM_LINKER_TABLES_BASE, VSM_TABLE_NAMES
from resources.files import LINKER_PATHS
from resources.functions import write_raw_and_collected_counts


METHOD = 'CPT'
linker = LINKER_PATHS['MISSENSE_ENST_TRANSCRIPT_AA_UNIPROT']

# load linker
linker_ht = hl.read_table(linker)
linker_ht = linker_ht.key_by('uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt')

# load method ht
vsm_ht = hl.read_table(FORMATTED_VSM_HT_PATHS[METHOD])
vsm_ht = vsm_ht.key_by('uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt')

ht = linker_ht.join(vsm_ht, how='right')
ht = ht.checkpoint(WRITE_VSM_LINKER_TABLES_PATH[METHOD], overwrite=True)

write_raw_and_collected_counts(vsm_ht, WRITE_VSM_LINKER_TABLES_PATH[METHOD], METHOD, 'uniprot_id')