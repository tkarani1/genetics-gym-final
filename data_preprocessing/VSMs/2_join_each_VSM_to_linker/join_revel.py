import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import FORMATTED_VSM_HT_PATHS, WRITE_VSM_LINKER_TABLES_BASE, VSM_TABLE_NAMES, VSM_COUNTS_BASE, LINKER_PATHS
from resources.functions import write_raw_and_collected_counts

METHOD = 'REVEL'
linker = LINKER_PATHS['MISSENSE_ENST_TRANSCRIPT']

linker_ht = hl.read_table(linker)
linker_ht = linker_ht.key_by('locus', 'alleles', 'enst')

vsm_ht = hl.read_table(FORMATTED_VSM_HT_PATHS[METHOD])
vsm_ht = vsm_ht.key_by('locus', 'alleles', 'enst')

ht = linker_ht.join(vsm_ht, how='right')
ht = ht.checkpoint(WRITE_VSM_LINKER_TABLES_BASE + VSM_TABLE_NAMES[METHOD] + '.ht', overwrite=True)

write_raw_and_collected_counts(vsm_ht, f'{VSM_COUNTS_BASE}/{VSM_TABLE_NAMES[METHOD]}.json', METHOD, 'SNP - enst vsm table')
write_raw_and_collected_counts(ht, f'{VSM_COUNTS_BASE}/{VSM_TABLE_NAMES[METHOD]}.json', METHOD, 'SNP - enst linker table')