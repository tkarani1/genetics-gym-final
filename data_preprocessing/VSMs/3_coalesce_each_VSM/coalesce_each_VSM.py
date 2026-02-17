import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import VSM_LINKER_TABLES, LINKER_PATHS, WRITE_SNP_COALESCED_VSM_TABLES_PATH
from resources.constants import HIGHER_IS_LESS_DELETERIOUS
from coalesce_VSM_functions import coalesce_VSM_by_SNP

methods = [method for method in VSM_LINKER_TABLES.keys()]


coalesce_type = 'SNP'
if coalesce_type == 'SNP':
    methods = methods.drop('CADD', 'GPN_MSA', 'MPC')    
    for method in methods:
        ht_path = VSM_LINKER_TABLES[method]
        ht, final_scores = coalesce_VSM_by_SNP(ht_path, ['locus', 'alleles'], HIGHER_IS_LESS_DELETERIOUS[method])
        ht.write(WRITE_SNP_COALESCED_VSM_TABLES_PATH[method])
        print(final_scores)


