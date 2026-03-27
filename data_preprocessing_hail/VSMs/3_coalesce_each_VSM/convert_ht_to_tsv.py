import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import WRITE_VSM_LINKER_TABLES_BASE, VSM_COUNTS_BASE, VSM_TABLE_NAMES, LINKER_PATHS, WRITE_SNP_COALESCED_VSM_BASE, FORMATTED_VSM_HT_PATHS
from resources.constants import SCORE_FIELDS
from resources.functions import write_tsv_bgz_from_ht_path, write_tsv_bgz_from_ht
from coalesce_each_VSM_functions import coalesce_VSM, coalesce_VSM_no_hierarchy

for method in VSM_TABLE_NAMES.keys():
    ht_path = WRITE_VSM_LINKER_TABLES_BASE + VSM_TABLE_NAMES[method] + '.ht'
    if not hl.hadoop_exists(ht_path):
        print(f"Table not found: {ht_path}, skipping method {method}")
        continue
    tsv_bgz_path = WRITE_VSM_LINKER_TABLES_BASE +  VSM_TABLE_NAMES[method] + '.tsv.bgz'
    write_tsv_bgz_from_ht_path(ht_path, tsv_bgz_path)