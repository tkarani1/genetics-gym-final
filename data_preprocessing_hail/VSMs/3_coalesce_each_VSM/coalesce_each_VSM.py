import hail as hl
import json
import os
hl.init(backend='spark', worker_memory="highmem", driver_memory='highmem') 
from resources.paths import WRITE_VSM_LINKER_TABLES_BASE, VSM_COUNTS_BASE, VSM_TABLE_NAMES, LINKER_PATHS, WRITE_SNP_COALESCED_VSM_BASE, FORMATTED_VSM_HT_PATHS
from resources.constants import SCORE_FIELDS
from resources.functions import write_tsv_bgz_from_ht, write_parquet_from_ht
from coalesce_each_VSM_functions import coalesce_VSM, coalesce_VSM_no_hierarchy

methods_all = [method for method in SCORE_FIELDS.keys()]

format_methods = ['CADD']
# methods = [x for x in methods_all if x not in format_methods]
# methods = ['ESM1B', 'CPT']

coalesce_type = 'SNP'
if coalesce_type == 'SNP':
    linker_ht_path = LINKER_PATHS['MISSENSE_ONLY_SNP']
    # linker_ht = hl.read_table(linker_ht_path)
    key_by = ['locus', 'alleles']

## Each column separately
# for method in methods:
#     for col in SCORE_FIELDS[method]:
#         col_dict = {col: SCORE_FIELDS[method][col]}
#         if method == 'AM':
#             ht_path = WRITE_VSM_LINKER_TABLES_BASE + VSM_TABLE_NAMES[method] + '_both_concat.ht'
#         ht_path = WRITE_VSM_LINKER_TABLES_BASE + VSM_TABLE_NAMES[method] + '.ht'
#         if not hl.hadoop_exists(ht_path):
#             print(f"Linker table not found: {method}")
#             continue
#         write_path_base = WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '_' + col
#         if hl.hadoop_exists(write_path_base + '.ht'):
#             print(f"Coalesced table already exists: {method}")
#             continue
#         ht = coalesce_VSM(ht_path, ['locus', 'alleles'], col_dict)
#         ht = ht.checkpoint(write_path_base + '.ht', overwrite=True)
#         # write_tsv_bgz_from_ht(ht, write_path_base + '.tsv.bgz')
#         write_parquet_from_ht(ht, write_path_base + '.parquet')
#         continue


## All columns together 
# for method in methods:    
#     ht_path = WRITE_VSM_LINKER_TABLES_BASE + VSM_TABLE_NAMES[method] + '.ht'
#     if not hl.hadoop_exists(ht_path):
#         print(f"Linker table not found: {method}")
#         continue

#     write_path = WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.ht'
#     if hl.hadoop_exists(write_path):
#         print(f"Coalesced table already exists: {method}")
#         continue

#     ht = coalesce_VSM(ht_path, ['locus', 'alleles'], SCORE_FIELDS[method])
#     ht = ht.checkpoint(WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.ht', overwrite=True)
#     score_cols = list(SCORE_FIELDS[method].keys())
#     # write_raw_and_collected_counts(ht, VSM_COUNTS_BASE, method, 'SNP collected', score_cols)
#     ht = locus_alleles_to_chr_pos_ref_alt(ht)
#     ht.export(WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.tsv.bgz', delimiter='\t')

# not joined to linker
for method in format_methods:
    ht_path = FORMATTED_VSM_HT_PATHS[method]
    print(ht_path)
    ht = coalesce_VSM_no_hierarchy(ht_path, linker_ht_path, ['locus', 'alleles'], SCORE_FIELDS[method])
    print(ht.describe())
    ht = ht.checkpoint(WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.ht', overwrite=True)
    print(ht.describe())
    write_parquet_from_ht(ht, WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.parquet')
    # write_tsv_bgz_from_ht(ht, WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.tsv.bgz')
