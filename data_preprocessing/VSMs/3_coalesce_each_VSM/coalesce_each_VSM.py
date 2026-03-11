import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import WRITE_VSM_LINKER_TABLES_BASE, VSM_COUNTS_BASE, VSM_TABLE_NAMES, LINKER_PATHS, WRITE_SNP_COALESCED_VSM_BASE, FORMATTED_VSM_HT_PATHS
from resources.constants import SCORE_FIELDS
from resources.functions import write_raw_and_collected_counts, locus_alleles_to_chr_pos_ref_alt
from coalesce_each_VSM_functions import coalesce_VSM, coalesce_VSM_no_hierarchy

methods = [method for method in SCORE_FIELDS.keys()]
format_methods = ['CADD', 'GPN_MSA']
methods.drop(format_methods)

coalesce_type = 'SNP'
if coalesce_type == 'SNP':
    linker_ht = hl.read_table(LINKER_PATHS['MISSENSE_ONLY_SNP'])
    key_by = ['locus', 'alleles']

for method in methods:
    ht_path = WRITE_VSM_LINKER_TABLES_BASE + VSM_TABLE_NAMES[method] + '.ht'
    ht = coalesce_VSM(ht_path, ['locus', 'alleles'], SCORE_FIELDS[method])
    ht = ht.checkpoint(WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.ht', overwrite=True)
    score_cols = list(SCORE_FIELDS[method].keys())
    write_raw_and_collected_counts(ht, VSM_COUNTS_BASE, method, 'SNP collected', score_cols)
    ht = locus_alleles_to_chr_pos_ref_alt(ht)
    ht.export(WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.tsv.bgz', delimiter='\t')

for method in format_methods:
    ht_path = FORMATTED_VSM_HT_PATHS[method]
    ht = coalesce_VSM_no_hierarchy(ht_path, ['locus', 'alleles'], SCORE_FIELDS[method])
    ht = ht.checkpoint(WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.ht', overwrite=True)
    score_cols = list(SCORE_FIELDS[method].keys())
    write_raw_and_collected_counts(ht, VSM_COUNTS_BASE, method, 'SNP collected', score_cols)
    ht = locus_alleles_to_chr_pos_ref_alt(ht)
    ht.export(WRITE_SNP_COALESCED_VSM_BASE + VSM_TABLE_NAMES[method] + '.tsv.bgz', delimiter='\t')


