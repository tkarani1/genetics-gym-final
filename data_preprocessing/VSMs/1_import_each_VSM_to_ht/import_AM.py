import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, FORMATTED_VSM_HT_PATHS

# AM
AM_ISOFORMS_PATH = MISSENSE_SCORE_RESOURCE_PATHS['AM_ISOFORMS_AA_SUBSTITUTIONS']
AM_CANONICAL_PATH = MISSENSE_SCORE_RESOURCE_PATHS['AM_CANONICAL']

# all isoforms table
am_table = hl.import_table("gs://dm_alphamissense/AlphaMissense_isoforms_hg38.tsv.gz", 
                           delimiter='\t', force_bgz=True, no_header=False, comment='#')
am_table = am_table.annotate(
    locus = hl.locus(am_table.f0, hl.int(am_table.f1), reference_genome='GRCh38'), 
    alleles = [am_table.f2, am_table.f3], 
    AM_iso_score = hl.float(am_table.f7), 
    enst = am_table.f5.split('\.')[0]
    )
am_table = am_table.select('locus', 'alleles', 'enst', 'AM_iso_score')
am_table.write(FORMATTED_VSM_HT_PATHS['AM_ISOFORMS_PATH'])

# # canonical table
am_table = hl.import_table(AM_CANONICAL_PATH, 
                           delimiter='\t', force_bgz=True, no_header=True, comment='#')
am_table_fixed = am_table.annotate(
    locus = hl.locus(am_table.f0, hl.int(am_table.f1), reference_genome='GRCh38'),
    alleles = [am_table.f2, am_table.f3], 
    AM_score = hl.float(am_table.f8), 
    enst_orig = am_table.f6, 
    consequence = am_table.f9, 
    genome = am_table.f4
    )
am_table = am_table_fixed.annotate(enst = am_table_fixed.enst_orig.split('\.')[0])
am_table_fixed = am_table_fixed.drop('f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9')
am_table_fixed_k = am_table_fixed.key_by('locus', 'alleles', 'enst')
am_table_fixed_k.write(FORMATTED_VSM_HT_PATHS['AM_CANONICAL_PATH'])  

