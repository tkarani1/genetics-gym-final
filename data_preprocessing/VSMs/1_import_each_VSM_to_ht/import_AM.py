import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, FORMATTED_VSM_HT_PATHS

# AM
AM_ISOFORMS_PATH = MISSENSE_SCORE_RESOURCE_PATHS['AM_ISOFORMS']
AM_CANONICAL_PATH = MISSENSE_SCORE_RESOURCE_PATHS['AM_CANONICAL']
# all isoforms table
am_table = hl.import_table(AM_ISOFORMS_PATH, 
                           delimiter='\t', force_bgz=True, no_header=False, comment='#')
am_table = am_table.annotate(
   enst = am_table.transcript_id.split('\.')[0], 
   aa_pos = hl.int(am_table.protein_variant[1:-1]),
   aa_ref = am_table.protein_variant[0],
   aa_alt = am_table.protein_variant[-1],
   AM_iso = hl.float(am_table.am_pathogenicity), 
)
am_table = am_table.select('enst', 'aa_pos', 'aa_ref', 'aa_alt', 'AM_iso')
am_table.write(FORMATTED_VSM_HT_PATHS['AM_ISOFORMS'], overwrite=True)

# # canonical table
am_table = hl.import_table(AM_CANONICAL_PATH, 
                           delimiter='\t', force_bgz=True, no_header=True, comment='#')
am_table = am_table.annotate(
    locus = hl.locus(am_table.f0, hl.int(am_table.f1), reference_genome='GRCh38'),
    alleles = [am_table.f2, am_table.f3], 
    enst = am_table.f6.split('\.')[0], 
    AM_canon = hl.float(am_table.f8), 
    )
am_table = am_table.select('locus', 'alleles', 'enst', 'AM_canon')
am_table.write(FORMATTED_VSM_HT_PATHS['AM_CANONICAL'], overwrite=True)  

