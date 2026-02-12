import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources import AM_ISOFORMS_PATH, AM_CANONICAL_PATH
from resources.paths import WRITE_VSM_TABLES_PATH

# AM
# all isoforms table
am_table = hl.import_table(AM_ISOFORMS_PATH, 
                           delimiter='\t', force_bgz=True, no_header=False, comment='#')
am_table = am_table.annotate(
   enst = am_table.transcript_id.split('\.')[0], 
   aa_pos = hl.int(am_table.protein_variant[1:-1]),
   aa_ref = am_table.protein_variant[0],
   aa_alt = am_table.protein_variant[-1],
   AM = hl.float(am_table.am_pathogenicity), 
   enst_orig = am_table.transcript_id
)
am_table = am_table.select('enst', 'aa_pos', 'aa_ref', 'aa_alt', 'AM', 'enst_orig')
am_table.write(f'{WRITE_VSM_TABLES_PATH}/AM_isoforms.ht')

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
am_table_fixed_k.write(f'{WRITE_VSM_TABLES_PATH}/AM_canonical.ht')

