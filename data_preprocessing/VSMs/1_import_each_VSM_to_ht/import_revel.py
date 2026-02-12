import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources import REVEL_PATH
from resources.paths import WRITE_VSM_TABLES_PATH
from resources.functions import locus_alleles_to_chr_pos_ref_alt

# Revel
revel_table = hl.import_table(REVEL_PATH, delimiter=',')
revel_table = revel_table.filter(~(revel_table.grch38_pos == '.'))
revel_table = revel_table.annotate(
    locus = hl.locus('chr' + revel_table.chr, hl.int(revel_table.grch38_pos), reference_genome='GRCh38'),
    alleles = [revel_table.ref, revel_table.alt], 
    aa_ref = revel_table.aaref, 
    aa_alt = revel_table.aaalt, 
    enst = revel_table.Ensembl_transcriptid.split(';'), 
    revel = hl.float(revel_table.REVEL)
)
revel_table = revel_table.explode('enst')
revel_table_fixed = revel_table.select('locus', 'alleles', 'enst', 'revel' )
revel_table_fixed.write(f'{WRITE_VSM_TABLES_PATH}/revel_enst.ht')
revel_table_fixed = locus_alleles_to_chr_pos_ref_alt(revel_table_fixed)
revel_table_fixed.export(f'{WRITE_VSM_TABLES_PATH}/revel_enst.tsv.bgz')