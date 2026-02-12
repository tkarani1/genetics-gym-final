import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources import CADD_PATH
from resources.paths import WRITE_VSM_TABLES_PATH

ht = hl.import_table(CADD_PATH, 
            filter = '##', delimiter = '\t', force_bgz =True)
ht = ht.select(
    chrom = ht['#Chrom'], 
    pos = hl.int(ht['Pos']),
    ref = ht['Ref'], 
    alt = ht['Alt'], 
    raw_score = hl.float(ht['RawScore']), 
    
)
ht = ht.select(
    locus = hl.locus('chr' + ht['chrom'],ht['pos'], reference_genome='GRCh38'),
    alleles = [ht['ref'], ht['alt']],
    cadd_score = ht['raw_score']
)
ht = ht.key_by('locus', 'alleles')
ht = ht.filter(hl.is_defined(ht.cadd_score))
ht.write(f'{WRITE_VSM_TABLES_PATH}/cadd.ht')