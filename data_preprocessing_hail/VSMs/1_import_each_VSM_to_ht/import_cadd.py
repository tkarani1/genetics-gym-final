import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, FORMATTED_VSM_HT_PATHS

METHOD = 'CADD'
ht = hl.import_table(MISSENSE_SCORE_RESOURCE_PATHS[METHOD], 
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
ht.write(FORMATTED_VSM_HT_PATHS[METHOD])    