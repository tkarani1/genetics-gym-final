import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, WRITE_VSM_TABLES_PATH, FORMATTED_VSM_HT_PATHS

METHOD = 'GPN_MSA'
## GPN-MSA
ht = hl.import_table(MISSENSE_SCORE_RESOURCE_PATHS[METHOD], delimiter='\t', force_bgz=True, no_header=True)
ht = ht.annotate(
    chr = 'chr'+ht.f0, 
    pos = hl.int(ht.f1), 
    alleles=[ht.f2, ht.f3], 
    gpn_msa_score = hl.float(ht.f4)
    )
ht = ht.annotate(
    locus = hl.locus(ht.chr, ht.pos, reference_genome='GRCh38'),
)
ht = ht.select(ht.locus, ht.alleles, ht.gpn_msa_score)
ht = ht.key_by('locus', 'alleles')
ht.write(FORMATTED_VSM_HT_PATHS[METHOD])    

