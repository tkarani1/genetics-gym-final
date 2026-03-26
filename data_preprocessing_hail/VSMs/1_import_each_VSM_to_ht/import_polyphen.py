import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, FORMATTED_VSM_HT_PATHS

METHOD = 'POLYPHEN'
## Polyphen2
ht = hl.read_table(MISSENSE_SCORE_RESOURCE_PATHS[METHOD])
ht = ht.explode(ht.vep.transcript_consequences)
ht = ht.filter(ht.vep.transcript_consequences.transcript_id.startswith('ENST') )
ht = ht.annotate(
    polyphen_score = ht.vep.transcript_consequences.polyphen_score, 
    enst = ht.vep.transcript_consequences.transcript_id, 
    ensg = ht.vep.transcript_consequences.gene_id, 
    gene_symbol = ht.ranscript_consequences.gene_symbol)
ht = ht.select(ht.enst, ht.polyphen_score, ht.ensg, ht.gene_symbol)
ht = ht.filter(hl.is_defined(ht.polyphen_score))
ht.write(FORMATTED_VSM_HT_PATHS[METHOD])    