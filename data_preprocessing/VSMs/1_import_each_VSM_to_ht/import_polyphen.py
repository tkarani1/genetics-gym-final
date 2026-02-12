import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources import POLYPHE2_PATH
from resources.paths import WRITE_VSM_TABLES_PATH


## Polyphen2
ht = hl.read_table(POLYPHE2_PATH)
ht = ht.explode(ht.vep.transcript_consequences)
ht = ht.filter(ht.vep.transcript_consequences.transcript_id.startswith('ENST') )
ht = ht.annotate(
    polyphen_score = ht.vep.transcript_consequences.polyphen_score, 
    enst = ht.vep.transcript_consequences.transcript_id, 
    ensg = ht.vep.transcript_consequences.gene_id, 
    gene_symbol = ht.ranscript_consequences.gene_symbol)
ht = ht.select(ht.enst, ht.polyphen_score, ht.ensg, ht.gene_symbol)
ht = ht.filter(hl.is_defined(ht.polyphen_score))
ht.write(f'{WRITE_VSM_TABLES_PATH}/polyphen.ht')