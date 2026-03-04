import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem')   
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, WRITE_VSM_TABLES_PATH, FORMATTED_VSM_HT_PATHS
METHOD = 'PROTEINMPNN'

## ProteinMPNN
ht = hl.read_table(MISSENSE_SCORE_RESOURCE_PATHS[METHOD])
ht = ht.key_by('uniprot')
# get isoform mapping from AF name
af_ht = hl.read_table(MISSENSE_SCORE_RESOURCE_PATHS['AF2_UNIPROT_ISOFORM_MAPPING_PATH'])
ht = ht.join(af_ht, how='left')
ht = ht.key_by('uniprot_isoform', 'aa_pos', 'aa_ref', 'aa_alt')
ht = ht.select(ht.proteinmpnn_llr)
ht.write(FORMATTED_VSM_HT_PATHS[METHOD])    
