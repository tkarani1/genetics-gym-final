import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, WRITE_VSM_TABLES_PATH, FORMATTED_VSM_HT_PATHS

METHOD = 'MPC'
# # MPC
# # original: 
ht = hl.read_table(MISSENSE_SCORE_RESOURCE_PATHS[METHOD])
ht = ht.rename({'mpc': 'mpc_score', 'transcript': 'enst'})
ht = ht.key_by('locus','alleles', 'enst')
ht.write(FORMATTED_VSM_HT_PATHS[METHOD], overwrite=True)    


