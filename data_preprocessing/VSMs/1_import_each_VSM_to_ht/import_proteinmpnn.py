import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem')   
from resources import PROTEINMPNN_PATH, AF2DB_PATH
from resources.paths import WRITE_VSM_TABLES_PATH

## ProteinMPNN
mpnn_ht = hl.read_table(PROTEINMPNN_PATH)
mpnn_ht = mpnn_ht.key_by('uniprot')
# get isoform mapping from AF name
af_ht = hl.read_table(AF2DB_PATH)
mpnn_ht = mpnn_ht.join(af_ht, how='left')
mpnn_ht = mpnn_ht.key_by('uniprot_isoform', 'aa_pos', 'aa_ref', 'aa_alt')
mpnn_ht = mpnn_ht.select()
mpnn_ht.write(f'{WRITE_VSM_TABLES_PATH}/protein_mpnn.ht')
## make sure to join this twice later