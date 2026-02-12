import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources import MPC_PATH
from resources.paths import WRITE_VSM_TABLES_PATH
from resources.functions import locus_alleles_to_chr_pos_ref_alt

# # MPC
# # original: 
ht = hl.read_table(MPC_PATH)
ht = ht.key_by('locus','alleles')
ht.write(f'{WRITE_VSM_TABLES_PATH}/mpc_grch38_deduped_with_outliers_2024-04-30.ht')
ht = locus_alleles_to_chr_pos_ref_alt(ht)
ht.export(f'{WRITE_VSM_TABLES_PATH}/mpc_grch38_deduped_with_outliers_2024-04-30.tsv.bgz')

