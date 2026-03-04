import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, WRITE_VSM_TABLES_PATH, FORMATTED_VSM_HT_PATHS
METHOD = 'RASP'

ht = hl.read_table(MISSENSE_SCORE_RESOURCE_PATHS[METHOD])
ht = ht.naive_coalesce(50)
ht = ht.annotate(
    uniprot_id = ht.uniprot,
    uniprot_aa_pos = hl.int(ht.variant[1:-1]),
    uniprot_aa_ref = ht.variant[0],
    uniprot_aa_alt = ht.variant[-1],
    rasp_score = ht.score_ml,
)
ht = ht.rename({'uniprot_aa_pos': 'aa_pos', 'uniprot_aa_ref': 'aa_ref', 'uniprot_aa_alt': 'aa_alt'})
ht = ht.key_by('uniprot_id')

af_ht = hl.read_table(MISSENSE_SCORE_RESOURCE_PATHS['AF2_UNIPROT_ISOFORM_MAPPING_PATH'])
ht = ht.join(af_ht, how='left')
ht = ht.key_by('uniprot_isoform', 'aa_pos', 'aa_ref', 'aa_alt')
ht = ht.select(ht.rasp_score)
ht.write(FORMATTED_VSM_HT_PATHS[METHOD])       
