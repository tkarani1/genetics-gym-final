import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, FORMATTED_VSM_HT_PATHS

METHOD = 'CPT'

cpt_ht = hl.import_table(MISSENSE_SCORE_RESOURCE_PATHS[METHOD], 
                        delimiter=',', force=True, source_file_field='filename')
cpt_ht = cpt_ht.annotate(
    aa_pos = hl.int(cpt_ht.mutant[1:-1]),
    aa_ref = cpt_ht.mutant[0],
    aa_alt = cpt_ht.mutant[-1],
    cpt_score = hl.float(cpt_ht.CPT1_score), 
    uniprot_entry = cpt_ht.filename.split('/')[-1].split('\.')[0]
)
mapping_ht = hl.import_table(MISSENSE_SCORE_RESOURCE_PATHS['CPT_UNIPROT_ENTRY_TO_ID_MAPPING_PATH'], 
                            delimiter='\t')
mapping_ht = mapping_ht.key_by('entry_name')
cpt_ht = cpt_ht.join(mapping_ht, how='left')

cpt_ht = cpt_ht.rename({'accession': 'uniprot_id'})
cpt_ht.write(FORMATTED_VSM_HT_PATHS[METHOD])    
