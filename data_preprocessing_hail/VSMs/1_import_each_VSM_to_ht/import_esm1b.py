import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, FORMATTED_VSM_HT_PATHS

METHOD = 'ESM1B'

# # ESM1b
esm1b_ht = hl.import_table(
        MISSENSE_SCORE_RESOURCE_PATHS[METHOD],
        types = {
            'aa_pos':hl.tint,
            'aa_ref':hl.tstr,
            'aa_alt':hl.tstr,
            'esm_score':hl.tfloat,
            'uniprot':hl.tstr,
        },
        force=True,
        missing='',
    )
esm1b_ht = esm1b_ht.rename({'uniprot': 'uniprot_isoform'})
esm1b_ht.write(FORMATTED_VSM_HT_PATHS[METHOD])    


