import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, FORMATTED_VSM_HT_PATHS

METHOD = 'POPEVE'
# PopEVE
ht = hl.import_table(
    MISSENSE_SCORE_RESOURCE_PATHS[METHOD],
    delimiter=',',
    source_file_field='filename',
    types = {
        'popEVE':hl.tfloat,
        'popped EVE':hl.tfloat,
        'popped ESM-1v':hl.tfloat,
        'EVE':hl.tfloat,
        'ESM-1v':hl.tfloat,
        
    },
    force=True,
    missing='',
)

ht2 = ht.select(
    aa_pos = hl.int(ht.mutant[1:-1]),
    aa_ref = ht.mutant[0],
    aa_alt = ht.mutant[-1],
    popEVE = ht['popEVE'], 
    popped_EVE = ht['popped EVE'], 
    popped_ESM_1v = ht['popped ESM-1v'], 
    EVE = ht['EVE'], 
    ESM_1v = ht['ESM-1v'], 
    NP = ht.filename.split('/')[-1][:-4]
    )
ht2.write(f'{WRITE_VSM_TABLES_PATH}/popeve.ht')

