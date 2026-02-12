import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources import POPEVE_PATH
from resources.paths import WRITE_VSM_TABLES_PATH
from resources.functions import locus_alleles_to_chr_pos_ref_alt
# PopEVE
popeve_raw_ht = hl.import_table(
    POPEVE_PATH,
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

popeve_ht = popeve_raw_ht.select(
    aa_pos = hl.int(popeve_raw_ht.mutant[1:-1]),
    aa_ref = popeve_raw_ht.mutant[0],
    aa_alt = popeve_raw_ht.mutant[-1],
    popEVE = popeve_raw_ht['popEVE'], 
    popped_EVE = popeve_raw_ht['popped EVE'], 
    popped_ESM_1v = popeve_raw_ht['popped ESM-1v'], 
    EVE = popeve_raw_ht['EVE'], 
    ESM_1v = popeve_raw_ht['ESM-1v'], 
    NP = popeve_raw_ht.filename.split('/')[-1][:-4]
    )
popeve_ht.write(f'{WRITE_VSM_TABLES_PATH}/popeve_2025.ht')
popeve_ht = locus_alleles_to_chr_pos_ref_alt(popeve_ht)
popeve_ht.export(f'{WRITE_VSM_TABLES_PATH}/popeve_2025.tsv.bgz')
