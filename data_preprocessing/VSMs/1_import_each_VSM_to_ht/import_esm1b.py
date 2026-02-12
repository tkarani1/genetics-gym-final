import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources import ESM1B_PATH
from resources.paths import WRITE_VSM_TABLES_PATH
from resources.functions import locus_alleles_to_chr_pos_ref_alt
# # ESM1b
esm1b_ht = hl.import_table(
        ESM1B_PATH,
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
esm1b_ht.write(f'{WRITE_VSM_TABLES_PATH}/esm1b_scores.ht')
esm1b_ht = locus_alleles_to_chr_pos_ref_alt(esm1b_ht)
esm1b_ht.export(f'{WRITE_VSM_TABLES_PATH}/esm1b_scores.tsv.bgz')


