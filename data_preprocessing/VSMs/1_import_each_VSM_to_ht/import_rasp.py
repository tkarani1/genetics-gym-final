import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources import RASP_PATH
from resources.paths import WRITE_VSM_TABLES_PATH
from resources.functions import locus_alleles_to_chr_pos_ref_alt

ht = hl.read_table(RASP_PATH)
ht = ht.naive_coalesce(50)
ht = ht.annotate(
    uniprot_id = ht.uniprot,
    uniprot_aa_pos = hl.int(ht.variant[1:-1]),
    uniprot_aa_ref = ht.variant[0],
    uniprot_aa_alt = ht.variant[-1],
    rasp_score = ht.score_ml,
)
ht = ht.key_by(
    'uniprot_id',
    'uniprot_aa_pos',
    'uniprot_aa_ref',
    'uniprot_aa_alt',
)
ht = ht.select(ht.rasp_score)
ht.write(f'{WRITE_VSM_TABLES_PATH}/rasp_scores.ht')
ht = locus_alleles_to_chr_pos_ref_alt(ht)
ht.export(f'{WRITE_VSM_TABLES_PATH}/rasp_scores.tsv.bgz')