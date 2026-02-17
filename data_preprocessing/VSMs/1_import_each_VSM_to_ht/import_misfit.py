import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, WRITE_VSM_TABLES_PATH, FORMATTED_VSM_HT_PATHS

METHOD = 'MISFIT'
MISFIT_MAPPING_PATH = MISSENSE_SCORE_RESOURCE_PATHS['MISFIT_MAPPING']

# # MisFit
misfit_ht = hl.import_table(
       MISSENSE_SCORE_RESOURCE_PATHS[METHOD],
        source_file_field='filename',
        types = {
            'Uniprot_position':hl.tint,
            'AA_alt':hl.tstr,
            'MisFit_D':hl.tfloat,
            'MisFit_S':hl.tfloat,
        },
        force=True,
        missing='',
    )
misfit_ht = misfit_ht.annotate(
    uniprot_id = misfit_ht.filename.split('/')[-1].split('\.')[0],
)
misfit_ht = misfit_ht.rename({'Uniprot_position': 'aa_pos', 'AA_alt': 'aa_alt'})
misfit_ht = misfit_ht.drop('filename')

mf_mapping = hl.import_table(MISFIT_MAPPING_PATH)
mf_mapping = mf_mapping.key_by('UniprotID')
misfit_ht = misfit_ht.annotate(
    enst = mf_mapping[misfit_ht.uniprot_id].TranscriptID, 
    ensg = mf_mapping[misfit_ht.uniprot_id].GeneID, 
    ensp = mf_mapping[misfit_ht.uniprot_id].ProteinID, 
    gene_symbol = mf_mapping[misfit_ht.uniprot_id].Symbol
)
misfit_ht.write(FORMATTED_VSM_HT_PATHS[METHOD])    