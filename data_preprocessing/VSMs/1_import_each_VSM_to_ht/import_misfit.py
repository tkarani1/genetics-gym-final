import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, WRITE_VSM_TABLES_PATH, FORMATTED_VSM_HT_PATHS

METHOD = 'MISFIT'

# # MisFit
misfit_ht = hl.import_table(MISSENSE_SCORE_RESOURCE_PATHS[METHOD], 
                types = {'Chrom': hl.tstr, 'Pos': hl.tint, 'Ref': hl.tstr, 'Alt': hl.tstr, 
                'MisFit_D': hl.tfloat, 'MisFit_S': hl.tfloat},
                missing='', force=True)
misfit_ht = misfit_ht.annotate(
    locus = hl.locus('chr'+misfit_ht.Chrom, misfit_ht.Pos, reference_genome='GRCh38'),
)
misfit_ht.describe()
# misfit_ht = misfit_ht.rename({'Chrom': 'chrom', 'Pos': 'pos', 'Ref': 'ref', 'Alt': 'alt', 
#                     'Symbol': 'gene_symbol', 'GeneID': 'ensg', 'TranscriptID': 'enst',
#                     'UniprotID': 'uniprot_id', 'Ensembl_protein_position': 'ens_pos',
#                     'Uniprot_position': 'uniprot_pos', 'AA_Ref': 'aa_ref', 'AA_alt': 'aa_alt'})

# misfit_ht.write(FORMATTED_VSM_HT_PATHS[METHOD], overwrite=True)