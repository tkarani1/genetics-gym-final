import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem')   
from resources.paths import MISSENSE_SCORE_RESOURCE_PATHS, WRITE_VSM_TABLES_PATH, FORMATTED_VSM_HT_PATHS

METHOD = 'PAI3D'

## PAI3D
ht = hl.import_csv(MISSENSE_SCORE_RESOURCE_PATHS[METHOD], min_partitions=100)
ht = ht.key_by(locus=hl.locus(ht.chr, hl.int(ht.pos), 'GRCh38'),
               alleles=[ht.non_flipped_ref, ht.non_flipped_alt],
               enst = ht.gene_name.split('\.')[0])
ht = ht.annotate(score_PAI3D=hl.float(ht.score_PAI3D))
ht = ht.select(ht.score_PAI3D)
ht.write(FORMATTED_VSM_HT_PATHS[METHOD])    