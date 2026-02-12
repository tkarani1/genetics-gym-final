import hail as hl
hl.init(worker_memory="highmem", driver_memory='highmem')   
from resources import PAI3D_PATH
from resources.paths import WRITE_VSM_TABLES_PATH

## PAI3D
ht = hl.import_csv(PAI3D_PATH, min_partitions=100)
ht = ht.key_by(locus=hl.locus(ht.chr, hl.int(ht.pos), 'GRCh38'),
               alleles=[ht.non_flipped_ref, ht.non_flipped_alt],
               enst = ht.gene_name.split('\.')[0])
ht = ht.annotate(score_PAI3D=hl.float(ht.score_PAI3D))
ht = ht.select(ht.score_PAI3D)
ht.write(f'{WRITE_VSM_TABLES_PATH}/PAI3D.ht')