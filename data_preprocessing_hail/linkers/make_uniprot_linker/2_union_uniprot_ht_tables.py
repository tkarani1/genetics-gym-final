import hail as hl
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from resources.files import ENST_TO_UNIPROT_FOLDER
hl.init(worker_memory="highmem", driver_memory='highmem')  

chrs = [str(x) for x in range(1, 23)] + ['X']

## union first two
c_0 = chrs[0]
ht_0 = hl.read_table(f'gs://genetics-gym/linkers/temp_uniprot/uniprot_mapping_chr_{c_0}.ht')

c_1 = chrs[1]
ht_1 = hl.read_table(f'gs://genetics-gym/linkers/temp_uniprot/uniprot_mapping_chr_{c_1}.ht')

union_ht = ht_0.union(ht_1)
union_ht = union_ht.checkpoint(f'gs://genetics-gym/linkers/temp_uniprot/uniprot_mapping_union_chr_{c_0}_and_{c_1}.ht', overwrite=True)

for c in chrs[2:]:
    ht = hl.read_table(f'gs://genetics-gym/linkers/temp_uniprot/uniprot_mapping_chr_{c}.ht')
    union_ht = union_ht.union(ht)
    union_ht = union_ht.checkpoint(f'gs://genetics-gym/linkers/temp_uniprot/uniprot_mapping_union_chr_{c}.ht', overwrite=True)

union_ht.write('gs://genetics-gym/linkers/enst_to_uniprot_mapping.ht', overwrite=True)