import hail as hl
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))
from resources.files import ENST_TO_UNIPROT_FOLDER
hl.init(worker_memory="highmem", driver_memory='highmem')   

uniprot_mapping_path = ENST_TO_UNIPROT_FOLDER

c= sys.argv[1]

# format uniprot mapping file
fname = f'chr{c}_ensembl_tid_to_uniprot.tsv'
ht_uni = hl.import_table(uniprot_mapping_path + fname, delimiter='\t', 
                                types = {'transcript_mane_select':hl.tbool},
                                force=True,
                                missing='')
ht_uni = ht_uni.annotate(
    uniprot_isoform = hl.coalesce(ht_uni.uniprot_isoform, ht_uni.uniprotswissprot, ht_uni.uniprotsptrembl)

)
ht_uni = ht_uni.annotate(
    uniprot_id = ht_uni.uniprot_isoform.split('-')[0]
)

ht_uni = ht_uni.key_by('ensembl_transcript_id')
ht_uni = ht_uni.select('uniprot_isoform', 'uniprot_id', 'transcript_mane_select')
ht_uni.write(f'gs://genetics-gym/linkers/temp_uniprot/uniprot_mapping_chr_{c}.ht', overwrite=True)
