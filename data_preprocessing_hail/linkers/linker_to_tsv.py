import hail as hl
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from resources.paths import ENST_TO_UNIPROT_HT
from resources.functions import write_tsv_bgz_from_ht

# ht = hl.read_table('gs://missense-scoring/mutation/everything_raw.ht')
# ht_2 = ht.explode(ht.vep.transcript_consequences)

# # # Missense
# ht_3 = ht_2.filter(ht_2.vep.transcript_consequences.most_severe_consequence == "missense_variant")
# ht_3 = ht_3.select(ht_3.vep.transcript_consequences)

chroms = [str(x) for x in range(1, 23)] + ['X', 'Y']
# for chrom in chroms:
#     ht_3_chrom = ht_3.filter(ht_3.locus.contig == 'chr' + chrom)
#     write_tsv_bgz_from_ht(ht_3_chrom, f'gs://genetics-gym/linkers/everything_raw_by_chrom_tsv/everything_raw_by_chrom_{chrom}.tsv.bgz')
#     print(f'{chrom} done')

# write_tsv_bgz_from_ht(ht_3, 'gs://genetics-gym/linkers/everything_raw.tsv.bgz')

LINKER_PATHS = {
    "MISSENSE_ENST_TRANSCRIPT": ["gs://genetics-gym/linkers/linker_missense_enst_transcript.ht", "gs://genetics-gym/linkers/linker_missense_enst_transcript_by_chrom_tsv/linker_missense_enst_transcript"],
    "MISSENSE_ENST_TRANSCRIPT_AA": ["gs://genetics-gym/linkers/linker_missense_enst_transcript_aa.ht", "gs://genetics-gym/linkers/linker_missense_enst_transcript_aa_by_chrom_tsv/linker_missense_enst_transcript_aa"],
    "MISSENSE_ENST_TRANSCRIPT_AA_UNIPROT": ["gs://genetics-gym/linkers/linker_missense_enst_transcript_aa_uniprot.ht", "gs://genetics-gym/linkers/linker_missense_enst_transcript_aa_uniprot_by_chrom_tsv/linker_missense_enst_transcript_aa_uniprot"],
    "MISSENSE_REFSEQ_TRANSCRIPT_AA": ["gs://genetics-gym/linkers/linker_missense_refseq_transcript_aa.ht", "gs://genetics-gym/linkers/linker_missense_refseq_transcript_aa_by_chrom_tsv/linker_missense_refseq_transcript_aa"],
}
def ht_by_chrom(ht, path):
    print(f'Writing {path} by chrom')
    for chrom in chroms:
        ht_chrom = ht.filter(ht.locus.contig == 'chr' + chrom)
        write_tsv_bgz_from_ht(ht_chrom, f'{path}_{chrom}.tsv.bgz')
        print(f'{chrom} done')
    print(f'{path} done')


for key, path in LINKER_PATHS.items():
    ht_by_chrom(hl.read_table(path[0]), path[1])

# for key, path in LINKER_PATHS.items():
#     to_write = path[0][:-3]
#     write_tsv_bgz_from_ht(hl.read_table(path[0]), f'{to_write}.tsv.bgz')
#     print(f'{path[1]} done')