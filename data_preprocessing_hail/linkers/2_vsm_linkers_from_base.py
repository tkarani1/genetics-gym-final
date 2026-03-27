import hail as hl
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from resources.paths import ENST_TO_UNIPROT_HT

hl.init(worker_memory="highmem", driver_memory='highmem')   

def convert_aa_three_to_one(ht):
    # Mapping dictionary from three-letter to one-letter amino acid codes
    aa_dict = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
        'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
        'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
        'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
        'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        'Sec': 'U', 'Pyl': 'O', 'Asx': 'B', 'Glx': 'Z',
        'Xaa': 'X', 'Ter': '*'
    }

    # Create a Hail dictionary expression
    aa_mapping = hl.dict(aa_dict)

    # Annotate the table by mapping the three-letter codes to one-letter codes
    ht = ht.annotate(
        aa_ref = aa_mapping.get(ht.aa_ref, ht.aa_ref),
        aa_alt = aa_mapping.get(ht.aa_alt, ht.aa_alt)
    )

    return ht

ht = hl.read_table('gs://missense-scoring/mutation/everything_raw.ht')
ht_2 = ht.explode(ht.vep.transcript_consequences)

# # Missense
ht_3 = ht_2.filter(ht_2.vep.transcript_consequences.most_severe_consequence == "missense_variant")
ht_4 = ht_3.filter(ht_3.vep.transcript_consequences.transcript_id.startswith('ENST') )

## Missense only SNP
# ht_a = ht_4.select()
# ht_a = ht_a.distinct()
# ht_a.write("gs://genetics-gym/linkers/linker_missense_only_snp.ht", overwrite=True)
# ht_a = hl.read_table("gs://genetics-gym/linkers/linker_missense_only_snp.ht")

## Missense SNP + ENSG + ENST
# ht_b = ht_4.select(
#     gene_symbol = ht_4.vep.transcript_consequences.gene_symbol,
#     enst = ht_4.vep.transcript_consequences.transcript_id,
#     ensp = ht_4.vep.transcript_consequences.protein_id,
#     ensg = ht_4.vep.transcript_consequences.gene_id, 
#     mane_select = hl.is_defined(ht_4.vep.transcript_consequences.mane_select),
#     canonical = hl.is_defined(ht_4.vep.transcript_consequences.canonical)
#     )
# ht_b = ht_b.key_by('locus', 'alleles', 'enst')
# ht_b = ht_b.distinct()
# ht_b.write("gs://genetics-gym/linkers/linker_missense_enst_transcript.ht", overwrite=True)

## Missense SNP + ENSG + ENST + AA position
pattern = r".*:p.(\D+)(\d+)(\D+)"
array_match = ht_4.vep.transcript_consequences.hgvsp.first_match_in(pattern)
ht_5 = ht_4.annotate(aa_ref = array_match[0], aa_pos = hl.int(array_match[1]), aa_alt = array_match[2])

# ht_c = ht_5.select(
#     aa_pos = ht_5.aa_pos,
#     aa_ref = ht_5.aa_ref,
#     aa_alt = ht_5.aa_alt,
#     gene_symbol = ht_5.vep.transcript_consequences.gene_symbol,
#     enst = ht_5.vep.transcript_consequences.transcript_id,
#     ensp = ht_5.vep.transcript_consequences.protein_id,
#     ensg = ht_5.vep.transcript_consequences.gene_id, 
#     mane_select = hl.is_defined(ht_5.vep.transcript_consequences.mane_select),
#     canonical = hl.is_defined(ht_5.vep.transcript_consequences.canonical)
# )
# ht_c = convert_aa_three_to_one(ht_c)    
# ht_c = ht_c.key_by('locus', 'alleles', 'enst')
# ht_c = ht_c.distinct()
# ht_c.write("gs://genetics-gym/linkers/linker_missense_enst_transcript_aa.ht", overwrite=True)

## Missense SNP + ENSG + ENST + AA position + Uniprot isoform
# ht_d = hl.read_table("gs://genetics-gym/linkers/linker_missense_enst_transcript_aa.ht")
# ht_d = ht_d.key_by('enst')
# enst_to_uniprot_ht = hl.read_table(ENST_TO_UNIPROT_HT)
# enst_to_uniprot_ht = enst_to_uniprot_ht.rename({'ensembl_transcript_id': 'enst'})
# # enst_to_uniprot_ht = enst_to_uniprot_ht.key_by('enst') # should already be keyed by enst
# ht_d = ht_d.join(enst_to_uniprot_ht, how='left')  #TODO: CHANGE TO OUTER???
# ht_d.write("gs://genetics-gym/linkers/linker_missense_enst_transcript_aa_uniprot.ht", overwrite=True)
# linker_enst_set = ht_d.aggregate(hl.agg.collect_as_set(ht_d.enst))
# mapping_enst_set = enst_to_uniprot_ht.aggregate(hl.agg.collect_as_set(enst_to_uniprot_ht.enst))
# print(len(linker_enst_set), len(mapping_enst_set))
# print('in linker but not in mapping: ', len(linker_enst_set - mapping_enst_set))
# print('in mapping but not in linker: ', len(mapping_enst_set - linker_enst_set))
# ht_d_f = ht_d.filter(ht_d.mane_select != ht_d.transcript_mane_select)
# print(ht_d_f.count())

## Missense SNP + ENSG + ENST + AA position  + Refseq transcript
# ht_4_v2 = ht_3.filter(ht_3.vep.transcript_consequences.transcript_id.startswith('NM') )
# pattern = r".*:p.(\D+)(\d+)(\D+)"
# array_match = ht_4_v2.vep.transcript_consequences.hgvsp.first_match_in(pattern)
# ht_5_v2 = ht_4_v2.annotate(aa_ref = array_match[0], aa_pos = hl.int(array_match[1]), aa_alt = array_match[2])

# ht_e = ht_5_v2.select(
#     aa_pos = ht_5_v2.aa_pos,
#     aa_ref = ht_5_v2.aa_ref,
#     aa_alt = ht_5_v2.aa_alt,
#     gene_symbol = ht_5_v2.vep.transcript_consequences.gene_symbol,
#     NM = ht_5_v2.vep.transcript_consequences.transcript_id,
#     NP = ht_5_v2.vep.transcript_consequences.protein_id,
#     ncbi = ht_5_v2.vep.transcript_consequences.gene_id, 
#     mane_select = hl.is_defined(ht_5_v2.vep.transcript_consequences.mane_select),
#     canonical = hl.is_defined(ht_5_v2.vep.transcript_consequences.canonical),
#     mane_ENST = ht_5_v2.vep.transcript_consequences.mane_select,
# )
# ht_e = convert_aa_three_to_one(ht_e)    
# ht_e = ht_e.key_by('locus', 'alleles', 'NM')
# ht_e = ht_e.distinct()
# ht_e.write("gs://genetics-gym/linkers/linker_missense_refseq_transcript_aa.ht", overwrite=True)
# ht_ef = ht_e.filter(~hl.is_defined(ht_e.gene_symbol))
# print(ht_ef.count())

# ht_ens = hl.read_table("gs://genetics-gym/linkers/linker_missense_enst_transcript_aa.ht")
# ht_ens_f = ht_ens.filter(~hl.is_defined(ht_ens.gene_symbol))
# print(ht_ens_f.count())

