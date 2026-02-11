import hail as hl
hl.init(worker_memory="highmem")    # driver memory to highmem

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
print(0)
ht_2 = ht.explode(ht.vep.transcript_consequences)

print(1)

# Missense
ht_3 = ht_2.filter(ht_2.vep.transcript_consequences.most_severe_consequence == "missense_variant")
ht_4 = ht_3.filter(ht_3.vep.transcript_consequences.transcript_id.startswith('ENST') )

## Missense only SNP
ht_a = ht_4.select()
ht_a = ht_a.distinct()
ht_a.write("gs://genetics-gym/linkers/linker_missense_only_snp.ht", overwrite=True)


## Missense SNP + ENSG + ENST
ht_b = ht_4.select(
    gene_symbol = ht_4.vep.transcript_consequences.gene_symbol,
    enst = ht_4.vep.transcript_consequences.transcript_id,
    ensp = ht_4.vep.transcript_consequences.protein_id,
    ensg = ht_4.vep.transcript_consequences.gene_id, 
    mane_select = hl.is_defined(ht_4.vep.transcript_consequences.mane_select),
    canonical = hl.is_defined(ht_4.vep.transcript_consequences.canonical)
    )
ht_b = ht_b.key_by('locus', 'allele', 'ensg', 'enst')
ht_b = ht_b.distinct()
ht_b.write("gs://genetics-gym/linkers/linker_missense_transcript.ht", overwrite=True)


