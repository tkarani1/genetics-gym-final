def locus_alleles_to_chr_pos_ref_alt(ht):
    ht = ht.annotate(
        chr = ht.locus.contig,
        pos = ht.locus.position,
        ref = ht.alleles[0],
        alt = ht.alleles[1]
    )
    ht = ht.key_by()
    ht.drop('locus', 'alleles')
    return ht

