#!/usr/bin/env python
# coding: utf-8

import hail as hl
from gnomad.utils.vep import process_consequences, CSQ_CODING
from gnomad.utils.constraint import annotate_mutation_type, trimer_from_heptamer, collapse_strand
import numpy as np
import pandas as pd
import os

def create_broadcast_dict(key, value = None):
    """
    Create broadcast join (local dictionary from key -> value)
    from a Hail Table.

    :param Expression key: Key Expression
    :param Expression value: Value Expression
    :return: Hail DictExpression (without an index)
    :rtype: DictExpression
    """
    if isinstance(key, hl.Table):
        key = key.key
    ht = key._indices.source
    if value is None:
        value = ht.row_value
    return hl.dict(ht.aggregate(hl.agg.collect((key, value)), _localize=False))

def annotate_with_mu(
    ht: hl.Table,
    mutation_ht: hl.Table,
    mu_annotation: str = "mu_snp",
) -> hl.Table:
    """
    Annotate SNP mutation rate for the input Table.

    .. note::

        Function expects that`ht` includes`mutation_ht`'s key fields. Note that these
        annotations don't need to be the keys of `ht`.

    :param ht: Input Table to annotate.
    :param mutation_ht: Mutation rate Table.
    :param mu_annotation: The name of mutation rate annotation in `mutation_ht`.
        Default is 'mu_snp'.
    :return: Table with mutational rate annotation added.
    """
    mu = create_broadcast_dict(mutation_ht, mutation_ht[mu_annotation]).get(hl.struct(**{k: ht[k] for k in mutation_ht.key}))
    return ht.annotate(
        **{mu_annotation: hl.case().when(hl.is_defined(mu), mu).or_error("Missing mu")}
    )

def step1():
    cov = hl.read_table("gs://gcp-public-data--gnomad/release/4.0/coverage/exomes/gnomad.exomes.v4.0.coverage.ht")
    meth = hl.read_table("gs://gcp-public-data--gnomad/resources/grch38/methylation_sites/methylation.ht")
    fi = hl.read_table("gs://gcp-public-data--gnomad/resources/context/grch38_context_vep_annotated.v105.ht/")
    mut = hl.read_table("gs://missense-scoring/mutation/gnomad.v4.0.mutation_rate.ht")

    fi = fi.filter(fi.locus.in_autosome() | (fi.locus.contig == 'chrX'))

    fi = fi.annotate(coverage=cov[fi.locus].median_approx,
    meth_levels=meth[fi.locus].methylation_level)

    fi = trimer_from_heptamer(fi)
    fi = fi.annotate(ref=fi.alleles[0], alt=fi.alleles[1])

    fi2 = fi.filter((hl.len(fi.ref) == 1) & (hl.len(fi.alt) == 1) & fi.context.matches('[ATCG]{3}') & (hl.len(fi.context) == 3))

    fi2 = annotate_mutation_type(fi2)
    fi2 = fi2.annotate(methylation_level = hl.coalesce(fi2.meth_levels, 0))

    fi2 = collapse_strand(fi2)

    methylation_expr = fi2.methylation_level
    methylation_cutoffs = hl.if_else(fi2.locus.contig != "chrX", (5, 0), (3, 0))

    fi2 = fi2.annotate(
            methylation_level=(
                hl.case()
                .when(fi2.cpg & (methylation_expr > methylation_cutoffs[0]), 2)
                .when(fi2.cpg & (methylation_expr > methylation_cutoffs[1]), 1)
                .default(0)
            ),
        )

    fi2_1 = annotate_with_mu(fi2, mut)

    fi3 = process_consequences(fi2_1)
    fi3 = fi3.filter(hl.literal(CSQ_CODING).contains(fi3.vep.worst_consequence_term))

    fi3.write('gs://missense-scoring/mutation/everything_raw.ht', overwrite = True)

def step2():
    fi3 = hl.read_table('gs://missense-scoring/mutation/everything_raw.ht')

    fi3 = fi3.annotate(tc_filtered=fi3.vep.transcript_consequences.filter(lambda x: hl.is_defined(x.mane_select) & x.transcript_id.startswith('ENST') & hl.literal(set(CSQ_CODING)).contains(x.most_severe_consequence)))

    fi3_1 = fi3.explode(fi3.tc_filtered)

    # print(fi3_1.aggregate(hl.agg.counter(hl.len((fi3_1.tc_filtered.uniprot_isoform)))))

    fi3_1 = fi3_1.annotate(uniprot = hl.or_missing(hl.len(fi3_1.tc_filtered.uniprot_isoform) > 0, fi3_1.tc_filtered.uniprot_isoform[0]))

    output = fi3_1.drop("vep")

    pattern = r".*:p.(\D+)(\d+)(\D+)"
    array_match = output.tc_filtered.hgvsp.first_match_in(pattern)
    output = output.annotate(aa_ref = array_match[0], aa_pos = hl.int(array_match[1]), aa_alt = array_match[2])

    output = output.naive_coalesce(5000)

    output.write('gs://missense-scoring/mutation/everything.ht', overwrite = True)


def main():
    step1()
    step2()

if __name__ == '__main__':
    main()