import hail as hl
hl.init(gcs_requester_pays_configuration='nnf-karczewski')
gnomad_v4 = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht')
gnomad_v2 = hl.read_table('gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht')
ukb = hl.read_table('gs://ukbb-exome-public/500k/results/vep.ht')
everything = hl.read_table('gs://missense-scoring/mutation/everything.ht')

gnomad_independent = everything.annotate(
    new_in_v4 = hl.is_missing(gnomad_v2[everything.key]) & hl.is_missing(ukb[everything.key]) & hl.is_defined(gnomad_v4[everything.key]),
    not_in_datasets = hl.is_missing(gnomad_v2[everything.key]) & hl.is_missing(ukb[everything.key]) & hl.is_missing(gnomad_v4[everything.key]),
    AF = gnomad_v4[everything.key].freq.AF[0],
    AC = gnomad_v4[everything.key].freq.AC[0],
    AN = gnomad_v4[everything.key].freq.AN[0]
)
gnomad_independent = gnomad_independent.filter(gnomad_independent.new_in_v4 | gnomad_independent.not_in_datasets)

ht = gnomad_independent.annotate(is_pos = ~gnomad_independent.new_in_v4)

ht.write('gs://nnf-parsa/gnomad-independent-set.ht',True)

ht = hl.read_table('gs://nnf-parsa/gnomad-independent-set.ht')
p = 1000000./ht.count()
ht = ht.sample(p)
ht.export('gs://missense-scoring/gnomad-independent-set-1M-subset.tsv.bgz')
