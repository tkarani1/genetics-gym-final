import hail as hl
import pandas as pd
import os
from resources.functions import check_duplicates, write_tsv_bgz_from_ht
from resources.paths import EVALUATION_TABLES_BASE

rg37 = hl.get_reference('GRCh37')  
rg38 = hl.get_reference('GRCh38')  
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38) 
def liftover_37_38(ht):
    ht = ht.annotate(
        new_locus=hl.liftover(ht.locus, 'GRCh38', include_strand=True),
        old_locus=ht.locus
    )
    ht = ht.filter(hl.is_defined(ht.new_locus) & ~ht.new_locus.is_negative_strand)  
    ht = ht.key_by()
    ht = ht.annotate(locus=ht.new_locus.result)
    ht = ht.drop('new_locus', 'old_locus')
    return ht


# extract missense variants
# extract locus, [alleles] from variant labeling
# label this dataset as pos_dd
ht_dd = hl.import_table('gs://genetics-gym/raw_data/kaplanis_variants_annotated_2024-05-15.txt')
ht_dd = ht_dd.filter(ht_dd.Consequence == 'missense_variant') 
ht_dd = ht_dd.annotate(**hl.parse_variant(ht_dd.Variant, reference_genome='GRCh37'))
ht_dd = liftover_37_38(ht_dd)
ht_dd.describe()
ht_dd = ht_dd.key_by('locus', 'alleles')
ht_dd = ht_dd.annotate(is_pos = True)
ht_dd = ht_dd.select('is_pos')
# check_duplicates(ht_dd, output_path='gs://trisha-tmp/genetics_gym_final_count_checks/eval_counts/DD_RAW_DUPLICATES_snp.tsv.bgz')
# check_duplicates(ht_dd, output_path='gs://trisha-tmp/genetics_gym_final_count_checks/eval_counts/DD_RAW_DUPLICATES_snp_gene.tsv.bgz')

# filtering to SNPs. hl.is_snp

ht_asd = hl.import_table('gs://genetics-gym/raw_data/Nat_Gen_published_autosomal_and_updated_XY_de_novo_calls_2024-05-13_no_cohort.txt')
ht_asd = ht_asd.filter(ht_asd.Consequence == 'missense_variant')
ht_asd = ht_asd.key_by(**hl.parse_variant('chr'+ht_asd.Variant, reference_genome='GRCh38'))
ht_asd = ht_asd.annotate(is_pos = ht_asd.Role == 'Proband')
ht_asd = ht_asd.select('is_pos')

# chd_cases = pd.read_csv('gs://genetics-gym/raw_data/NIHMS906719-chd-dnm-cases.csv')
# chd_cases = chd_cases.rename(columns={"AA change": "AA_change", "pLI score": "pLI Score"})
# chd_cases = chd_cases.drop(['Cardiac Category','EM','NDD'],axis=1)
# chd_cases.to_csv('gs://genetics-gym/temp/NIHMS906719-chd-dnm-cases-processed.csv', index=False)

ht_chd = hl.import_table(
    [
        'gs://genetics-gym/temp/NIHMS906719-chd-dnm-cases-processed.csv',
        'gs://genetics-gym/raw_data/NIHMS906719-chd-dnm-controls.csv'
    ],
    source_file_field='case-control',
    delimiter = ','
)
ht_chd = ht_chd.annotate(locus = hl.locus(ht_chd.CHROM, hl.int(ht_chd.POS), reference_genome='GRCh37'))
ht_chd = ht_chd.annotate(alleles = [ht_chd.REF,ht_chd.ALT])
ht_chd = liftover_37_38(ht_chd)
ht_chd = ht_chd.key_by('locus','alleles')
ht_chd = ht_chd.annotate(is_pos = hl.if_else(ht_chd['case-control'].contains('case'), True, False))
ht_chd = ht_chd.select('is_pos')

# ht_asd.write('gs://trisha-tmp/intermediates/asd-no-combined-controls.ht', overwrite=True)
# ht_dd.write('gs://trisha-tmp/intermediates/dd-no-combined-controls.ht', overwrite=True)
# ht_chd.write('gs:/trisha-tmp/intermediates/chd-no-combined-controls.ht', overwrite=True)

ht_controls_from_asd = ht_asd.filter(~ht_asd.is_pos)
ht_controls_from_asd = ht_controls_from_asd.select('is_pos')
ht_dd = ht_dd.union(ht_controls_from_asd)

ht_controls_from_chd = ht_chd.filter(~ht_chd.is_pos)
ht_controls_from_chd = ht_controls_from_chd.select('is_pos')
ht_dd = ht_dd.union(ht_controls_from_chd)

ht_asd = ht_asd.union(ht_controls_from_chd)
ht_chd = ht_chd.union(ht_controls_from_asd)

ht_asd = ht_asd.checkpoint(EVALUATION_TABLES_BASE + 'asd-with-combined-controls.ht', overwrite=True)
ht_dd = ht_dd.checkpoint(EVALUATION_TABLES_BASE + 'dd-with-combined-controls.ht', overwrite=True)
ht_chd = ht_chd.checkpoint(EVALUATION_TABLES_BASE + 'chd-with-combined-controls.ht', overwrite=True)

write_tsv_bgz_from_ht(ht_asd, EVALUATION_TABLES_BASE + 'asd-with-combined-controls.tsv.bgz')
write_tsv_bgz_from_ht(ht_dd, EVALUATION_TABLES_BASE + 'dd-with-combined-controls.tsv.bgz')
write_tsv_bgz_from_ht(ht_chd, EVALUATION_TABLES_BASE + 'chd-with-combined-controls.tsv.bgz')

# drop duplicate locus-allele keys
ht_asd = ht_asd.distinct()
ht_dd = ht_dd.distinct()
ht_chd = ht_chd.distinct()

ht_asd = ht_asd.checkpoint(EVALUATION_TABLES_BASE + 'asd-with-combined-controls-distinct.ht', overwrite=True)
ht_dd = ht_dd.checkpoint(EVALUATION_TABLES_BASE + 'dd-with-combined-controls-distinct.ht', overwrite=True)
ht_chd = ht_chd.checkpoint(EVALUATION_TABLES_BASE + 'chd-with-combined-controls-distinct.ht', overwrite=True)

write_tsv_bgz_from_ht(ht_asd, EVALUATION_TABLES_BASE + 'asd-with-combined-controls-distinct.tsv.bgz')
write_tsv_bgz_from_ht(ht_dd, EVALUATION_TABLES_BASE + 'dd-with-combined-controls-distinct.tsv.bgz')
write_tsv_bgz_from_ht(ht_chd, EVALUATION_TABLES_BASE + 'chd-with-combined-controls-distinct.tsv.bgz')

# # case getting extracted from dataset, all controls combined
