import hail as hl
from resources.functions import liftover_37_38
from resources.paths import LINKER_PATHS, EVALUATION_TABLE_NAMES, EVALUATION_TYPE_NAMES, EVALUATION_TABLES_BASE

linker = LINKER_PATHS['MISSENSE_ONLY_SNP']
NAME = 'SCHEMA'
write_path = f"{EVALUATION_TABLES_BASE}/{EVALUATION_TYPE_NAMES[NAME]}_{EVALUATION_TYPE_NAMES['missense']}"

def munge_rvas(name, input_path, liftover=True):
    ht_linker = hl.read_table(final_linker)
    ht_rvas = hl.import_vcf(input_path, array_elements_required=False,
                            reference_genome='GRCh37' if liftover else 'GRCh38'
                            ).rows()
    if name == 'asc':
        ht_rvas = ht_rvas.filter(~ht_rvas.info.groups.contains('ASC_DN'))
        ht_rvas = ht_rvas.annotate(
            info=ht_rvas.info.annotate(
                ac_case=[hl.sum(ht_rvas.info.ac_case)],
                ac_ctrl=[hl.sum(ht_rvas.info.ac_ctrl)]
            )
        )
    # case or control has to be 0
    ht_rvas = ht_rvas.filter((ht_rvas.info.ac_case[-1] == 0) | (ht_rvas.info.ac_ctrl[-1] == 0))
    if liftover:
        ht_rvas = liftover_37_38(ht_rvas).key_by('locus', 'alleles')
    rvas = ht_rvas[ht_linker.key].info
    ht = ht_linker.annotate(
        is_pos=rvas.ac_case[-1] > 0
    )
    ht = ht.filter(hl.is_defined(ht.is_pos))
    ht.naive_coalesce(10).write(f'gs://missense-scoring/mutation/{name}_evaluation_table.ht', True)



if __name__ == '__main__':
    # munge_ukb()
    # munge_rvas('schema', paths['schema'])
    munge_rvas('asc', paths['asc'])
    munge_rvas('epi25', paths['epi25'], False)
