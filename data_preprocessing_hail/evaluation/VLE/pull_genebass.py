import hail as hl
from resources.functions import liftover_37_38
from resources.paths import LINKER_PATHS, EVALUATION_TABLE_NAMES, EVALUATION_TYPE_NAMES, EVALUATION_TABLES_BASE, EVALUATION_RESOURCE_PATHS

linker = LINKER_PATHS['MISSENSE_ONLY_SNP']
NAME = 'GENEBASS'
write_path = f"{EVALUATION_TABLES_BASE}/{EVALUATION_TABLE_NAMES[NAME]}_{EVALUATION_TYPE_NAMES['missense']}"

def create_matched_sets(ht):
    ht = ht.annotate(matching_freq_bin=10 ** (hl.floor(hl.log10(ht.matching_freq) * 3) / 3))
    res = ht.aggregate(hl.agg.group_by(ht.matching_freq_bin, hl.struct(pos=hl.agg.count_where(ht.original_is_pos),
                                                                       total=hl.agg.count_where(~ht.original_is_pos))),
                       _localize=False)
    ht = ht.annotate(is_pos=hl.case().when(ht.original_is_pos, True)
                     .when(hl.rand_bool(res[ht.matching_freq_bin].pos / res[ht.matching_freq_bin].total), False)
                     .or_missing()
                     )
    return ht


linker_ht = hl.read_table(linker) 
ht = hl.read_table('gs://missense-scoring/mutation/genebass.ht')
ht = ht.filter(ht.keep_var_expected_ac & ht.keep_var_annt)
ht = ht[linker_ht.key]    # no filtering for missense bc the linker does filtering (linker is only missense)
ht = linker_ht.annotate(
    original_is_pos=ht.all_sig_pheno_cnt > 0,  # how many traits are associated with this variant / >0 means associated with at least 1 trait
    matching_freq=ht.call_stats.AF,
)
ht = ht.filter(hl.is_defined(ht.original_is_pos))
ht = create_matched_sets(ht)
ht = ht.naive_coalesce(10)