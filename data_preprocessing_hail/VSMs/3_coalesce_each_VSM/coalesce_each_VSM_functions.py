import hail as hl
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from resources.constants import Direction
# hl.init(backend='spark', worker_memory="highmem", driver_memory='highmem') 

def coalesce_VSM(ht_path, key_by, score_dict, keep_max_per_category=True): 
    ht = hl.read_table(ht_path)
    score_cols = list(score_dict.keys())

    for s in score_cols:
        if score_dict[s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS:
            ht = ht.annotate(**{f'{s}_neg': -1 * ht[s]})

    grouped = ht.group_by(*[ht[k] for k in key_by])
    agg_exprs = {}

    max_field_names = []
    for s in score_cols:
        value_field = f'{s}_neg' if score_dict[s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS else s

        mane_max = hl.agg.filter(ht.mane_select & hl.is_defined(ht[value_field]), hl.agg.max(ht[value_field]))
        canon_max = hl.agg.filter(ht.canonical & hl.is_defined(ht[value_field]), hl.agg.max(ht[value_field]))
        any_max = hl.agg.filter(hl.is_defined(ht[value_field]), hl.agg.max(ht[value_field]))

        agg_exprs[s] = hl.coalesce(mane_max, canon_max, any_max)
        if keep_max_per_category:
            mane_name = f'{s}_mane_max'
            canon_name = f'{s}_canon_max'
            any_name = f'{s}_any_max'
            agg_exprs[mane_name] = mane_max
            agg_exprs[canon_name] = canon_max
            agg_exprs[any_name] = any_max
            max_field_names.extend([mane_name, canon_name, any_name])

    ht = grouped.aggregate(**agg_exprs)
    ht = ht.key_by(*key_by)
    if keep_max_per_category:
        return ht.select(*score_cols, *max_field_names)
    return ht.select(*score_cols)

def coalesce_VSM_no_hierarchy(ht_path, linker_ht_path, key_by, score_dict):
    ht = hl.read_table(ht_path)
    linker_ht = hl.read_table(linker_ht_path)
    ht = ht.join(linker_ht, how='inner')    # keep only missense SNPs
    score_cols = list(score_dict.keys())
    for s in score_cols:
        if score_dict[s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS:
            ht = ht.annotate(**{f'{s}_neg': -1 * ht[s]})

    grouped = ht.group_by(*[ht[k] for k in key_by])
    agg_exprs = {}

    for s in score_cols:
        value_field = f'{s}_neg' if score_dict[s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS else s
        agg_exprs[s] = hl.agg.filter(hl.is_defined(ht[value_field]), hl.agg.max(ht[value_field]))

    ht = grouped.aggregate(**agg_exprs)
    return ht.key_by(*key_by)

    # should_negate = score_dict[original_s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS
