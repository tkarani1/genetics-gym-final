import hail as hl
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from resources.constants import Direction
from resources.paths import TEMP_PATH
# hl.init(worker_memory="highmem", driver_memory='highmem') 

def coalesce_VSM(ht_path, key_by, score_dict, keep_max_per_category=True): 
    ht = hl.read_table(ht_path)
    score_cols = list(score_dict.keys())
    ht = ht.select(*key_by, *score_cols)

    for s in score_cols:
        if score_dict[s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS:
            ht = ht.annotate(**{f'{s}_neg': -1 * ht[s]})
    ht = ht.key_by(*key_by)
    print(1)
    print(ht.describe())
    ht = ht.collect_by_key()
    print(2)
    print(ht.describe())
    ht = ht.checkpoint(TEMP_PATH + 'coalesce_VSM_temp.ht', overwrite=True)
    print(3)
    print(ht.describe())
    for s in score_cols:
        if score_dict[s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS:
            s = f'{s}_neg'
        
        ht = ht.annotate(
            **{f'{s}_mane_scores':  hl.filter(lambda x: x.mane_select == True, ht.values).map(lambda x: x[s])},
            **{f'{s}_canon_scores':  hl.filter(lambda x: x.canonical == True, ht.values).map(lambda x: x[s])},
            **{f'{s}_any_scores':  hl.filter(lambda x: hl.is_defined(x[s]), ht.values).map(lambda x: x[s])},
        )

        ht = ht.annotate(
            **{f'{s}_mane_max': hl.max(ht[f'{s}_mane_scores'])},
            **{f'{s}_canon_max': hl.max(ht[f'{s}_canon_scores'])},
            **{f'{s}_any_max': hl.max(ht[f'{s}_any_scores'])},
        )
        ht = ht.annotate(
            **{f'{s}': hl.coalesce(ht[f'{s}_mane_max'], ht[f'{s}_canon_max'], ht[f'{s}_any_max'])}
        )
        to_drop = [f'{s}_mane_scores', f'{s}_canon_scores', f'{s}_any_scores', 'values']
        if not keep_max_per_category:
            to_drop.extend([f'{s}_mane_max', f'{s}_canon_max', f'{s}_any_max'])
        ht = ht.drop(*to_drop)
    ht = ht.key_by(*key_by)
    return ht

def coalesce_VSM_no_hierarchy(ht_path, key_by, score_dict):
    ht = hl.read_table(ht_path)
    score_cols = list(score_dict.keys())
    for s in score_cols:
        if score_dict[s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS:
            ht = ht.annotate(**{f'{s}_neg': -1 * ht[s]})

    ht = ht.key_by(*key_by)
    ht = ht.collect_by_key()
    ht = ht.checkpoint(TEMP_PATH + 'coalesce_VSM_no_hierarchy_temp.ht', overwrite=True)

    for s in score_cols:
        if score_dict[s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS:
            s = f'{s}_neg'
        ht = ht.annotate(**{f'{s}': hl.max(ht.values.map(lambda x: x[s]))})
    return ht

    # should_negate = score_dict[original_s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS
