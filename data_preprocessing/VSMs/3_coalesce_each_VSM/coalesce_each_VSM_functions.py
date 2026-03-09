import hail as hl
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from resources.constants import Direction
from resources.paths import TEMP_PATH
# hl.init(worker_memory="highmem", driver_memory='highmem') 


def coalesce_VSM(ht_path, key_by, score_dict):
    ht = hl.read_table(ht_path)
    score_cols = list(score_dict.keys())

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
        to_drop = [f'{s}_mane_scores', f'{s}_canon_scores', f'{s}_any_scores']
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

def coalesce_VSM_old(ht_path, key_by, score_dict):
    ht = hl.read_table(ht_path)
    ht = ht.key_by(*key_by)
    ht = ht.collect_by_key()
    score_cols = list(score_dict.keys())
    final_scores = []
    for s in score_cols:
        ht = ht.annotate(
            **{f'{s}_mane_scores':  hl.filter(lambda x: x.mane_select == True, ht.values).map(lambda x: x[s])},
            **{f'{s}_canon_scores':  hl.filter(lambda x: x.canonical == True, ht.values).map(lambda x: x[s])},
            **{f'{s}_any_scores':  hl.filter(lambda x: hl.is_defined(x[s]), ht.values).map(lambda x: x[s])},
        )
        if score_dict[s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS:
            ht = ht.annotate(
                **{f'{s}_mane_scores_neg': ht[f'{s}_mane_scores'].map(lambda x: 1 - x)},
                **{f'{s}_canon_scores_neg': ht[f'{s}_canon_scores'].map(lambda x: 1 - x)},
                **{f'{s}_any_scores_neg': ht[f'{s}_any_scores'].map(lambda x: 1 - x)},
            )
            ht = ht.annotate(
                **{f'{s}_mane_max': hl.max(ht[f'{s}_mane_scores_neg'])},
                **{f'{s}_canon_max': hl.max(ht[f'{s}_canon_scores_neg'])},
                **{f'{s}_any_max': hl.max(ht[f'{s}_any_scores_neg'])},
            )
            ht = ht.annotate(
                **{f'{s}_neg': hl.coalesce(ht[f'{s}_mane_max'], ht[f'{s}_canon_max'], ht[f'{s}_any_max'])}
            )
            to_drop = [f'{s}_mane_scores_neg', f'{s}_canon_scores_neg', f'{s}_any_scores_neg', f'{s}_mane_max', f'{s}_canon_max', f'{s}_any_max']
            ht = ht.drop(*to_drop)
            final_scores.append(f'{s}_neg')
        else: 
            # append
            ht = ht.annotate(
                **{f'{s}_mane_max': hl.max(ht[f'{s}_mane_scores'])},
                **{f'{s}_canon_max': hl.max(ht[f'{s}_canon_scores'])},
                **{f'{s}_any_max': hl.max(ht[f'{s}_any_scores'])},
            )
            ht = ht.annotate(
                **{f'{s}': hl.coalesce(ht[f'{s}_mane_max'], ht[f'{s}_canon_max'], ht[f'{s}_any_max'])}
            )
            to_drop = [f'{s}_mane_max', f'{s}_canon_max', f'{s}_any_max']
            ht = ht.drop(*to_drop)
            final_scores.append(f'{s}')
    return ht, final_scores



def coalesce_VSM_wrong(ht_path, key_by, score_dict):
    ht = hl.read_table(ht_path)
    ht = ht.key_by(*key_by)
    ht = ht.collect_by_key()
    
    score_cols = list(score_dict.keys())
    final_scores = []
    
    # Build all annotations in a single annotate call to avoid stream reuse
    annotations = {}
    for original_s in score_cols:
        output_s = f'{original_s}_neg' if score_dict[original_s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS else original_s
        should_negate = score_dict[original_s]['sense'] == Direction.HIGHER_IS_LESS_DELETERIOUS
        
        # Build single expression for this score
        if should_negate:
            score_expr = hl.coalesce(
                hl.max(ht.values.filter(lambda x: x.mane_select == True).map(lambda x: -1 * x[original_s])),
                hl.max(ht.values.filter(lambda x: x.canonical == True).map(lambda x: -1 * x[original_s])),
                hl.max(ht.values.filter(lambda x: hl.is_defined(x[original_s])).map(lambda x: -1 * x[original_s]))
            )
        else:
            score_expr = hl.coalesce(
                hl.max(ht.values.filter(lambda x: x.mane_select == True).map(lambda x: x[original_s])),
                hl.max(ht.values.filter(lambda x: x.canonical == True).map(lambda x: x[original_s])),
                hl.max(ht.values.filter(lambda x: hl.is_defined(x[original_s])).map(lambda x: x[original_s]))
            )
        
        annotations[output_s] = score_expr
        final_scores.append(output_s)
    
    # Apply all annotations at once
    ht = ht.annotate(**annotations)
    
    return ht, final_scores