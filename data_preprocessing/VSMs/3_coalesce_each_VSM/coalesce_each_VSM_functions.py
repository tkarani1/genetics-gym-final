import hail as hl
import json
hl.init(worker_memory="highmem", driver_memory='highmem') 

def coalesce_VSM(ht_path, key_by, higher_is_less_deleterious):
    ht = hl.read_table(ht_path)
    ht = ht.key_by(*key_by)
    ht = ht.collect_by_key()
    score_cols = higher_is_less_deleterious.keys()
    final_scores = []
    for s in score_cols:
        ht = ht.annotate(
            **{f'{s}_mane_scores':  hl.filter(lambda x: x.mane_select == True, ht.values).map(lambda x: x[s])},
            **{f'{s}_canon_scores':  hl.filter(lambda x: x.canonical == True, ht.values).map(lambda x: x[s])},
            **{f'{s}_any_scores':  hl.filter(lambda x: hl.is_defined(x[s]), ht.values).map(lambda x: x[s])},
        )
        if higher_is_less_deleterious[s]:
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