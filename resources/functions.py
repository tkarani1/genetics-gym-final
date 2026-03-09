import hail as hl
import json

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

def write_raw_and_collected_counts(ht, counts_file_path, method_name, keyed_by, filter_cols):
    if filter_cols:
        for c in filter_cols:
            ht = ht.filter(hl.is_defined(ht[c]))
    raw_count = ht.count()
    print(ht.describe())
    ht_collected = ht.collect_by_key()
    collected_count = ht_collected.count()

    # Write counts to JSON
    counts_file = f'{counts_file_path}/{method_name}_counts.json'
    # Load existing counts if file exists
    if hl.hadoop_exists(counts_file):
        with hl.hadoop_open(counts_file, 'r') as f:
            counts = json.load(f)
    else:
        counts = {}

    # Update with AM counts
    counts[method_name] = {
        f'raw isoforms {keyed_by}': raw_count,
        f'collected {keyed_by}': collected_count,
    }

    # Write updated counts to cloud
    with hl.hadoop_open(counts_file, 'w') as f:
        json.dump(counts, f, indent=2)
        
    print(f"Counts written to: {counts_file}")
