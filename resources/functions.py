import hail as hl
import json

def locus_alleles_to_chr_pos_ref_alt(ht):
    ht = ht.annotate(
        chrom = ht.locus.contig,
        pos = ht.locus.position,
        ref = ht.alleles[0],
        alt = ht.alleles[1]
    )
    ht = ht.key_by()
    ht.drop('locus', 'alleles')
    return ht

def write_tsv_bgz_from_ht_path(ht_path, tsv_bgz_path):
    ht = hl.read_table(ht_path)
    ht = locus_alleles_to_chr_pos_ref_alt(ht)
    if not tsv_bgz_path.endswith('.tsv.bgz'):
        tsv_bgz_path += '.tsv.bgz'
    ht.export(tsv_bgz_path, delimiter='\t')

def write_tsv_bgz_from_ht(ht, tsv_bgz_path):
    ht = locus_alleles_to_chr_pos_ref_alt(ht)
    if not tsv_bgz_path.endswith('.tsv.bgz'):
        tsv_bgz_path += '.tsv.bgz'
    ht.export(tsv_bgz_path, delimiter='\t')


def check_duplicates(ht, key_by=None, output_path=None):
    if key_by is not None:  # assume ht is already keyed 
        ht = ht.key_by(*key_by)
    ht_grouped = ht.collect_by_key()
    print('HT before grouping: ', ht.count())
    print('HT after grouping: ', ht_grouped.count())
    ht_grouped_filtered = ht_grouped.filter(hl.len(ht_grouped.values) > 1)
    print('Number of duplicated keys: ', ht_grouped_filtered.count())

    if output_path is not None:
        write_tsv_bgz_from_ht(ht_grouped_filtered, output_path)
    return

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



# def write_raw_and_collected_counts(ht, counts_file_path, method_name, keyed_by, Optional[filter_methods] = None):
#     ## TODO: neg score handling
#     ## TODO: filter to method is defined
#     # if filter_cols is not None:
#     #     for c in filter_cols:
#     #         ht = ht.filter(hl.is_defined(ht[c]))
#     raw_count = ht.count()
#     print(ht.describe())
#     ht_collected = ht.collect_by_key()
#     collected_count = ht_collected.count()

#     # Write counts to JSON
#     counts_file = f'{counts_file_path}/{method_name}_counts.json'
#     # Load existing counts if file exists
#     if hl.hadoop_exists(counts_file):
#         with hl.hadoop_open(counts_file, 'r') as f:
#             counts = json.load(f)
#     else:
#         counts = {}

#     # Update with AM counts
#     counts[method_name] = {
#         f'raw isoforms {keyed_by}': raw_count,
#         f'collected {keyed_by}': collected_count,
#     }

#     # Write updated counts to cloud
#     with hl.hadoop_open(counts_file, 'w') as f:
#         json.dump(counts, f, indent=2)
        
#     print(f"Counts written to: {counts_file}")
