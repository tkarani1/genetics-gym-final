import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import FORMATTED_VSM_HT_PATHS, WRITE_VSM_LINKER_TABLES_PATH, LINKER_PATHS

linker = LINKER_PATHS['MISSENSE_ENST_TRANSCRIPT_AA']

# Read and count raw data
am1 = hl.read_table(FORMATTED_VSM_HT_PATHS['AM_ISOFORMS_PATH'])
am1_raw_count = am1.count()
am1 = am1.key_by('enst', 'aa_pos', 'aa_ref', 'aa_alt')
am1 = am1.select(am1.AM)

am2 = hl.read_table(FORMATTED_VSM_HT_PATHS['AM_CANONICAL_PATH'])
am2_raw_count = am2.count()
am2 = am2.key_by('locus', 'alleles', 'enst')


linker_ht = hl.read_table(linker)
linker_ht = linker_ht.key_by('enst', 'aa_pos', 'aa_ref', 'aa_alt')
ht = linker_ht.join(am1, how='right')

ht = ht.key_by('locus', 'alleles', 'enst')
ht = ht.join(am2, how='right')
ht = ht.checkpoint(f'{WRITE_VSM_LINKER_TABLES_PATH}/temp/vsm_am_temp.ht', overwrite=True)

ht = ht.annotate(
    AM = hl.coalesce(ht.AM_score, ht.AM) # joining the canonical one first
)
ht.write(f'{WRITE_VSM_LINKER_TABLES_PATH}/linker_AM.ht')

am2_collected = am2.collect_by_key()
am2_collected_count = am2_collected.count()
am2_key2 = am2.key_by('locus', 'alleles')
am2_key2_collected = am2_key2.collect_by_key()
am2_key2_collected_count = am2_key2_collected.count()

# Write AM counts to JSON
am_counts_file = f'{WRITE_VSM_LINKER_TABLES_PATH}/AM_counts.json'
# Load existing counts if file exists
if hl.hadoop_exists(am_counts_file):
    with hl.hadoop_open(am_counts_file, 'r') as f:
        counts = json.load(f)
else:
    counts = {}

# Update with AM counts
counts['AM'] = {
    'raw_isoforms': am1_raw_count,
    'raw_canonical': am2_raw_count,
    'collected_canonical_SNP_enst': am2_collected_count,
    'collected_canonical_SNP_only': am2_key2_collected_count,
}

# Write updated counts to cloud
with hl.hadoop_open(am_counts_file, 'w') as f:
    json.dump(counts, f, indent=2)
    
print(f"Counts written to: {am_counts_file}")







