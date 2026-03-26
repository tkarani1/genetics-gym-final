import hail as hl
import json
import os
hl.init(worker_memory="highmem", driver_memory='highmem') 
from resources.paths import FORMATTED_VSM_HT_PATHS, WRITE_VSM_LINKER_TABLES_BASE, LINKER_PATHS, VSM_TABLE_NAMES
from resources.functions import write_tsv_bgz_from_ht

METHOD = 'AM'
linker = LINKER_PATHS['MISSENSE_ENST_TRANSCRIPT']

# Read and count raw data
am_iso = hl.read_table(FORMATTED_VSM_HT_PATHS['AM_ISOFORMS'])
am_iso = am_iso.repartition(500)
am_iso = am_iso.key_by('locus', 'alleles', 'enst')
am_iso = am_iso.select(am_iso.AM_iso_score)

am_canon = hl.read_table(FORMATTED_VSM_HT_PATHS['AM_CANONICAL'])
am_canon = am_canon.repartition(500)
am_canon = am_canon.key_by('locus', 'alleles', 'enst')
am_canon = am_canon.select(am_canon.AM_canon)


linker_ht = hl.read_table(linker)
linker_ht = linker_ht.repartition(500)
linker_ht = linker_ht.key_by('locus', 'alleles', 'enst')

# join both

# ht = linker_ht.join(am_iso, how='right')
# ht = ht.key_by('locus', 'alleles', 'enst')
# ht = ht.join(am_canon, how='right')
# ht = ht.checkpoint(f'{WRITE_VSM_LINKER_TABLES_BASE}/temp/vsm_am_temp.ht', overwrite=True)
# ht = ht.annotate(
#     AM = hl.coalesce(ht.AM_canon, ht.AM_iso_score) # joining the canonical one first
# )
# ht.write(f'{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[METHOD]}_both_joined.ht')
# write_tsv_bgz_from_ht(ht, f'{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[METHOD]}_both_joined.tsv.bgz')


# Separate isoforms and canonical
ht1 = linker_ht.join(am_iso, how='right')
ht1 = ht1.rename({'AM_iso_score': 'AM'})
# ht1.write(f'{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[METHOD]}_iso.ht')
# write_tsv_bgz_from_ht(ht, f'{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[METHOD]}_iso.tsv.bgz')

ht2 = linker_ht.join(am_canon, how='right')
ht2 = ht2.rename({'AM_canon': 'AM'})
# ht2.write(f'{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[METHOD]}_canon.ht')
# write_tsv_bgz_from_ht(ht, f'{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[METHOD]}_canon.tsv.bgz')

# concat
ht1 = ht1.repartition(500)
ht2 = ht2.repartition(500)
ht = ht2.union(ht1)
ht.write(f'{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[METHOD]}_both_concat.ht', overwrite=True)
write_tsv_bgz_from_ht(ht, f'{WRITE_VSM_LINKER_TABLES_BASE}{VSM_TABLE_NAMES[METHOD]}_both_concat.tsv.bgz')


