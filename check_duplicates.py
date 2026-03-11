#!/usr/bin/env python3
import gzip

print("Processing file in chunks, looking for duplicates with DIFFERENT uniprot_ids...")

# Track statistics
total_lines = 0
duplicates_found = []
max_duplicates = 5
total_duplicates_checked = 0

# Store current key and its entries
current_key = None
current_entries = []

with gzip.open('/Users/tk508/Work/GenGym/genetics-gym-final/data/VSM/AlphaMissense_hg38.tsv.gz', 'rt') as f:
    for line in f:
        # Skip comment lines
        if line.startswith('#'):
            continue
        
        total_lines += 1
        if total_lines % 1000000 == 0:
            print(f"Processed {total_lines:,} lines, found {len(duplicates_found)} with different uniprot_ids (checked {total_duplicates_checked} total duplicates)...")
        
        parts = line.strip().split('\t')
        if len(parts) < 10:
            continue
            
        chrom, pos, ref, alt, genome, uniprot_id, transcript_id, protein_variant, am_pathogenicity, am_class = parts
        
        # Create key
        key = f"{chrom}-{pos}-{ref}-{alt}"
        
        # If we encounter a new key
        if key != current_key:
            # Check if previous key had duplicates
            if current_key is not None and len(current_entries) > 1:
                total_duplicates_checked += 1
                
                # Check if uniprot_ids are different
                uniprot_ids = set(e['uniprot_id'] for e in current_entries)
                if len(uniprot_ids) > 1:
                    duplicates_found.append((current_key, current_entries))
                    print(f"Found duplicate #{len(duplicates_found)} with DIFFERENT uniprot_ids: {current_key} ({len(current_entries)} entries, {len(uniprot_ids)} different proteins)")
                    
                    # Stop if we've found enough duplicates
                    if len(duplicates_found) >= max_duplicates:
                        print(f"\nReached {max_duplicates} duplicates with different uniprot_ids, stopping...")
                        break
            
            # Start tracking new key
            current_key = key
            current_entries = []
        
        # Add entry to current key
        current_entries.append({
            'uniprot_id': uniprot_id,
            'transcript_id': transcript_id,
            'protein_variant': protein_variant,
            'am_pathogenicity': am_pathogenicity,
            'am_class': am_class
        })

    # Check last key if we didn't hit the limit
    if len(duplicates_found) < max_duplicates and current_key is not None and len(current_entries) > 1:
        uniprot_ids = set(e['uniprot_id'] for e in current_entries)
        if len(uniprot_ids) > 1:
            total_duplicates_checked += 1
            duplicates_found.append((current_key, current_entries))

print(f"\n=== SUMMARY ===")
print(f"Total lines processed: {total_lines:,}")
print(f"Duplicate keys found: {len(duplicates_found)}")

if duplicates_found:
    print("\n=== ANALYZING DUPLICATES ===")
    
    different_uniprot = 0
    different_transcript = 0
    different_transcript_same_uniprot = 0
    different_am_score = 0
    
    for key, entries in duplicates_found:
        uniprot_ids = set(e['uniprot_id'] for e in entries)
        transcript_ids = set(e['transcript_id'] for e in entries)
        am_scores = set(e['am_pathogenicity'] for e in entries)
        
        has_diff_uniprot = len(uniprot_ids) > 1
        has_diff_transcript = len(transcript_ids) > 1
        has_diff_am_score = len(am_scores) > 1
        
        if has_diff_uniprot:
            different_uniprot += 1
        if has_diff_transcript:
            different_transcript += 1
        if has_diff_transcript and not has_diff_uniprot:
            different_transcript_same_uniprot += 1
        if has_diff_am_score:
            different_am_score += 1
    
    print(f"\n=== RESULTS ===")
    print(f"1) Different uniprot_id: {different_uniprot}/{len(duplicates_found)} keys")
    print(f"2) Different transcript_id: {different_transcript}/{len(duplicates_found)} keys")
    print(f"3) Different transcript_id but same uniprot_id: {different_transcript_same_uniprot}/{len(duplicates_found)} keys")
    print(f"4) Different am_pathogenicity scores: {different_am_score}/{len(duplicates_found)} keys")
    
    print(f"\n=== DETAILED EXAMPLES ===")
    for key, entries in duplicates_found:
        print(f"\nKey: {key} (count: {len(entries)})")
        for i, entry in enumerate(entries, 1):
            print(f"  Entry {i}:")
            print(f"    uniprot_id: {entry['uniprot_id']}")
            print(f"    transcript_id: {entry['transcript_id']}")
            print(f"    protein_variant: {entry['protein_variant']}")
            print(f"    am_pathogenicity: {entry['am_pathogenicity']}")
            print(f"    am_class: {entry['am_class']}")
else:
    print("\nNo duplicate keys found!")
