#!/usr/bin/env python3
# ================================================================================
# SCRIPT: step6_update_summary_files.py
# ================================================================================
# Description: 
#   Cross-references the FASTA extraction summary TSV against the master BED file
#   to correct any strand mislabeling introduced during initial extraction.
#
# Usage:
#   python3 step6_update_summary_files.py <step4_summary_tsv> <bed_file> <reference_species>
# ================================================================================

import sys
import csv
import os

def get_valid_chromosomes(bed_file):
    """Reads the BED file to learn valid chromosomes and build a coordinate-to-strand map."""
    valid_chroms = set()
    bed_map = {}
    try:
        with open(bed_file, 'r') as b:
            for line in b:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) >= 6:
                    chrom = parts[0].strip()
                    start = parts[1].strip()
                    end = parts[2].strip()
                    strand = parts[5].strip()
                    
                    valid_chroms.add(chrom)
                    bed_map[f"{chrom}_{start}_{end}"] = strand
    except Exception as e:
        print(f"Error reading BED file: {e}")
        sys.exit(1)
    return valid_chroms, bed_map

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 step6_update_summary_files.py <step4_summary_tsv> <bed_file> <reference_species>")
        sys.exit(1)

    tsv_path = sys.argv[1]
    bed_path = sys.argv[2]
    reference_species = sys.argv[3]

    if not os.path.exists(tsv_path):
        print(f"Error: TSV file '{tsv_path}' not found.")
        sys.exit(1)

    print(f"Loading BED file rules for {reference_species}...")
    valid_chroms, bed_map = get_valid_chromosomes(bed_path)
    print(f"Learned {len(valid_chroms)} unique chromosomes and {len(bed_map)} coordinate rules.")

    updated_rows = []
    fixed_count = 0
    total_rows = 0

    print("Scanning and updating TSV...")
    # Read the TSV
    with open(tsv_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames

        if not fieldnames or 'file_name' not in fieldnames or 'strand' not in fieldnames:
            print("Error: TSV must contain 'file_name' and 'strand' headers.")
            sys.exit(1)

        for row in reader:
            total_rows += 1
            base_name = row['file_name'].replace(".aln", "").replace(".txt", "")
            
            # 1. Strip species prefix
            base_name_no_sp = base_name
            if base_name.startswith(f"{reference_species}."):
                base_name_no_sp = base_name[len(f"{reference_species}."):]
            
            # 2. Smart match the chromosome
            matched_chrom = None
            for chrom in sorted(list(valid_chroms), key=len, reverse=True):
                if base_name_no_sp.startswith(f"{chrom}_"):
                    matched_chrom = chrom
                    break
            
            # 3. Extract coordinates and look up true strand
            if matched_chrom:
                coord_str = base_name_no_sp[len(f"{matched_chrom}_"):]
                try:
                    m_start, m_end = coord_str.split('_')
                    lookup_key = f"{matched_chrom}_{m_start}_{m_end}"
                    
                    correct_strand = bed_map.get(lookup_key)
                    
                    # Update if the bash script recorded it incorrectly
                    if correct_strand and row['strand'] != correct_strand:
                        row['strand'] = correct_strand
                        fixed_count += 1
                        
                except ValueError:
                    pass # Ignore if the filename is malformed
            
            updated_rows.append(row)

    # Write the updated TSV (Non-destructive: saves as _fixed.tsv)
    out_path = tsv_path.replace(".tsv", "_fixed.tsv")
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(updated_rows)

    print("-" * 50)
    print(f"Summary TSV Update Complete!")
    print(f"Total rows processed: {total_rows}")
    print(f"Rows updated with correct strand: {fixed_count}")
    print(f"Saved corrected summary to: {out_path}")

if __name__ == "__main__":
    main()