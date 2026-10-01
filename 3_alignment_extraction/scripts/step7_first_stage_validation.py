#!/usr/bin/env python3
"""
==============================================================================
Script: step7_first_stage_validation.py
Description: End-to-end data integrity validation. Checks non-dash nucleotide counts across Raw MAF, Filtered MAF, and final extracted ALN to catch dropped or mangled sequences.
Usage: python3 step7_first_stage_validation.py <raw_mafextractor_dir> <step3_filtered_dir> <step4_output_dir> <reference_species>
==============================================================================
"""

import os
import sys

def count_maf_bases(filepath, reference_species_name):
    """Counts non-dash reference_specieserence bases in a MAF file."""
    count = 0
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("s ") and reference_species_name in line.split()[1]:
                    # Extract the sequence and count the actual characters (no dashes)
                    count += len(line.split()[-1].replace("-", ""))
    except Exception:
        pass
    return count

def count_aln_bases(filepath, reference_species_name):
    """Counts non-dash reference_specieserence bases in a FASTA/ALN file."""
    count = 0
    is_reference_species = False
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    is_reference_species = (reference_species_name in line)
                elif is_reference_species:
                    count += len(line.strip().replace("-", ""))
    except Exception:
        pass
    return count

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 stepwise_validation.py <raw_mafextractor_dir> <step3_filtered_dir> <step4_output_dir> <reference_species>")
        sys.exit(1)

    raw_mafextractor_dir, step3_filtered_dir, step4_output_dir, reference_species = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    
    print(f"{'Filename':<55} | {'RAW bases':<10} | {'FILT bases':<10} | {'ALN bases':<10} | {'STATUS'}")
    print("-" * 115)
    
    raw_files = sorted([f for f in os.listdir(raw_mafextractor_dir) if f.endswith('.txt')])
    
    mismatch_count = 0
    for f in raw_files:
        raw_path = os.path.join(raw_mafextractor_dir, f)
        filt_path = os.path.join(step3_filtered_dir, f)
        
        base_name = f.replace(".txt", "")
        aln_path = os.path.join(step4_output_dir, base_name, base_name + ".aln")

        raw_cnt = count_maf_bases(raw_path, reference_species)
        filt_cnt = count_maf_bases(filt_path, reference_species)
        aln_cnt = count_aln_bases(aln_path, reference_species) if os.path.exists(aln_path) else "MISSING"

        status = "PASS"
        if raw_cnt != filt_cnt:
            status = "FAIL: Filter Script Dropped Bases!"
        elif str(aln_cnt) != "MISSING" and filt_cnt != aln_cnt:
            status = "FAIL: FASTA Extractor Dropped Bases!"
        
        if status != "PASS" or str(aln_cnt) == "MISSING":
            mismatch_count += 1
            print(f"{f:<55} | {raw_cnt:<10} | {filt_cnt:<10} | {aln_cnt:<10} | {status}")

    print("-" * 115)
    print(f"Stepwise Validation complete! Total files checked: {len(raw_files)}")
    print(f"Files failing stepwise validation: {mismatch_count}")

if __name__ == "__main__":
    main()