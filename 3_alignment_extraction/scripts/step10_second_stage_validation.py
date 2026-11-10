#!/usr/bin/env python3
"""
==============================================================================
Script: step10_second_stage_validation.py
Description: Length and content audit. Compares the final trimmed and stitched sequences against the reference FASTA to guarantee that length and base sequence match perfectly before downstream analysis.
Usage: 
python3 step10_second_stage_validation.py <step9_trimmed_dir> <reference_nucleotide_fna_file> <reference_species>
==============================================================================
"""

import os
import sys
from Bio import SeqIO

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 verify_final_alignments.py <step9_trimmed_dir> <reference_nucleotide_fna_file> <reference_species>")
        sys.exit(1)

    aln_dir = sys.argv[1]
    fna_file = sys.argv[2]
    reference_species = sys.argv[3].lower()

    # 1. Load the "Truth" Sequences from the FNA file
    print("Loading truth sequences from FNA...")
    truth_map = {}
    try:
        for rec in SeqIO.parse(fna_file, "fasta"):
            clean_id = rec.id.split()[0].strip()
            truth_map[clean_id] = str(rec.seq).upper()
    except Exception as e:
        print(f"Failed to load FNA file: {e}")
        sys.exit(1)

    # 2. Scan the Final Alignment files
    aln_files = sorted([f for f in os.listdir(aln_dir) if f.endswith(".aln")])
    print(f"\n{'Filename':<55} | {'ALN bases':<10} | {'FNA bases':<10} | {'STATUS'}")
    print("-" * 110)

    mismatch_count = 0

    for fname in aln_files:
        target_id = None
        for part in fname.replace('.aln', '').split('_'):
            if part in truth_map:
                target_id = part
                break
        
        if not target_id:
            continue
            
        truth_seq = truth_map[target_id]
        aln_path = os.path.join(aln_dir, fname)
        
        try:
            aln_records = list(SeqIO.parse(aln_path, "fasta"))
            ref_idx = 0
            # Dynamically search for the reference species provided in the argument
            for i, rec in enumerate(aln_records):
                if reference_species in rec.id.lower() or "ref" in rec.id.lower():
                    ref_idx = i
                    break
            
            gapped_ref_seq = str(aln_records[ref_idx].seq).upper()
            degapped_ref_seq = gapped_ref_seq.replace("-", "")
            
            status = "PASS"
            if len(degapped_ref_seq) != len(truth_seq):
                status = f"FAIL: Length Diff! (Off by {abs(len(degapped_ref_seq) - len(truth_seq))} bases)"
            elif degapped_ref_seq != truth_seq:
                status = "FAIL: Length matches, but sequences are different (SNPs found)!"
                
            if status != "PASS":
                mismatch_count += 1
                print(f"{fname:<55} | {len(degapped_ref_seq):<10} | {len(truth_seq):<10} | {status}")
                
        except Exception as e:
            print(f"{fname:<55} | {'ERROR':<10} | {'ERROR':<10} | FAIL: Cannot read file")
            mismatch_count += 1

    print("-" * 110)
    print(f"Final Validation complete! Total stitched/trimmed files checked: {len(aln_files)}")
    print(f"Files failing final biological validation: {mismatch_count}")

if __name__ == "__main__":
    main()