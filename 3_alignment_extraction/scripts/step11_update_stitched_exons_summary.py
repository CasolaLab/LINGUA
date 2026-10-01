#!/usr/bin/env python3
"""
==============================================================================
Script: step11_update_stitched_exons_summary.py
Description: Sequence integrity auditor. Calculates expected lengths from BED and actual nucleotide lengths to flag frameshifts.
Usage: python3 step11_update_stitched_exons_summary.py <step8_output_summary_tsv> <step9_trimmed_dir> <bed_file> <reference_nucleotide_fna_file> <reference_species>
==============================================================================
"""

import pandas as pd
from Bio import SeqIO
import os
import sys

def update_and_validate(step8_output_summary_tsv, step9_trimmed_dir, bed_file, reference_nucleotide_fna_file, reference_species):
    # 1. Load the Summary TSV
    print(f"Reading summary: {step8_output_summary_tsv}")
    df = pd.read_csv(step8_output_summary_tsv, sep='\t')
    
    # 2. Aggregate BED lengths safely
    print(f"Reading BED and summing exons: {bed_file}")
    bed = pd.read_csv(bed_file, sep='\t', comment='#', header=None)
    bed['length'] = bed[2] - bed[1]
    
    # Safely find the ID column (Usually index 3 in standard BED, or 6 if custom BED12)
    id_col = 3 if 3 in bed.columns else (6 if 6 in bed.columns else bed.columns[-2])
    expected_map = bed.groupby(id_col)['length'].sum().to_dict()
    expected_map = {str(k).strip(): v for k, v in expected_map.items()}

    # 3. Load Truth Sequence lengths from FNA
    print(f"Reading Reference FNA: {reference_nucleotide_fna_file}")
    truth_map = {}
    for record in SeqIO.parse(reference_nucleotide_fna_file, "fasta"):
        clean_id = record.id.split()[0].strip()
        seq_str = str(record.seq).upper().replace("-", "")
        truth_map[clean_id] = len(seq_str)

    # 4. Iterate over the summary and validate
    results = []
    print("Validating stitched alignments against FNA truth...")
    
    for index, row in df.iterrows():
        transcript_id = str(row['gff_transcriptID']).strip()
        fname = str(row['file_name']).strip()
        aln_path = os.path.join(step9_trimmed_dir, fname + ".aln")

        expected_bed = expected_map.get(transcript_id, "N/A")
        expected_truth = truth_map.get(transcript_id, 0)
        actual_len = 0
        status = "PENDING"

        if expected_truth == 0:
            status = "FAIL_MISSING_FNA_TRUTH"
        elif not os.path.exists(aln_path):
            status = "FAIL_MISSING_ALN_FILE"
        else:
            try:
                # Parse the final alignment file to find the species sequence length
                found_species = False
                for rec in SeqIO.parse(aln_path, "fasta"):
                    if reference_species.lower() in rec.id.lower():
                        actual_len = len(str(rec.seq).replace("-", ""))
                        found_species = True
                        break

                if not found_species:
                    status = "FAIL_SPECIES_NOT_IN_ALN"
                else:
                    diff_truth = actual_len - expected_truth
                    is_in_frame = (actual_len % 3 == 0)

                    if diff_truth == 0:
                        status = "PASS"
                    elif actual_len > 0 and is_in_frame:
                        # Shorter but frame is preserved
                        status = f"PASS_IN_FRAME_SHORT({diff_truth})"
                    elif not is_in_frame:
                        # Critical Error: Broke the 3-bp codon frame
                        status = f"FAIL_FRAME_SHIFT({actual_len % 3})"
                    else:
                        status = f"FAIL({diff_truth})"

            except Exception as e:
                status = "ERROR_PARSING_ALN"

        # Write to our row
        row['expected_bed'] = expected_bed
        row['expected_fna'] = expected_truth
        row['actual_nt'] = actual_len
        row['validation'] = status
        results.append(row)

    # 5. Save output as TSV
    out_path = step8_output_summary_tsv.replace(".tsv", "_final_validation.tsv")
    final_df = pd.DataFrame(results)
    final_df.to_csv(out_path, sep='\t', index=False)
    
    # Print summary statistics
    pass_count = sum(1 for r in results if "PASS" in r['validation'])
    print("-" * 50)
    print(f"Validation complete! {pass_count} / {len(results)} files passed perfectly.")
    print(f"File saved: {out_path}")

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python3 step11_update_stitched_exons_summary.py <step8_output_summary_tsv> <step9_trimmed_dir> <bed_file> <reference_nucleotide_fna_file> <reference_species>")
        sys.exit(1)
        
    update_and_validate(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])