#!/usr/bin/env python3
"""
==============================================================================
Script: step11_update_stitched_exons_summary.py
Description: Sequence integrity auditor. Calculates expected lengths from BED and actual nucleotide lengths to flag frameshifts.
  expected_bed = sum of the transcript's BED exon lengths (convention detected as in step 9); expected_fna = CDS
  length in the FNA; actual_nt = reference length in the step 9 alignment (gaps removed); validation compares
  actual_nt with expected_fna; bed_fna_check = MATCH, MISMATCH(bed=..,fna=..) or NA (BED vs FNA).
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
    
    # 2. Aggregate BED exon lengths per transcript
    print(f"Reading BED and summing exons: {bed_file}")
    bed = pd.read_csv(bed_file, sep='\t', comment='#', header=None, dtype=str)
    # a header line without '#' is skipped: keep only rows whose start and end are numbers
    numeric = bed[1].str.strip().str.fullmatch(r"\d+") & bed[2].str.strip().str.fullmatch(r"\d+")
    bed = bed[numeric.fillna(False)].copy()
    bed[1] = bed[1].astype(int)
    bed[2] = bed[2].astype(int)
    # Group by the transcript ID: column 7 (index 6, ProteinID) as in steps 8-14. The former code grouped by
    # column 4 (index 3, the per-exon name), which never matches a transcript ID, so expected_bed was always N/A.
    id_col = 6 if 6 in bed.columns else 3
    bed['tid'] = bed[id_col].astype(str).str.strip()

    # 3. Load Truth Sequence lengths from FNA
    print(f"Reading Reference FNA: {reference_nucleotide_fna_file}")
    truth_map = {}
    for record in SeqIO.parse(reference_nucleotide_fna_file, "fasta"):
        clean_id = record.id.split()[0].strip()
        seq_str = str(record.seq).upper().replace("-", "")
        truth_map[clean_id] = len(seq_str)

    # BED convention, decided as step 9 does: 0-based half-open (exon length = end - start) or 1-based inclusive
    # (end - start + 1), whichever makes the BED exon lengths equal the FNA length for more transcripts.
    sums = {conv: (bed[2] - bed[1] + extra).groupby(bed['tid']).sum().to_dict()
            for conv, extra in (("0-based", 0), ("1-based", 1))}
    fits = {conv: sum(1 for t, L in m.items() if truth_map.get(t) == L) for conv, m in sums.items()}
    convention = "1-based" if fits["1-based"] > fits["0-based"] else "0-based"
    expected_map = {t: int(L) for t, L in sums[convention].items()}
    print(f"BED convention: {convention} (transcripts whose exon lengths match the FNA: "
          f"0-based {fits['0-based']}, 1-based {fits['1-based']})")

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

        # BED vs FNA: does the annotation's exon structure describe the same CDS length as the FNA?
        # (validation itself compares the alignment with the FNA only; this catches a BED/FNA mismatch)
        if expected_bed == "N/A" or expected_truth == 0:
            bed_fna_check = "NA"
        elif expected_bed == expected_truth:
            bed_fna_check = "MATCH"
        else:
            bed_fna_check = f"MISMATCH(bed={expected_bed},fna={expected_truth})"

        # Write to our row
        row['expected_bed'] = expected_bed
        row['expected_fna'] = expected_truth
        row['actual_nt'] = actual_len
        row['validation'] = status
        row['bed_fna_check'] = bed_fna_check
        results.append(row)

    # 5. Save output as TSV
    out_path = step8_output_summary_tsv.replace(".tsv", "_final_validation.tsv")
    final_df = pd.DataFrame(results)
    final_df.to_csv(out_path, sep='\t', index=False)
    
    # Print summary statistics
    pass_count = sum(1 for r in results if "PASS" in r['validation'])
    checks = {k: sum(1 for r in results if r['bed_fna_check'].split("(")[0] == k) for k in ("MATCH", "MISMATCH", "NA")}
    print("-" * 50)
    print(f"Validation complete! {pass_count} / {len(results)} files passed perfectly.")
    print(f"BED exon lengths vs FNA ({convention}): match {checks['MATCH']}, mismatch {checks['MISMATCH']}, "
          f"not comparable {checks['NA']}")
    print(f"File saved: {out_path}")

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python3 step11_update_stitched_exons_summary.py <step8_output_summary_tsv> <step9_trimmed_dir> <bed_file> <reference_nucleotide_fna_file> <reference_species>")
        sys.exit(1)
        
    update_and_validate(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])