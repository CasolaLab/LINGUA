================================================================================
SCRIPT: step11_update_stitched_exons_summary.py
================================================================================
Description: Adds to the stitched-alignment master_summary.tsv (step8c/master_summary.tsv in run_stage3.sh) the CDS length from the FNA (expected_fna), the reference length in the step 9 alignment with gaps removed (actual_nt) and a validation status (PASS, PASS_IN_FRAME_SHORT, FAIL_FRAME_SHIFT, FAIL_MISSING_FNA_TRUTH, FAIL_MISSING_ALN_FILE, ...), and writes <step8_summary>_final_validation.tsv. Steps 12-15 read this file; steps 12-14 take reference_length from its actual_nt column.

Usage:
  python3 step11_update_stitched_exons_summary.py <step8_output_summary_tsv> <step9_trimmed_dir> <bed_file> <reference_nucleotide_fna_file> <reference_species>

Inputs:
  - <step8_output_summary_tsv>: master_summary.tsv from step 8c (run_stage3.sh) or step 8; the output is written next to it.
  - <step9_trimmed_dir>: Trimmed alignments from step 9.
  - <bed_file>: Reference BED file containing exon coordinates and transcript IDs.
  - <reference_nucleotide_fna_file>: The same (header-cleaned) CDS FASTA given to step 9; with another file every row fails with FAIL_MISSING_FNA_TRUTH and actual_nt = 0.
  - <reference_species>: The reference species name prefix.