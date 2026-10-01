================================================================================
SCRIPT: step10_second_stage_validation.py
================================================================================
Description: An audit script that checks the final trimmed/stitched sequences against the reference FASTA to guarantee that length and base content match perfectly.

Usage:
  python3 step10_second_stage_validation.py <step9_trimmed_dir> <reference_nucleotide_fna_file> <reference_species>

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <reference_nucleotide_fna_file>: Fasta file of true reference nucleotide sequences.
  - <reference_species>: The reference species name prefix.
