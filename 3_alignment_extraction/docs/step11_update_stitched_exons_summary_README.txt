================================================================================
SCRIPT: step11_update_stitched_exons_summary.py
================================================================================
Description: Calculates expected lengths from the BED file and compares them to actual nucleotide lengths in the alignment.

Usage:
  python3 step11_update_stitched_exons_summary.py <step8_output_summary_tsv> <step9_trimmed_dir> <bed_file> <reference_nucleotide_fna_file> <reference_species>

Inputs:
  - <step8_output_summary_tsv>: Summary TSV file from stitched exons in step 8.
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <bed_file>: Reference BED file containing exon coordinates and Transcript IDs.
  - <reference_nucleotide_fna_file>: Fasta file of true reference nucleotide sequences.
  - <reference_species>: The reference species name prefix.

