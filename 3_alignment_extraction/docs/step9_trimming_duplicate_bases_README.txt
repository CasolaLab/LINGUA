================================================================================
SCRIPT: step9_trimming_duplicate_bases.py
================================================================================
Description: Trims extraneous padding and intron bleed-over from the mafextractions. It compares the alignment to the true Reference FNA sequence, strictly forcing the sequence to adhere to the correct biological reading frame.
**Before running step9, run clean_nucleotide_fasta_headers.py on the reference nucleotide fasta file first**

Usage:
  python3 step9_trimming_duplicate_bases.py <step8_stitched_dir> <output_dir> <reference_nucleotide_fna_file> <bed_file> <threads>

Inputs:
  - <step8_stitched_dir>: Stitched alignments from step 8.
  - <output_dir>: Destination for the trimmed .aln files.
  - <reference_nucleotide_fna_file>: Fasta file of true reference nucleotide sequences.
  - <bed_file>: Reference BED file containing exon coordinates and Transcript IDs.
  - <threads>: Number of threads to use.
