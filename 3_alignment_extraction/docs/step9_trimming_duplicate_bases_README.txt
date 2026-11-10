================================================================================
SCRIPT: step9_trimming_duplicate_bases.py  
================================================================================
Description: Removes from every stitched alignment the one extra reference base that step 1 extracted per exon, at its known position from the BED (last base of each exon block on plus-strand genes, first base on minus-strand genes; no extra base on an exon that ends at the end of its scaffold), and removes columns where the reference has a gap, so every species is projected onto the reference CDS. Because the extra base is located by position rather than by matching sequence, every species is cut at the same correct column even when the extra base equals its neighbor. A file is written only if the trimmed reference row equals its CDS in the FNA (N or IUPAC ambiguity codes may match a base they stand for); every other file is refused and listed with its reason (e.g. REFERENCE_MISMATCH_AFTER_TRIM, EXTRACTION_SHIFTED, NOT_IN_FNA) in <output_dir>/step9_trimming_report.tsv. The BED convention (0-based or 1-based) is detected from the data.

Usage:
  python3 step9_trimming_duplicate_bases.py <step8_stitched_dir> <output_dir> <reference_nucleotide_fna_file> <bed_file> <threads>

Inputs:
  - <step8_stitched_dir>: Stitched alignments from step 8c (run_stage3.sh passes step8c/; plain step 8 output also works). Searched recursively, so keep other .aln files out of it.
  - <output_dir>: Destination for the trimmed .aln files; must be new or empty.
  - <reference_nucleotide_fna_file>: FASTA of the reference CDS sequences, named by the BED transcript IDs.
  - <bed_file>: Reference BED file containing exon coordinates and transcript IDs.
  - <threads>: Number of threads to use.

Notes: before running step 9 on its own, run clean_nucleotide_fasta_headers.py on the reference nucleotide FASTA (run_stage3.sh does this).
