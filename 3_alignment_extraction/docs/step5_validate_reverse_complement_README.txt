================================================================================
SCRIPT: step5_validate_reverse_complement.py
================================================================================
Description: Uses local sequence alignment (dynamic programming) to compare the original raw MAF sequence against the final .aln FASTA sequence to verify strand orientation.

Usage:
  python3 step5_validate_reverse_complement.py <step4_output_dir> <mafextractor_dir> <bed_file> <reference_species> <threads>

Inputs:
  - <step4_output_dir>: Directory containing extracted .aln FASTAs.
  - <mafextractor_dir>: Directory containing original unfiltered MAF files.
  - <bed_file>: Reference BED file with strand targets.
  - <reference_species>: Name of the reference species.
  - <threads>: Number of threads to use.
