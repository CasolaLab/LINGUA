================================================================================
SCRIPT: step7_first_stage_validation.py
================================================================================
Description: A data-integrity checker that counts valid (non-dash) nucleotides across the Raw MAF, Filtered MAF, and fasta extracted files.

Usage:
  python3 step7_first_stage_validation.py <raw_mafextractor_dir> <step3_filtered_dir> <step4_output_dir> <reference_species>

Inputs:
  - <raw_mafextractor_dir>: Directory containing original unfiltered MAF files.
  - <step3_filtered_dir>: Directory of filtered MAF files.
  - <step4_output_dir>: Directory containing extracted .aln FASTAs.
  - <reference_species>: The reference species name prefix.
