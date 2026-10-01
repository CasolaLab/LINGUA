================================================================================
SCRIPT: step4_reverse_complement_add_missingbp.sh
================================================================================
Description: Converts filtered MAF files into aligned FASTA format. It dynamically fills authorized gaps (using the CSV from Step 2) and reverse-complements based on reference strand data.

Usage:
  ./step4_reverse_complement_add_missingbp.sh -m <step3_filtered_dir> -g <genome_fasta_dir> -o <output_dir> -b <bed_file> -c <step2_combined_summary_file> [-t <threads>]

Inputs:
  - -m: Directory of filtered MAF files.
  - -g: Directory containing reference genome fasta files.
  - -o: Output directory for .aln files.
  - -b: BED file of target regions (used for strand lookup).
  - -c: The combined_all_species.csv file generated in Step 2.
