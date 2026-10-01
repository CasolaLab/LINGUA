================================================================================
SCRIPT: step3_filtering_species.py
================================================================================
Description: An auditing tool that drops species sequences from MAF files if they exhibit unnatural genomic jumps (changing chromosomes or jumping >50kb on the same chromosome).

Usage:
  python3 step3_filtering_species.py -i <mafextractor_dir> -o <output_dir> -r <reference_species> -s <output_summary_csv_name> [-j <jump_threshold>] [-t <threads>]

Inputs:
  - -i / --mafextractor_dir: Directory with raw MAF files.
  - -o / --output: Directory for the cleansed MAF files.
  - -r / --ref: The reference species name to protect.
  - -s / --summary: Path to output the purge log CSV.
  - -j / --jump: (Optional) Max allowed base pair jump (Default: 50000).
  - -t / --threads: (Optional) Number of threads to use (Default: 8).