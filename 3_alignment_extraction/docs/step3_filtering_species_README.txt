================================================================================
SCRIPT: step3_filtering_species.py  (strict step 3; run_stage3.sh uses it with STEP3=strict)
================================================================================
Description: Removes a species from an exon's MAF file if, between consecutive blocks of that species, it changes chromosome or jumps more than the jump threshold (|start - previous end| > 50 kb by default). The species is removed from every block of that exon; the reference species is never removed, and blocks without the reference are not written. Strand is not considered. Files and summary rows are processed and written in sorted order, so reruns give identical files.

Usage:
  python3 step3_filtering_species.py -i <mafextractor_dir> -o <output_dir> -r <reference_species> -s <output_summary_csv_name> [-j <jump_threshold>] [-t <threads>]

Inputs:
  - -i / --mafextractor_dir: Directory with step 1 MAF files.
  - -o / --output_dir: Directory for the filtered MAF files.
  - -r / --ref: The reference species name to protect.
  - -s / --summary: Path to output the purge log CSV (file_name, species_purged, status).
  - -j / --jump: (Optional) Max allowed base pair jump (Default: 50000).
  - -t / --threads: (Optional) Number of threads to use (Default: 8).

Note: the default step 3 in run_stage3.sh is step3_dominant_chain.py (see its README).
