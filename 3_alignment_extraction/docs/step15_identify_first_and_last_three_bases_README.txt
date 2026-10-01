================================================================================
SCRIPT: step15_identify_first_and_last_three_bases.sh
================================================================================
Description: Extracts the exact first and last 3 base pairs for a species from the multiple sequence alignment.

Usage:
  ./step15_identify_first_and_last_three_bases.sh <step9_trimmed_dir> <output_dir> <reference_species> <step11_summary_tsv>

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <output_dir>: Destination for the summary file.
  - <reference_species>: The reference species name(or other species name of interest).
  - <step11_summary_tsv>: Summary TSV file from stitched and trimmed exons in step 11.