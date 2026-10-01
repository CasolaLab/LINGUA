================================================================================
SCRIPT: step17_combine_all_species_info.py
================================================================================
Description: The final aggregation script. It compiles the mapped genomic coordinates of the functional exons for ALL species relative to the reference sequence into a master TSV.

Usage:
  python3 step17_combine_all_species_info.py <step9_trimmed_dir> <step3_filtered_dir> <bed_file> <reference_species> <output.tsv>

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <step3_filtered_dir>: Directory of filtered MAF files from step 3.
  - <bed_file>: Reference BED file containing exon coordinates and Transcript IDs.
  - <reference_species>: The reference species name(or other species name of interest).
  - <output.tsv>: Output summary TSV file name.
