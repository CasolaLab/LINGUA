================================================================================
SCRIPT: step16_extract_start_and_stops.py
================================================================================
Description: A fast directory-wide parser that pulls the first and last 3 bases of sequences, standardizing gap removal and validating if they match known translation start/stop motifs.

Usage:
  python3 step16_extract_start_and_stops.py <step9_trimmed_dir> <output.tsv>

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <output.tsv>: Output summary TSV file name.
