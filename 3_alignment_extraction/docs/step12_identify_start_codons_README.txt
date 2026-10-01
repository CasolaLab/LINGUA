================================================================================
SCRIPT: step12_identify_start_codons.sh
================================================================================
Description: Scans the stitched and trimmed multi-species alignment to dynamically locate the biological start codon (ATG) in frame, outputting coordinates mapped back to the reference.

Usage:
  ./step12_identify_start_codons.sh <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv>

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <output_dir>: Destination for the start codon summary files.
  - <bed_file>: Reference BED file containing exon coordinates and Transcript IDs.
  - <reference_species>: The reference species name prefix.
  - <step11_summary_tsv>: Summary TSV file from stitched and trimmed exons in step 11.
