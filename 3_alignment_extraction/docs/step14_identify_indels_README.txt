================================================================================
SCRIPT: step14_identify_indels.sh
================================================================================
Description: Identifies biological insertions and deletions. Crucially, it isolates specific 3-bp insertions that themselves encode a premature stop codon.

Usage:
  ./step14_identify_indels.sh <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv>

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <output_dir>: Destination for the indel summary files.
  - <bed_file>: Reference BED file containing exon coordinates and Transcript IDs.
  - <reference_species>: The reference species name prefix.
  - <step11_summary_tsv>: Summary TSV file from stitched and trimmed exons in step 11.
