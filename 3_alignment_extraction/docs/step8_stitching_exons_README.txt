--------------------------------------------------------------------------------
8. step8_stitching_exons.sh
--------------------------------------------------------------------------------
Description: Stitches individual exon alignment blocks into full continuous transcript models based on reference coordinates.

Usage:
  export TREE_FILE=/path/to/tree.nwk
  export REF_SPECIES=reference_species_name
  export SPECIES_NAME="reference_species_name"
  ./step8_stitching_exons.sh <step4_output_dir> <bed_file> <output_dir> <reference_species>

Inputs:
  - <step4_output_dir>: Directory containing the individual .aln files.
  - <bed_file>: Reference BED file containing exon coordinates and Transcript IDs.
  - <output_dir>: Destination for the stitched .aln files.
  - <reference_species>: The reference species name.