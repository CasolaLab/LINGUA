================================================================================
SCRIPT: step4_reverse_complement_add_missingbp.sh
================================================================================
Description: Converts the filtered MAF files into aligned FASTA (.aln), one per exon: concatenates each species' blocks in MAF order, pads a species that is absent from a block with dashes so all rows keep the same length, and reverse-complements every row when the reference exon is on the minus strand in the BED. Bases a species skips between blocks are unaligned sequence and are not added (no genome lookup; genome FASTAs and samtools are not used). Writes <output_dir>/<exon>/<exon>.aln and <output_dir>/extraction_summary.tsv (rows in sorted order). Requires GNU parallel.

Usage:
  ./step4_reverse_complement_add_missingbp.sh -m <step3_filtered_dir> -o <output_dir> -b <bed_file> [-t <threads>]

Inputs:
  - -m: Directory of filtered MAF files from step 3.
  - -o: Output directory for .aln files.
  - -b: BED file of target regions (used for strand lookup).
  - -t: (Optional) Number of parallel jobs (Default: 4).
