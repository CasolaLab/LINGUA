================================================================================
SCRIPT: step17_combine_all_species_info.py
================================================================================
Description: The final aggregation script. For every species (living or ancestral) in each step 9 alignment, writes the genomic span of the exons it kept: the smallest start and largest end of the step 3 MAF blocks that supplied its bases in that alignment. Only blocks on the species' kept contig (from its header in the step 9 alignment) are used; with --mask, blocks whose reference range lies inside a range step 8c masked for that species (exons gapped as other contig, opposite strand, beyond the span limit or not collinear) are skipped; blocks that cover only the extra base step 1 extracts past the exon end are skipped. The reference row is the BED span of the transcript. Each row is typed from the tree: living (a tip), ancestor (a labeled internal node, whatever its naming scheme), not_in_tree, or NA without a tree.

Usage:
  python3 step17_combine_all_species_info.py <step9_trimmed_dir> <step3_filtered_dir> <bed_file> <reference_species> <output.tsv> [--mask <step8c_missing_data.tsv>] [--tree <newick>]

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <step3_filtered_dir>: Directory of filtered MAF files from step 3.
  - <bed_file>: Reference BED file containing exon coordinates and Transcript IDs ('#' and header lines skipped).
  - <reference_species>: The reference species name.
  - <output.tsv>: Output summary TSV file name.
  - --mask: (Optional) step8c/step8c_missing_data.tsv. run_stage3.sh always passes it.
  - --tree: (Optional, default $TREE_FILE) Newick tree with the MAF species and ancestor names; run_stage3.sh passes it.

Output columns: file_name, species_name, chromosome, start, end (0-based half-open, plus strand), species_type.

Checks: with --mask the run stops if a block outside the mask lies on a contig other than the species' kept contig, or if a species' kept blocks lie on more than one contig (both mean step 8c and step 9 disagree). The run prints how many blocks were used and skipped, by reason.

Notes:
  - Alignment files are read in sorted order, so reruns write rows in the same order.
  - Without --mask, exons step 8c gapped on the kept contig cannot be recognized, and a species with several rows
    (step 8c not run) keeps a comma-joined list of its row contigs.
  - The transcript is matched to the file name by the longest BED transcript ID that prefixes it (before, the first
    match in BED order).
