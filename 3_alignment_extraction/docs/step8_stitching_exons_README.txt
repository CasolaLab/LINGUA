================================================================================
SCRIPT: step8_stitching_exons.py
================================================================================

Description: Stitches the per-exon alignments of step 4 into one alignment per transcript, in transcript order (exons by chromStart, descending on the minus strand), and writes master_summary.tsv. A species missing from an exon gets dashes of that exon's length. Rows are keyed by the full header (species.chromosome), so a species whose exons come from different chromosomes gets one partial row per chromosome. run_stage3.sh runs this step without options and then step 8c (step8c_one_collinear_row.py), which keeps one row per species from its dominant collinear chain and needs every species.chromosome row. With --one-row-per-species (not used by run_stage3.sh; do not feed its output to step 8c) each species keeps only its row with the most aligned bases (ties: the first seen), without checking strand, exon order or distance; the exons of its other rows become gaps, every species that had more than one row is logged in <output_dir>/one_row_per_species_log.tsv, and the dropped rows are written, in the same columns, to <output_dir>_one_row_dropped/ for inspection. A header line in the BED is detected and skipped (a BED without a header is read completely).

Usage:
  TREE_FILE=/path/to/tree.nwk python3 step8_stitching_exons.py <step4_output_dir> <bed_file> <output_dir> <reference_species> [--one-row-per-species]

Inputs:
  - TREE_FILE (environment): Newick tree whose tip names are the MAF species names; used for distant_species and distant_node.
  - <step4_output_dir>: Directory containing the per-exon .aln files from step 4.
  - <bed_file>: Reference BED file containing exon coordinates and transcript IDs.
  - <output_dir>: Destination for the stitched .aln files and master_summary.tsv.
  - <reference_species>: The reference species name.

master_summary.tsv columns: file_name, gff_transcriptID, chromosome (BED name with a leading Chr written as chr), strand, length (sum of end - start + 1, i.e. the stitched reference length before step 9 removes one extra base per exon), number_of_exons, aligned_species, distant_species, distant_node, start, end.
