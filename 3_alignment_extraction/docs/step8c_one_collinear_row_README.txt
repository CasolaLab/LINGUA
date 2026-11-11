================================================================================
SCRIPT: step8c_one_collinear_row.py  
================================================================================
Description: Gives every species and every ancestral node one row per stitched alignment, from one place in its
genome. Step 8 keys rows by species.contig, so a species whose exons were extracted from different contigs has
several partial rows, and exons on one contig are merged however far apart or on whichever strand they lie. For
each species, step 8c reads the source of every exon from its step 3 MAF (contig, strand, start, end; minus-strand
coordinates converted to the plus strand), groups the exons by contig and strand, chains them in reference exon
order, and breaks a chain when the species' coordinates stop moving in one direction or the chain's span would
exceed --maxspan. The chain with the most aligned bases is kept (ties: more exons, shorter span, contig name); every
other exon of that species becomes gaps. An exon whose blocks lie on both strands is split into one piece per
strand, and only the pieces on the kept strand stay. The reference row is never changed.

Usage:
  python3 step8c_one_collinear_row.py <workdir> <output_dir> <log_tsv> [--maxspan 50000] [--dry-run] [--mask PATH] [--tree <newick>]
  python3 step8c_one_collinear_row.py --self-test

Inputs:
  - <workdir>: Folder holding <reference_species>.bed (the only .bed there; the reference species is taken from its
    name), step3/ (step 3 MAFs), step4/ (per-exon alignments) and step8/ (stitched alignments and
    master_summary.tsv). Step 8 must have been run WITHOUT --one-row-per-species.
  - <output_dir>: New or empty folder for the rebuilt alignments (must differ from step8/).
  - <log_tsv>: One line per (alignment, species) changed: rows before, kept contig and strand, kept and gapped bases,
    and every gapped exon piece with its reason (other contig, opposite strand, beyond <maxspan> bp, not collinear).
  - --maxspan: Largest genomic span of a kept chain (default 50000, the jump limit step 3 uses within an exon).
  - --dry-run: Write only the log (and the mask when --mask is given), to test other --maxspan values.
  - --tree: (Default $TREE_FILE) Newick tree used for the log's type column: living (a tree tip), ancestor (a labeled
    internal node, so any naming scheme works: N12, fullTreeAnc235, PrimatesAnc045, ...) or not_in_tree; without a
    tree the type is 'unknown (no tree)'. The type is a label only; living and ancestral rows are processed the same.
  - --self-test: Run a synthetic case in a temporary folder and check the expected result.

Outputs:
  - <output_dir>/*.aln: Every alignment (unchanged ones copied byte for byte).
  - <output_dir>/master_summary.tsv: The step 8 summary with aligned_species updated to the new row count.
  - <output_dir>/step8c_missing_data.tsv: The missing-data mask: alignment, species, exon, ref_start, ref_end, reason.
    Lists every reference range (0-based half-open BED coordinates) in which a non-reference species has no aligned
    sequence after this step: exons it is absent from ("absent", "absent (padding only)") and every piece gapped here.
    Step 14 reads it so these gaps are not reported as deletions.

Checks (each stops the run): one row per species afterwards; species set unchanged; row lengths equal; every kept
base identical to step 8 at the same column; kept + gapped bases equal the step 8 bases of each species; reference
row unchanged; per-exon widths sum to the alignment length; a species on two contigs within one exon (step 3
prevents this) or bases without a step 3 MAF block.

