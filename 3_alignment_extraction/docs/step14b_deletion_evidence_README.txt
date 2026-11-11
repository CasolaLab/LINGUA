================================================================================
SCRIPT: step14b_deletion_evidence.py  (runs after step 14; step 14's own files are unchanged)
================================================================================
Description: For every deletion step 14 reported, records where its gap came from and whether the species still has
its own sequence in that spot, in one table (one row per deletion, same order as step14/master_indels.tsv). Nothing
is removed or relabeled in step 14's output: the table is evidence for deciding which deletions to trust.

Usage:
  python3 step14b_deletion_evidence.py <step14_master_indels.tsv> <step1_dir> <step3_dir> <step9_dir> <bed_file> <reference_species> <output.tsv> [--tree <newick>]

Inputs:
  - <step14_master_indels.tsv>: step14/master_indels.tsv.
  - <step1_dir>, <step3_dir>: step 1 and step 3 MAFs (what each species had before and after step 3).
  - <step9_dir>: step 9 alignments (which relatives are present, for the sharing columns).
  - <bed_file>: Reference BED (0-based half-open), to turn CDS positions into genomic positions.
  - <reference_species>: The reference species name.
  - <output.tsv>: Must be outside step14/: stage 4 reads every .tsv there as indel calls.
  - --tree: (Default $TREE_FILE) Newick tree; gives species_type and the parent and sister for the sharing columns.

Columns (step 14's columns first, then):
  species_type            living (tree tip), ancestor (labeled internal node), not_in_tree, or NA without a tree
  gap_source              cactus_gap    every position is a '-' the aligner placed inside a block step 3 kept
                          no_alignment  some position had no block of this species at all, even before step 3
                          dropped_piece some position was covered only by a block step 3 dropped from this species
  bases_cactus_gap, bases_no_alignment, dropped_covers_bp
                          how many positions of the deletion fall in each kind of gap
  dropped_from            for dropped_piece: where the dropped block came from, compared with the blocks the species
                          kept in this exon: other chromosome / same chromosome, other strand / same chromosome and
                          strand (several joined by ';')
  dropped_distance_bp     for a dropped block on the kept chromosome: its distance from the kept blocks (0 = overlap)
  reference_gap_bp        reference bases between the kept blocks just before and after the gap
  species_bases_in_place  species bases between those same two kept blocks, in the species' own genome
  sides_touch             yes when species_bases_in_place is 0: the species sequence runs straight across the
                          missing stretch, which is what a real deletion looks like
  in_place_difference_bp  reference_gap_bp - species_bases_in_place: the likely real size of the change
                          (negative when the species has more bases there than the reference)
  in_place_frameshift     yes when in_place_difference_bp is not a multiple of 3
  exon_edge               start / end (the gap reaches the exon's first / last base, in gene direction), both,
                          none, or "spans N exons"
  cant_tell_reason        why the in-place check could not run: no kept piece before / after the gap in this exon,
                          no kept piece anywhere in this exon, gap crosses into the next exon, the two sides are on
                          different chromosomes or strands, or the species' two sides overlap
  shared_with_parent      yes / no: an overlapping deletion is reported for the tree parent in the same alignment;
                          no data (parent not in the alignment); no parent (tree root). Parents are reconstructed
                          partly from these same species, so this is weak evidence.
  shared_with_sister      the same for any other child of the parent (closer to independent evidence)
  evidence                aligner gap (kept as called)  gap_source cactus_gap
                          likely real deletion          sides_touch = yes
                          possible smaller frameshift   sides do not touch; in_place_difference_bp not a multiple of 3
                          probably not a frameshift     sides do not touch; in_place_difference_bp 0 or a multiple of 3
                          can't tell                    the in-place check could not run (see cant_tell_reason)

Checks (each stops the run): a deletion position where a kept block gives the species a base; a step 3 block not
found in step 1; rows that do not add up to step 14's deletions.

Notes:
  - The in-place check needs a kept block on both sides of the gap in the same exon; gaps at an exon edge or across
    exons are "can't tell". For dropped_piece rows, dropped_from, dropped_distance_bp and the sharing columns are
    still filled for those.
  - The species' actual bases in the gap are not shown: that would need each species' genome FASTA (ancestral
    genomes are not available as FASTA here).
  - Capra hircus test run (2026-11-10, 81,481 deletions): cactus_gap 78,808; no_alignment 2,335 (793 likely real,
    370 possible smaller frameshift, 193 probably not, 979 can't tell); dropped_piece 338 (1 likely real, 54 possible
    smaller frameshift, 79 probably not, 204 can't tell).
