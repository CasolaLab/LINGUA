================================================================================
SCRIPT: step14_identify_indels.sh  (LINGUA masked version, with sorted input and a masked-deletions table)
================================================================================
Description: For every non-reference species in each step 9 alignment, finds runs of columns where the species has a gap and the reference a base (deletion) or the species a base and the reference a gap (insertion), and reports a run if its length is not a multiple of 3 (frameshift), or if it is a 3-bp insertion whose bases are a stop codon. With the step 8c missing-data mask, a deletion lying wholly inside a masked range of that species is missing data and is not reported, and a deletion overlapping a masked range is trimmed to its part(s) outside it, each reported only if its own length is not a multiple of 3 (remainders of 1 bp at a mask boundary are dropped). Transcript ID, chromosome and strand come from the step 11 summary; reference_length is the step 11 actual_nt column; cds_length is that species' own length with gaps removed. Writes one <file>_indels.tsv per alignment and master_indels.tsv; alignment files are read in sorted order.

Usage:
  ./step14_identify_indels.sh <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv> [missing_data_mask_tsv [masked_deletions_tsv]]

Inputs:
  - <step9_trimmed_dir>: Trimmed alignments from step 9.
  - <output_dir>: Destination for the indel tables.
  - <bed_file>: Reference BED file, 0-based half-open (used to convert the mask to CDS positions).
  - <reference_species>: The reference species name prefix.
  - <step11_summary_tsv>: The step 11 *_final_validation.tsv (must contain actual_nt).
  - [missing_data_mask_tsv]: (Optional) step8c/step8c_missing_data.tsv. Without it every gap run is tested as before.
  - [masked_deletions_tsv]: (Optional, needs the mask) Every deletion the mask removed (action masked_out) or
    trimmed (action trimmed), with the remainders reported, the in-frame remainders dropped and the 1-bp boundary
    remainders dropped. Must be outside <output_dir>: stage 4 reads every .tsv in the step 14 folder as indel calls.

Notes:
  - run_stage3.sh passes step8c/step8c_missing_data.tsv and writes step14_masked_deletions.tsv in the workdir.
  - Insertions are not reported in practice: the alignments reaching step 14 contain no column where the reference
    has a gap. In the Capra hircus test run (2026-11-10) the reference rows had no gaps already in the step 1 MAFs
    (13,071 reference rows, 0 gap characters), and step 9 removes any such column. Frameshift insertions and 3-bp
    stop-codon insertions relative to the reference are therefore not detected by stage 3 as it stands.
  - Gaps from a species missing part of an exon (absent from some MAF blocks, padded with dashes by step 4) are not
    in the mask and are still tested as deletions.
