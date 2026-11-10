================================================================================
SCRIPT: step12_identify_start_codons.sh
================================================================================
Description: For every non-reference species in each step 9 alignment, removes that species' gaps and reads its bases in groups of three from its own first aligned base, reporting every codon that is ATG. Positions are reported in reference coordinates: the reference position of the alignment column holding the codon's first and third base. Transcript ID, chromosome and strand come from the step 11 summary (matched by file name); reference_length is the step 11 actual_nt column (found by name); cds_length is that species' own length with gaps removed. Writes one <file>_starts.tsv per alignment and master_start.tsv; alignment files are read in sorted order.

Usage:
  ./step12_identify_start_codons.sh <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv>

Inputs:
  - <step9_trimmed_dir>: Trimmed alignments from step 9.
  - <output_dir>: Destination for the start codon tables.
  - <bed_file>: Reference BED file (still a required argument).
  - <reference_species>: The reference species name prefix.
  - <step11_summary_tsv>: The step 11 *_final_validation.tsv (must contain actual_nt; the script stops otherwise).

Note: a species split over several rows (exons on different chromosomes) is read from its last row only. run_stage3.sh runs step 8c before step 9, so every species has one row; when running steps on their own, run step 8c first or remove such species.
