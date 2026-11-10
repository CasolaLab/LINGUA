================================================================================
SCRIPT: step3_dominant_chain.py  (default step 3 in run_stage3.sh; STEP3=dominant)
================================================================================
Description: For every exon (one step 1 MAF file) and every non-reference species, including ancestral nodes: groups the species' blocks by chromosome and strand; within each group finds the best collinear chain (blocks in MAF order whose species coordinates increase without overlapping, at most the jump threshold apart, maximizing aligned bases); keeps the group whose chain has the most aligned bases and only that chain; and removes the species from all its other blocks in that exon. Unlike step3_filtering_species.py, the species stays in the exon with its main chain instead of being removed from the whole exon. Reference blocks are never touched; blocks that lose the reference are not written. Ties go to the earliest blocks and the first group seen.

Usage:
  python3 step3_dominant_chain.py -i <mafextractor_dir> -o <output_dir> -r <reference_species> -s <summary_csv> [-j 50000] [-t 8]

Inputs:
  - -i / --mafextractor_dir: Directory with step 1 MAF files.
  - -o / --output_dir: Directory for the filtered MAF files.
  - -r / --ref: The reference species name.
  - -s / --summary: Output CSV: file_name, species, status (clean = nothing dropped, trimmed = blocks dropped), blocks_kept, blocks_dropped, bases_kept, bases_dropped.
  - -j / --jump: (Optional) Max bp between consecutive blocks of a chain (Default: 50000).
  - -t / --threads: (Optional) Number of threads (Default: 8).

