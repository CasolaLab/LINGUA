================================================================================
GENOMES2GENES README
================================================================================

--------------------------------------------------------------------------------
0. run_stage3.sh
--------------------------------------------------------------------------------
Description: Runs stage 3 steps 1-17 for one focal (reference) species in one working folder: input preparation
(exon list and header-cleaned reference CDS from the BED), step 1 (mafExtractor, or a copy of existing step 1
files), steps 2-8, step 8c (one collinear row per species and the missing-data mask, from LINGUA), and steps 9-17,
with step 14 reading the step 8c mask. Stops at the first failed step.

Usage:
  bash run_stage3.sh <workdir> <focal_species> <bed> <cds_fasta> <maf> <tree_newick>
  Run detached and low-priority, e.g.
  setsid nohup nice -n 19 ionice -c3 env THREADS=8 bash run_stage3.sh ... > <workdir>.log 2>&1 < /dev/null &

Inputs:
  - <workdir>: A new folder (created if needed); it must not already contain step outputs unless START_STEP is used,
    and must not contain any .bed file other than <focal_species>.bed (step 8c reads the only .bed there).
  - <focal_species>: The species prefix of the MAF sequence names (e.g. Arabidopsis_thaliana, Homo_sapiens).
  - <bed>: Tab-separated: chrom, start (0-based), end, name, score, strand, transcript ID. A header line and '#'
    lines are skipped; a trailing '+' on chromosome names is removed. The BED must be 0-based half-open: the step
    8c mask is in BED coordinates and step 14 converts it as 0-based, so the run stops before step 14 if step 9
    detects a 1-based BED.
  - <cds_fasta>: Reference CDS nucleotide sequences named by the BED transcript IDs (a '.p' suffix is accepted).
    The run stops if a BED transcript has no CDS.
  - <maf>: The reference-anchored MAF. Give '-' instead when STEP1_SRC is set
    or START_STEP > 1; PROVENANCE.txt then records "MAF not given (not used: ...)". A path, when given, must exist
    and is recorded with its size and date (not checksummed).
  - <tree_newick>: Newick tree with the MAF species names as tips (step 8).

Environment (all optional):
  - ENV_BIN: bin folder of a Python env with Biopython and pandas (default /opt/miniconda3/envs/busco_env/bin).
  - MAFEXTRACTOR_BIN: folder holding mafExtractor (default /opt/cactus/cactus/bin).
  - GNU parallel must already be on PATH (step 4).
  - THREADS (8); JOBS (parallel step 1 chunks, default THREADS; each chunk reads the whole MAF).
  - STEP3=dominant|strict (default dominant: step3_dominant_chain.py; strict: step3_filtering_species.py); JUMP (50000).
  - MAXSPAN (50000): largest genomic span of the exon chain step 8c keeps per species.
  - ONE_ROW is no longer used (step 8c always keeps one row per species); the run stops if it is set.
  - STEP1_SRC: folder of existing step 1 files (<species>.<chrom>_<start>_<end>.txt) to copy instead of running
    mafExtractor; mafExtractor is then not needed.
  - ALLOW_MISSING=1: continue when some exons have no step 1 file (default: stop). Their transcripts fail step 9.
  - START_STEP=1-17 (default 1): resume in an existing workdir (see below).
  - FLOOR_GB (50): refuse to start with less free space in the workdir.

Checks before anything runs (the run stops, without writing STATUS, if one fails):
  inputs exist; mafExtractor on PATH if step 1 runs without STEP1_SRC; GNU parallel if step 4 runs; Biopython and
  pandas; the tree can be read (warning if the focal species is not a tip); free space; no other .bed in the
  workdir; no outputs of steps >= START_STEP already in the workdir; outputs of steps < START_STEP present.

Step 1 checks (after extraction or copying): every exon must have a non-empty file containing reference lines
(problems listed in step1_problems.tsv, run stops); exons without a file are listed in step1_missing.txt and stop
the run unless ALLOW_MISSING=1. Extraction failures are listed in step1/step1_failed.tsv.

Step 14 check: step 9's report must say "# bed_convention: 0-based" (run stops at step14_bed_check otherwise).

Resuming (START_STEP=N): input preparation reruns (cheap and deterministic); steps before N are skipped; outputs of
step N and later must be moved out of the workdir first (the script never deletes anything). Step 8c counts as part
of step 8: START_STEP=8 reruns steps 8 and 8c, START_STEP=9 needs step8c/.
Outputs per step: 1 step1/; 2 step2/; 3 step3/, step3_summary.csv; 4 step4/; 6 step4/extraction_summary_fixed.tsv;
8 step8/, step8c/ (with step8c_missing_data.tsv), step8c_log.tsv; 9 step9/; 11 step8c/master_summary_final_validation.tsv;
12 step12/; 13 step13/; 14 step14/, step14_masked_deletions.tsv; 15 step15/; 16 step16_starts_stops.tsv;
17 step17_all_species.tsv. Steps 5, 7 and 10 only report (in the log).

Status and records in the workdir:
  - STATUS: "RUNNING since <date>", then "SUCCEEDED <date>" or "FAILED at <step> (exit N) <date>"
    (read by tools/run_monitor.py).
  - PROVENANCE.txt: appended on every run: inputs with md5s, options, repo commit, tool paths.
  - coords.tsv, N_EXONS, ref_raw.fna, ref.fna, <focal_species>.bed: prepared inputs.

The tree is also passed to step 8c and step 17, which type every species as living (a tip) or ancestor (a labeled
internal node), whatever the ancestor naming scheme. Step 17 gets the step 8c mask, so its coordinates cover only
the exons each species kept.

What to check after a run: STATUS; step1_missing.txt and step1_problems.tsv empty; step8c_log.tsv (species rows
rebuilt by step 8c, typed living or ancestor, and every exon gapped, with its reason); no not_in_tree types in
step8c_log.tsv or step17_all_species.tsv (a species missing from the tree); refused files in step9/step9_trimming_report.tsv;
step14_masked_deletions.tsv (deletions not reported or trimmed because of the mask); audits/audit_stage3_endtoend.py
and, with STEP3=dominant, audits/audit_optionA.py (these audit scripts are not part of this repository).

Known limitation: steps 12-17 never see insertions relative to the reference; step 14 reports deletions only (see
its README). In the Capra hircus test run (2026-11-10) the reference rows had no gaps already in the step 1 MAFs (13,071 reference rows, 0 gap characters), so insertions never reached the alignments; step 9 would also remove any column where
the reference has a gap.


--------------------------------------------------------------------------------
1. step1_mafextractor.sh
--------------------------------------------------------------------------------
Description: A wrapper for mafExtractor that reads a list of target coordinates and extracts localized MAF blocks from a master MAF file, one file per exon. Each exon is written to a temporary file and kept as <species>.<chrom>_<start>_<stop>.txt only if mafExtractor exits successfully and the output contains at least one reference line (s <species>.). Failed exons are listed on stderr and in <output_dir>/step1_failed.tsv, and the script exits with status 1; it also stops at once if mafExtractor is not on PATH or the MAF or coordinate file cannot be read. The coordinate file is read, not modified (Windows line endings are removed while reading).

Usage:
  ./step1_mafextractor.sh <species_name> <reference_species_coordinates_file> <maf_file> <output_dir>

Inputs:
  - <species_name>: The reference species prefix as used in the MAF sequence names (e.g., Arabidopsis_thaliana).
  - <reference_species_coordinates_file>: A tab-separated file with a header line, columns sequence_name (chromosome ID), start, end (BED coordinates).
  - <maf_file>: The master reference-anchored MAF file.
  - <output_dir>: Destination folder for the extracted localized MAF files.

Notes:
  - mafExtractor treats --stop as inclusive, so each file holds end - start + 1 reference bases: one base past the BED exon end. Step 9 removes it.
  - run_stage3.sh can copy existing step 1 files instead (STEP1_SRC) and then checks that every exon has a non-empty file with reference lines.


--------------------------------------------------------------------------------
2. step2_missingbp_summaries.sh
--------------------------------------------------------------------------------
Description: Parses localized MAF files to map species presence/absence and calculate missing base pairs between alignment blocks.

Usage:
  ./step2_missingbp_summaries.sh <mafextractor_dir> <output_dir>

Inputs:
  - <mafextractor_dir>: Directory containing the .txt MAF files from Step 1.
  - <output_dir>: Directory to save the generated CSV summaries.

Notes:
  - optional step to record unaligned regions of non-reference species between maf blocks.
--------------------------------------------------------------------------------
3. step3_dominant_chain.py  (default step 3 in run_stage3.sh; STEP3=dominant)
--------------------------------------------------------------------------------
Description: For every exon (one step 1 MAF file) and every non-reference species, including ancestral nodes: groups the species' blocks by chromosome and strand; within each group finds the best collinear chain (blocks in MAF order whose species coordinates increase without overlapping, at most the jump threshold apart, maximizing aligned bases); keeps the group whose chain has the most aligned bases and only that chain; and removes the species from all its other blocks in that exon. Reference blocks are never touched; blocks that lose the reference are not written. Ties go to the earliest blocks and the first group seen.

Usage:
  python3 step3_dominant_chain.py -i <mafextractor_dir> -o <output_dir> -r <reference_species> -s <summary_csv> [-j 50000] [-t 8]

Inputs:
  - -i / --mafextractor_dir: Directory with step 1 MAF files.
  - -o / --output_dir: Directory for the filtered MAF files.
  - -r / --ref: The reference species name.
  - -s / --summary: Output CSV: file_name, species, status (clean = nothing dropped, trimmed = blocks dropped), blocks_kept, blocks_dropped, bases_kept, bases_dropped.
  - -j / --jump: (Optional) Max bp between consecutive blocks of a chain (Default: 50000).
  - -t / --threads: (Optional) Number of threads (Default: 8).


--------------------------------------------------------------------------------
3 (alternative). step3_filtering_species.py  (strict step 3; run_stage3.sh uses it with STEP3=strict)
--------------------------------------------------------------------------------
Description: Removes a species from an exon's MAF file if, between consecutive blocks of that species, it changes chromosome or jumps more than the jump threshold (|start - previous end| > 50 kb by default). The species is removed from every block of that exon; the reference species is never removed, and blocks without the reference are not written. Strand is not considered. Files and summary rows are processed and written in sorted order, so reruns give identical files.

Usage:
  python3 step3_filtering_species.py -i <mafextractor_dir> -o <output_dir> -r <reference_species> -s <output_summary_csv_name> [-j <jump_threshold>] [-t <threads>]

Inputs:
  - -i / --mafextractor_dir: Directory with step 1 MAF files.
  - -o / --output_dir: Directory for the filtered MAF files.
  - -r / --ref: The reference species name to protect.
  - -s / --summary: Path to output the purge log CSV (file_name, species_purged, status).
  - -j / --jump: (Optional) Max allowed base pair jump (Default: 50000).
  - -t / --threads: (Optional) Number of threads to use (Default: 8).

Note: the default step 3 in run_stage3.sh is step3_dominant_chain.py (see its README).


--------------------------------------------------------------------------------
4. step4_reverse_complement_add_missingbp.sh
--------------------------------------------------------------------------------
Description: Converts the filtered MAF files into aligned FASTA (.aln), one per exon: concatenates each species' blocks in MAF order, pads a species that is absent from a block with dashes so all rows keep the same length, and reverse-complements every row when the reference exon is on the minus strand in the BED. Bases a species skips between blocks are unaligned sequence and are not added. Writes <output_dir>/<exon>/<exon>.aln and <output_dir>/extraction_summary.tsv (rows in sorted order). Requires GNU parallel.

Usage:
  ./step4_reverse_complement_add_missingbp.sh -m <step3_filtered_dir> -o <output_dir> -b <bed_file> [-t <threads>]

Inputs:
  - -m: Directory of filtered MAF files from step 3.
  - -o: Output directory for .aln files.
  - -b: BED file of target regions (used for strand lookup).
  - -t: (Optional) Number of parallel jobs (Default: 4).


--------------------------------------------------------------------------------
5. step5_validate_reverse_complement.py
--------------------------------------------------------------------------------
Description: Uses local sequence alignment (dynamic programming) to compare the original raw MAF sequence against the final .aln FASTA sequence to verify strand orientation.

Usage:
  python3 step5_validate_reverse_complement.py <step4_output_dir> <mafextractor_dir> <bed_file> <reference_species> <threads>

Inputs:
  - <step4_output_dir>: Directory containing extracted .aln FASTAs.
  - <mafextractor_dir>: Directory containing original unfiltered MAF files.
  - <bed_file>: Reference BED file with strand targets.
  - <reference_species>: Name of the reference species.
  - <threads>: Number of threads to use.


--------------------------------------------------------------------------------
6. step6_update_summary_files.py
--------------------------------------------------------------------------------
Description: A post-processing audit script that synchronizes the FASTA extraction summary TSV with the original BED file.

Usage:
  python3 step6_update_summary_files.py <step4_summary_tsv> <bed_file> <reference_species>

Inputs:
  - <step4_summary_tsv>: The extraction_summary.tsv generated by script #4.
  - <bed_file>: The master reference BED file.
  - <reference_species>: The reference species name prefix.


--------------------------------------------------------------------------------
7. step7_first_stage_validation.py
--------------------------------------------------------------------------------
Description: A data-integrity checker that counts valid (non-dash) nucleotides across the Raw MAF, Filtered MAF, and fasta extracted files.

Usage:
  python3 step7_first_stage_validation.py <raw_mafextractor_dir> <step3_filtered_dir> <step4_output_dir> <reference_species>

Inputs:
  - <raw_mafextractor_dir>: Directory containing original unfiltered MAF files.
  - <step3_filtered_dir>: Directory of filtered MAF files.
  - <step4_output_dir>: Directory containing extracted .aln FASTAs.
  - <reference_species>: The reference species name prefix.


--------------------------------------------------------------------------------
8. step8_stitching_exons.py
--------------------------------------------------------------------------------
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


--------------------------------------------------------------------------------
8c. step8c_one_collinear_row.py  (from LINGUA stage 3; runs after step 8, before step 9)
--------------------------------------------------------------------------------
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

Notes:
  - Change from the LINGUA copy: BED header and '#' lines are skipped by content (non-numeric start), as in step 8,
    instead of always skipping the first line, so a BED without a header keeps its first exon. The log's type
    column comes from the tree, not from the name pattern N<number>, which typed every mammal ancestor (e.g.
    fullTreeAnc235) as living: 777 of the 1,561 rows rebuilt in the Capra hircus test run (2026-11-10).


--------------------------------------------------------------------------------
9. step9_trimming_duplicate_bases.py  (LINGUA 69ca427, per-exon padding removal)
--------------------------------------------------------------------------------
Description: Removes from every stitched alignment the one extra reference base that step 1 extracted per exon, at its known position from the BED (last base of each exon block on plus-strand genes, first base on minus-strand genes; no extra base on an exon that ends at the end of its scaffold), and removes columns where the reference has a gap, so every species is projected onto the reference CDS. Because the extra base is located by position rather than by matching sequence, every species is cut at the same correct column even when the extra base equals its neighbor. A file is written only if the trimmed reference row equals its CDS in the FNA (N or IUPAC ambiguity codes may match a base they stand for); every other file is refused and listed with its reason (e.g. REFERENCE_MISMATCH_AFTER_TRIM, EXTRACTION_SHIFTED, NOT_IN_FNA) in <output_dir>/step9_trimming_report.tsv. The BED convention (0-based or 1-based) is detected from the data.

Usage:
  python3 step9_trimming_duplicate_bases.py <step8_stitched_dir> <output_dir> <reference_nucleotide_fna_file> <bed_file> <threads>

Inputs:
  - <step8_stitched_dir>: Stitched alignments from step 8c (run_stage3.sh passes step8c/; plain step 8 output also works). Searched recursively, so keep other .aln files out of it.
  - <output_dir>: Destination for the trimmed .aln files; must be new or empty.
  - <reference_nucleotide_fna_file>: FASTA of the reference CDS sequences, named by the BED transcript IDs.
  - <bed_file>: Reference BED file containing exon coordinates and transcript IDs.
  - <threads>: Number of threads to use.

Notes: before running step 9 on its own, run clean_nucleotide_fasta_headers.py on the reference nucleotide FASTA (run_stage3.sh does this).


--------------------------------------------------------------------------------
10. step10_second_stage_validation.py
--------------------------------------------------------------------------------
Description: An audit script that checks the final trimmed/stitched sequences against the reference FASTA to guarantee that length and base content match perfectly.

Usage:
  python3 step10_second_stage_validation.py <step9_trimmed_dir> <reference_nucleotide_fna_file> <reference_species>

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <reference_nucleotide_fna_file>: Fasta file of true reference nucleotide sequences.
  - <reference_species>: The reference species name prefix.


--------------------------------------------------------------------------------
11. step11_update_stitched_exons_summary.py
--------------------------------------------------------------------------------
Description: Adds to the stitched-alignment master_summary.tsv (step8c/master_summary.tsv in run_stage3.sh) the CDS length from the FNA (expected_fna), the reference length in the step 9 alignment with gaps removed (actual_nt) and a validation status (PASS, PASS_IN_FRAME_SHORT, FAIL_FRAME_SHIFT, FAIL_MISSING_FNA_TRUTH, FAIL_MISSING_ALN_FILE, ...), and writes <step8_summary>_final_validation.tsv. Steps 12-15 read this file; steps 12-14 take reference_length from its actual_nt column.

Usage:
  python3 step11_update_stitched_exons_summary.py <step8_output_summary_tsv> <step9_trimmed_dir> <bed_file> <reference_nucleotide_fna_file> <reference_species>

Inputs:
  - <step8_output_summary_tsv>: master_summary.tsv from step 8c (run_stage3.sh) or step 8; the output is written next to it.
  - <step9_trimmed_dir>: Trimmed alignments from step 9.
  - <bed_file>: Reference BED file containing exon coordinates and transcript IDs.
  - <reference_nucleotide_fna_file>: The same (header-cleaned) CDS FASTA given to step 9; with another file every row fails with FAIL_MISSING_FNA_TRUTH and actual_nt = 0.
  - <reference_species>: The reference species name prefix.


--------------------------------------------------------------------------------
12. step12_identify_start_codons.sh
--------------------------------------------------------------------------------
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


--------------------------------------------------------------------------------
13. step13_identify_stop_codons.sh
--------------------------------------------------------------------------------
Description: For every non-reference species in each step 9 alignment, removes that species' gaps and reads its bases in groups of three from its own first aligned base, reporting every codon that is TAA, TAG, TGA. Positions are reported in reference coordinates: the reference position of the alignment column holding the codon's first and third base. Transcript ID, chromosome and strand come from the step 11 summary (matched by file name); reference_length is the step 11 actual_nt column (found by name); cds_length is that species' own length with gaps removed. Writes one <file>_stops.tsv per alignment and master_stop.tsv; alignment files are read in sorted order.

Usage:
  ./step13_identify_stop_codons.sh <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv>

Inputs:
  - <step9_trimmed_dir>: Trimmed alignments from step 9.
  - <output_dir>: Destination for the stop codon tables.
  - <bed_file>: Reference BED file (still a required argument).
  - <reference_species>: The reference species name prefix.
  - <step11_summary_tsv>: The step 11 *_final_validation.tsv (must contain actual_nt; the script stops otherwise).

Note: a species split over several rows (exons on different chromosomes) is read from its last row only. run_stage3.sh runs step 8c before step 9, so every species has one row; when running steps on their own, run step 8c first or remove such species.


--------------------------------------------------------------------------------
14. step14_identify_indels.sh  (LINGUA masked version)
--------------------------------------------------------------------------------
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


--------------------------------------------------------------------------------
15. step15_identify_first_and_last_three_bases.sh
--------------------------------------------------------------------------------
Description: Extracts the exact first and last 3 base pairs for a species from the multiple sequence alignment.

Usage:
  ./step15_identify_first_and_last_three_bases.sh <step9_trimmed_dir> <output_dir> <reference_species> <step11_summary_tsv>

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <output_dir>: Destination for the summary file.
  - <reference_species>: The reference species name(or other species name of interest).
  - <step11_summary_tsv>: Summary TSV file from stitched and trimmed exons in step 11.

Note: alignment files are read in sorted order, so reruns write rows in the same order.


--------------------------------------------------------------------------------
16. step16_extract_start_and_stops.py
--------------------------------------------------------------------------------
Description: A fast directory-wide parser that pulls the first and last 3 bases of sequences, standardizing gap removal and validating if they match known translation start/stop motifs.

Usage:
  python3 step16_extract_start_and_stops.py <step9_trimmed_dir> <output.tsv>

Inputs:
  - <step9_trimmed_dir>: Final stitched and trimmed alignments from step 9.
  - <output.tsv>: Output summary TSV file name.


--------------------------------------------------------------------------------
17. step17_combine_all_species_info.py
--------------------------------------------------------------------------------
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
    (step 8c not run) keeps a comma-joined list of its row contigs, as before 2026-11-10.
  - Before 2026-11-10 every block of every exon was used. In the Capra hircus test run (2026-11-10, 632 alignments,
    15,408 rows) the fix changed 1,965 rows: 432 that spanned several contigs (all now on one), 1,123 other rows of
    species step 8c rebuilt (span narrower), and 410 rows of species it did not rebuild, all narrowed by skipping the
    extra base past an exon end (154 by 1 bp; 256 by more, where that base was aligned a few bases away).
  - The transcript is matched to the file name by the longest BED transcript ID that prefixes it (before, the first
    match in BED order).
