================================================================================
SCRIPT: run_stage3.sh
================================================================================
Description: Runs stage 3 steps 1-17 for one focal (reference) species in one working folder. Stops at the first failed step.

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
12 step12/; 13 step13/; 14 step14/, step14_masked_deletions.tsv, step14_deletion_evidence.tsv (step 14b); 15 step15/; 16 step16_starts_stops.tsv;
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
step14_masked_deletions.tsv (deletions not reported or trimmed because of the mask); step14_deletion_evidence.tsv
(step 14b: where each reported deletion's gap came from and whether the species has sequence in place); audits/audit_stage3_endtoend.py
and, with STEP3=dominant, audits/audit_optionA.py (these audit scripts are not part of this repository).

