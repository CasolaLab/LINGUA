# 08 Troubleshooting

Errors seen in real runs, what they mean and what to do. Every LINGUA driver stops at the first problem and says why,
so read the last lines of the log before anything else.

## Any stage

| Message or symptom | Cause | Fix |
|---|---|---|
| `declare: -A: invalid option` | the script ran under macOS's bash 3.2 | activate the tutorial environment and run scripts as `bash script.sh`, or install bash 5 |
| `readlink: illegal option -- m`, or `df` errors | macOS tools instead of GNU coreutils | activate the tutorial environment (it includes coreutils) |
| a stage refuses to start because its output exists | every stage protects its output | move the old output away, or give a new output folder |
| the disk fills during a long run, or a remote session closes right after login | temporary files on a small root filesystem (`/tmp`) | set `TMPDIR` to a folder on a data disk before the run, and check `df -h` first |

## Stage 0

| Message or symptom | Cause | Fix |
|---|---|---|
| a species is missing from `master_summary_report.tsv` | its GFF or protein file was not found under an accepted extension | check `pipeline.log` for `[WARNING] Missing input files`, and the file names |
| `N_Cleaned_Proteins` differs from `N_Isoform_Step_Input` | the isoform step read a different file than the cleaning step wrote | do not use that species until the cause is found (see the stage README) |

## Stage 1

| Message or symptom | Cause | Fix |
|---|---|---|
| warnings about protein files not found, and few unclustered candidates | `--prt_dir` points at the stage 0 folder instead of OrthoFinder's `protein/` copy | point `--prt_dir` at the `protein/` folder made by `run_orthofinder.sh` |
| exit status 1 and violations in `clsg_verification_report.tsv` | candidates that disagree with the HOG table | do not use the candidates; check that the HOG file is the one of the ingroup node |

## Stage 2

| Message or symptom | Cause | Fix |
|---|---|---|
| exit status 2, `FAILED at step N`, genes marked `NOT_RUN` | a search did not finish for some genes (killed, out of time, out of memory) | resume or rerun that step; unfinished searches are never counted as "no hit" |
| a warning that an analysis species is inside a database | the database holds the species under study | rebuild it with `make-db --exclude-inputs`; otherwise every candidate finds itself |
| a species with almost every candidate removed | usually the same self-hit problem, or wrong species tags in the database | check `<db>.species.tsv` written by `make-db` |

## Stage 3

| Message or symptom | Cause | Fix |
|---|---|---|
| `GNU parallel not on PATH` | step 4 needs GNU parallel | activate the tutorial environment |
| `lacks Biopython or pandas (set ENV_BIN)` | the driver uses the Python in `ENV_BIN`, whose default is a path on the lab workstation | `ENV_BIN="$(dirname "$(command -v python3)")"` |
| the run refuses to start for lack of space | free space below `FLOOR_GB` (50 GB by default) | free space, or lower `FLOOR_GB` for a small run |
| `FAILED at step14_bed_check` | the BED is 1-based | write 0-based starts (BED convention) |
| refused files in `step9/step9_trimming_report.tsv` | the reference row of the alignment does not equal the CDS | check that the CDS FASTA and the BED describe the same transcript |
| a run picks up the wrong BED | the work folder holds another `.bed` file | use a new, empty work folder |

## Stage 4

| Message or symptom | Cause | Fix |
|---|---|---|
| `no query FASTA for <prefix> in <folder>` at step 4 | the run covers all species in `species.tsv`, and one has no protein file | restrict the run with `ONLY="<prefix> <prefix>"` or add the missing species |
| loci that should be there are missing from the summaries | `REMOVE_LIST` still points at the published removal list | set `REMOVE_LIST=` (empty) after a LINGUA stage 2 |

## Summary stage

| Message or symptom | Cause | Fix |
|---|---|---|
| `input not found: .../<path>/<path>/annotation.gff3` | GFF paths in the GFF table are read relative to the table | give paths relative to the table's folder, or absolute paths |
| a stop on gene ID forms or on an ID collision | the same gene appears under two forms (`.p` suffix, isoform number) or two species share an ID | make the IDs agree across the inputs; the message names the IDs |
| the driver refuses to start because `OUTDIR` contains a space | BLAST+ at step 4 cannot use such paths (checked only when step 4 is part of the run) | use an output path without spaces |
| `FAILED at step N` | see the end of `OUTDIR/stepN.log`, which the driver also prints | |
