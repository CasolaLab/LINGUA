#!/usr/bin/env bash
# run_stage3.sh
#
# LINGUA stage 3 (alignment extraction), steps 1-17, for one focal species of a Progressive Cactus alignment.
# Corrected pipeline (not the Brassicaceae paper reproduction):
#   step 3   dominant collinear chain per species and exon (default), keeping the orthologous blocks and dropping
#            blocks on other chromosomes/strands, more than JUMP bp apart, out of order or overlapping;
#            STEP3=strict uses the original filter, which removes the species from the whole exon instead
#   step 4   no genome lookup: bases a species skips between MAF blocks are not added
#   step 8   Python port, run without --one-row-per-species (every species.contig row is kept for step 8c)
#   step 8c  step8c_one_collinear_row.py: one row per species and ancestor, from its dominant collinear chain
#            of exons (one contig and strand, collinear with the reference exon order, within MAXSPAN bp); every
#            other exon becomes gaps, and every reference range a species has no sequence in (exons it is absent
#            from, exons gapped here) goes to the missing-data mask step8c/step8c_missing_data.tsv
#   step 9   per-exon padding removal, which fixes the exon-junction errors; reads step8c/
#   step 14  masked version: deletions inside the step 8c mask are missing data and not reported; deletions
#            overlapping it are trimmed to their unmasked parts, each tested as a frameshift; every masked or
#            trimmed deletion is listed in step14_masked_deletions.tsv (outside step14/, which stage 4 reads)
#   step 14b step14b_deletion_evidence.py: for every reported deletion, where its gap came from (aligner gap, no
#            alignment, piece dropped by step 3) and whether the species has its own sequence in place, in
#            step14_deletion_evidence.tsv; step 14's own files are unchanged
#
# Usage:
#   bash run_stage3.sh <workdir> <focal_species> <bed> <cds_fasta> <maf> <tree_newick>
#   <focal_species> must match the species prefix of the MAF sequence names (e.g. Homo_sapiens).
#   <bed>: tab-separated, columns chrom, start (0-based), end, name, score, strand, transcript ID (a header line
#          and '#' comment lines are skipped).
#   <cds_fasta>: CDS nucleotide sequences named by the BED transcript IDs.
#   <maf>: read only by step 1 (mafExtractor); give '-' instead when STEP1_SRC is set or START_STEP > 1.
# Environment:
#   ENV_BIN      python env bin with Biopython and pandas (default busco_env); put GNU parallel on PATH too
#   MAFEXTRACTOR_BIN  folder holding mafExtractor (default /opt/cactus/cactus/bin)
#   THREADS (8), JOBS (parallel step 1 chunks, default THREADS), STEP3 (dominant|strict), JUMP (50000),
#   MAXSPAN (50000, step 8c chain span), FLOOR_GB (50, refuse to start below this free space)
#   STEP1_SRC    folder of existing step 1 files to copy instead of running mafExtractor
#   ALLOW_MISSING=1  continue when some exons have no step 1 file (default: stop); they are listed in
#                step1_missing.txt and their transcripts will not pass step 9
#   START_STEP   1-17 (default 1): resume in an existing workdir. Input preparation always reruns (it is cheap and
#                deterministic); steps before START_STEP are skipped and their outputs must exist; outputs of
#                START_STEP and later must NOT exist (move them away first; nothing is deleted by this script).
#                Step 8c belongs to step 8: START_STEP=8 reruns both, START_STEP=9 needs step8c/.
# Status: <workdir>/STATUS says RUNNING, then "SUCCEEDED <date>" or "FAILED at <step> (exit N) <date>"
#   (tools/run_monitor.py reads it).
# Run detached and low-priority, e.g.
#   setsid nohup nice -n 19 ionice -c3 bash run_stage3.sh ... > <workdir>.log 2>&1 < /dev/null &
set -euo pipefail

WD="${1:?usage: run_stage3.sh <workdir> <focal_species> <bed> <cds_fasta> <maf> <tree_newick>}"
SP="${2:?focal species}"
# readlink -m, not -f: -f fails (and set -e exits silently) when a parent folder is missing, before the input check
BED_SRC=$(readlink -m "${3:?bed}")
CDS_SRC=$(readlink -m "${4:?cds fasta}")
MAF="${5:?maf (or - when STEP1_SRC is set or START_STEP > 1)}"
[ "$MAF" = - ] || MAF=$(readlink -m "$MAF")
export TREE_FILE=$(readlink -m "${6:?tree}")
REPO=$(cd "$(dirname "$0")" && pwd)
export PATH=${ENV_BIN:-/opt/miniconda3/envs/busco_env/bin}:${MAFEXTRACTOR_BIN:-/opt/cactus/cactus/bin}:$PATH
S=$REPO/scripts
T=${THREADS:-8}
JOBS=${JOBS:-$T}
STEP3=${STEP3:-dominant}
JUMP=${JUMP:-50000}
MAXSPAN=${MAXSPAN:-50000}
FLOOR_GB=${FLOOR_GB:-50}
START_STEP=${START_STEP:-1}
B=$SP.bed

# ---------------------------------------------------------------- checks before anything runs
die() { echo "ERROR: $*" >&2; exit 2; }
for f in "$BED_SRC" "$CDS_SRC" "$TREE_FILE"; do [ -e "$f" ] || die "missing input $f"; done
[ -z "${STEP1_SRC:-}" ] || [ -d "$STEP1_SRC" ] || die "missing STEP1_SRC $STEP1_SRC"
# the MAF is read only when step 1 runs mafExtractor; otherwise '-' may be given instead of a path
if [ "$MAF" = - ]; then
    [ "${START_STEP:-1}" != 1 ] || [ -n "${STEP1_SRC:-}" ] || die "a MAF is needed when step 1 runs mafExtractor (no STEP1_SRC)"
else
    [ -e "$MAF" ] || die "missing input $MAF"
fi
case $STEP3 in dominant|strict) ;; *) die "STEP3 must be dominant or strict";; esac
[ -z "${ONE_ROW:-}" ] || die "ONE_ROW is no longer used: step 8c always keeps one collinear row per species (unset it)"
[[ "$MAXSPAN" =~ ^[0-9]+$ ]] || die "MAXSPAN must be a whole number of bp"
[[ "$START_STEP" =~ ^([1-9]|1[0-7])$ ]] || die "START_STEP must be 1-17"
if [ "$START_STEP" -eq 1 ] && [ -z "${STEP1_SRC:-}" ]; then
    command -v mafExtractor >/dev/null || die "mafExtractor not on PATH (set MAFEXTRACTOR_BIN), needed for step 1"
fi
if [ "$START_STEP" -le 4 ]; then
    command -v parallel >/dev/null || die "GNU parallel not on PATH, needed for step 4"
fi
python3 -c "import Bio, pandas" 2>/dev/null || die "$(command -v python3) lacks Biopython or pandas (set ENV_BIN)"
python3 - "$TREE_FILE" "$SP" <<'EOF' || die "cannot read the tree $TREE_FILE"
import sys
from Bio import Phylo
tips = {t.name for t in Phylo.read(sys.argv[1], "newick").get_terminals() if t.name}
if sys.argv[2] not in tips:
    print(f"WARNING: {sys.argv[2]} is not a tip of the tree; step 8's distant_species will be a species pair", file=sys.stderr)
EOF

mkdir -p "$WD" && cd "$WD"
free=$(df -BG --output=avail . | tail -1 | tr -dc 0-9)
[ "$free" -ge "$FLOOR_GB" ] || die "only $free GB free (< $FLOOR_GB GB)"
# step 8c takes the reference species from the only .bed in the workdir, so no other .bed may sit next to $B
others=$(find . -maxdepth 1 -name '*.bed' ! -name "$B" | head -3)
[ -z "$others" ] || die "step 8c needs $B to be the only .bed in $WD; move these away: $(echo $others)"

# outputs of each step (used to refuse overwriting and to check what a resumed run needs); step 8 includes step 8c
declare -A OUT=(
  [1]="step1" [2]="step2" [3]="step3 step3_summary.csv" [4]="step4" [5]="" [6]="step4/extraction_summary_fixed.tsv"
  [7]="" [8]="step8 step8c step8c_log.tsv" [9]="step9" [10]="" [11]="step8c/master_summary_final_validation.tsv"
  [12]="step12" [13]="step13" [14]="step14 step14_masked_deletions.tsv step14_deletion_evidence.tsv" [15]="step15"
  [16]="step16_starts_stops.tsv" [17]="step17_all_species.tsv")
for k in $(seq "$START_STEP" 17); do
    for p in ${OUT[$k]}; do
        [ -e "$p" ] && die "$WD/$p already exists (output of step $k >= START_STEP=$START_STEP); move it away first"
    done
done
for k in $(seq 1 $((START_STEP - 1))); do
    for p in ${OUT[$k]}; do
        [ -e "$p" ] || die "START_STEP=$START_STEP needs the output of step $k, but $WD/$p does not exist"
    done
done

# ---------------------------------------------------------------- status file and step runner
CURRENT=setup
echo "RUNNING since $(date -Is)" > STATUS
finish() {
    rc=$?
    if [ "$rc" -eq 0 ] && [ "$CURRENT" = done ]; then echo "SUCCEEDED $(date -Is)" > STATUS
    else echo "FAILED at $CURRENT (exit $rc) $(date -Is)" > STATUS; fi
}
trap finish EXIT
run() {   # run <step number> <name> <command...>: skipped when <step number> < START_STEP
    local n="$1" name="$2"; shift 2
    if [ "$n" -lt "$START_STEP" ]; then echo "== $name skipped (START_STEP=$START_STEP)"; return 0; fi
    CURRENT="$name"
    echo "== $name $(date -Is)"
    "$@"
}

{
    echo "run_stage3.sh $(date -Is)"
    echo "repo: $REPO ($(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo 'no git'))"
    echo "focal species: $SP; STEP3=$STEP3 JUMP=$JUMP MAXSPAN=$MAXSPAN THREADS=$T JOBS=$JOBS START_STEP=$START_STEP ALLOW_MISSING=${ALLOW_MISSING:-0}"
    for f in "$BED_SRC" "$CDS_SRC" "$TREE_FILE"; do echo "input: $f md5 $(md5sum < "$f" | cut -c1-32)"; done
    if [ "$MAF" = - ] && [ -n "${STEP1_SRC:-}" ]; then
        echo "input: MAF not given (not used: step 1 copied from $STEP1_SRC)"
    elif [ "$MAF" = - ]; then
        echo "input: MAF not given (not used: step 1 skipped, START_STEP=$START_STEP)"
    else
        echo "input: $MAF ($(stat -c %s "$MAF") bytes, $(stat -c %y "$MAF"))"
    fi
    [ -n "${STEP1_SRC:-}" ] && echo "step 1 copied from: $STEP1_SRC"
    echo "mafExtractor: $(command -v mafExtractor || echo 'not found'); parallel: $(command -v parallel || echo 'not found'); python3: $(command -v python3)"
} >> PROVENANCE.txt

# ---------------------------------------------------------------- input preparation (always)
CURRENT=prepare
tr -d '\r' < "$BED_SRC" | awk 'BEGIN{FS=OFS="\t"} !/^#/ && $2 ~ /^[0-9]+$/ {sub(/\+$/,"",$1)} {print}' > "$B"
echo "== prepare $(date -Is)"
python3 "$S/prepare_species_inputs.py" "$B" "$CDS_SRC" coords.tsv N_EXONS ref_raw.fna
CURRENT=clean_headers
echo "== clean_headers $(date -Is)"
python3 "$S/clean_nucleotide_fasta_headers.py" -f ref_raw.fna -b "$B" -o ref.fna

# ---------------------------------------------------------------- step 1
check_step1() {   # every exon has a non-empty step 1 file with reference lines
    local want n missing problems
    want=$(cat N_EXONS)
    : > step1_missing.txt; : > step1_problems.tsv
    tail -n +2 coords.tsv | while IFS=$'\t' read -r c s e; do
        f="step1/$SP.${c}_${s}_${e}.txt"
        if [ ! -e "$f" ]; then echo "$SP.${c}_${s}_${e}.txt" >> step1_missing.txt
        elif [ ! -s "$f" ]; then echo -e "$f\tempty file" >> step1_problems.tsv
        elif ! grep -q "^s $SP\." "$f"; then echo -e "$f\tno reference (s $SP.) lines" >> step1_problems.tsv
        fi
    done
    n=$(find step1 -maxdepth 1 -name '*.txt' | wc -l)
    missing=$(wc -l < step1_missing.txt); problems=$(wc -l < step1_problems.tsv)
    echo "step 1: $n files for $want exons; $missing missing (step1_missing.txt); $problems empty or without reference lines (step1_problems.tsv)"
    [ "$problems" -eq 0 ] || { echo "ERROR: step 1 has empty or invalid files; see step1_problems.tsv" >&2; return 1; }
    if [ "$missing" -gt 0 ] && [ "${ALLOW_MISSING:-0}" != 1 ]; then
        echo "ERROR: $missing exon(s) have no step 1 file (step1_missing.txt); set ALLOW_MISSING=1 to continue without them" >&2
        return 1
    fi
    [ "$n" -gt 0 ] || { echo "ERROR: no step 1 files" >&2; return 1; }
}

step1_copy() {
    mkdir -p step1
    tail -n +2 coords.tsv | while IFS=$'\t' read -r c s e; do
        f="$SP.${c}_${s}_${e}.txt"
        if [ -f "$STEP1_SRC/$f" ]; then cp -p "$STEP1_SRC/$f" step1/; fi
    done
}

step1_extract() {
    mkdir -p step1 chunks
    tail -n +2 coords.tsv | split -n l/"$JOBS" -d -a 3 - chunks/c_
    local pids=() fail=0
    for c in chunks/c_???; do
        { printf 'seq\tstart\tstop\n'; cat "$c"; } > "$c.tsv"
        bash "$S/step1_mafextractor.sh" "$SP" "$c.tsv" "$MAF" step1 > "$c.log" 2>&1 &
        pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p" || fail=1; done
    [ "$fail" -eq 0 ] || { echo "ERROR: step 1 failed for some exons (step1/step1_failed.tsv, chunks/*.log)" >&2; return 1; }
}

if [ -n "${STEP1_SRC:-}" ]; then run 1 step1 step1_copy; else run 1 step1 step1_extract; fi
CURRENT=step1_check
check_step1

# ---------------------------------------------------------------- steps 2-17
run 2 step2 bash "$S/step2_missingbp_summaries.sh" step1 step2
if [ "$STEP3" = dominant ]; then
    run 3 step3 python3 "$S/step3_dominant_chain.py" -i step1 -o step3 -r "$SP" -s step3_summary.csv -j "$JUMP" -t "$T"
else
    run 3 step3 python3 "$S/step3_filtering_species.py" -i step1 -o step3 -r "$SP" -s step3_summary.csv -t "$T"
fi
run 4 step4 bash "$S/step4_reverse_complement_add_missingbp.sh" -m step3 -o step4 -b "$B" -t "$T"
run 5 step5 python3 "$S/step5_validate_reverse_complement.py" step4 step1 "$B" "$SP" "$T"
run 6 step6 python3 "$S/step6_update_summary_files.py" step4/extraction_summary.tsv "$B" "$SP"
run 7 step7 python3 "$S/step7_first_stage_validation.py" step1 step3 step4 "$SP"
run 8 step8 python3 "$S/step8_stitching_exons.py" step4 "$B" step8 "$SP"
run 8 step8c python3 "$S/step8c_one_collinear_row.py" . step8c step8c_log.tsv --maxspan "$MAXSPAN" --tree "$TREE_FILE"
run 9 step9 python3 "$S/step9_trimming_duplicate_bases.py" step8c step9 ref.fna "$B" "$T"
run 10 step10 python3 "$S/step10_second_stage_validation.py" step9 ref.fna "$SP"
run 11 step11 python3 "$S/step11_update_stitched_exons_summary.py" step8c/master_summary.tsv step9 "$B" ref.fna "$SP"
V=step8c/master_summary_final_validation.tsv
run 12 step12 bash "$S/step12_identify_start_codons.sh" step9 step12 "$B" "$SP" "$V"
run 13 step13 bash "$S/step13_identify_stop_codons.sh" step9 step13 "$B" "$SP" "$V"
check_bed_convention() {   # the step 8c mask holds 0-based half-open BED coordinates; step 14 converts them as such
    local conv
    conv=$(grep -m1 '^# bed_convention:' step9/step9_trimming_report.tsv | awk '{print $3}')
    [ "$conv" = 0-based ] || { echo "ERROR: step 9 found a '${conv:-unknown}' BED; the step 8c mask needs 0-based coordinates" >&2; return 1; }
}
run 14 step14_bed_check check_bed_convention
run 14 step14 bash "$S/step14_identify_indels.sh" step9 step14 "$B" "$SP" "$V" step8c/step8c_missing_data.tsv step14_masked_deletions.tsv
run 14 step14b python3 "$S/step14b_deletion_evidence.py" step14/master_indels.tsv step1 step3 step9 "$B" "$SP" step14_deletion_evidence.tsv --tree "$TREE_FILE"
run 15 step15 bash "$S/step15_identify_first_and_last_three_bases.sh" step9 step15 "$SP" "$V"
run 16 step16 python3 "$S/step16_extract_start_and_stops.py" step9 step16_starts_stops.tsv
run 17 step17 python3 "$S/step17_combine_all_species_info.py" step9 step3 "$B" "$SP" step17_all_species.tsv --mask step8c/step8c_missing_data.tsv --tree "$TREE_FILE"
CURRENT=done
echo "== done $(date -Is)"
