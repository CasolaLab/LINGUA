#!/usr/bin/env bash
# run_stage5.sh - LINGUA stage 5: births, deaths and their rates on every branch.
#
# Usage: bash 5_birth_death/run_stage5.sh stage5.conf [first_step] [last_step]
#   stage5.conf is a bash file of variables (see config/stage5.conf.example):
#     TREE_CONFIG     stage 4 tree_config.json
#     COORDS_DIR      stage 3 step 17 tables (*_all_species_summary.tsv)
#     NODES_DIR       stage 4 final classification (*_nodes_summary_filtered.tsv)
#     TBLASTN_DIR     stage 4 tBLASTn summaries (*_cdngs_tblastn_*_summary.tsv)
#     HOGS            OrthoFinder HOG table (N0.tsv, or the node holding all species)
#     PROTEINS_DIR    candidate protein FASTAs, one per species
#     BRANCH_LENGTHS  TSV branch,length (or node,age)
#     GFF_TABLE       TSV species,gff[,fai], one GFF3 per focal species (step 3, required)
#     NO_GFF_SPECIES  optional comma-separated species with no annotation; their loci are counted without step 3
#     OUTDIR          output folder (must not exist when starting at step 1)
#     THREADS         BLASTp threads (default 4); PYTHON, BLASTP, MAKEBLASTDB optional
#     HOG_ID_MATCH    normalized (default) or exact (reproduces the Brassicaceae paper)
#     EXCLUDE_OVERLAP same_strand_cds (default: only loci whose CDS overlaps a same-strand gene's CDS are left out;
#                     antisense and intronic loci are counted) or any (every overlap left out, the paper)
#     HOG_SPAN        split (default: an orthogroup counts once per lineage) or merge (once per orthogroup, the paper)
#     PRESENCE        tblastn (default: the stage 4 tBLASTn summaries in TBLASTN_DIR) or orf (step 1a: presence in every
#                     extant genome by the stage 4 ORF rule, for stage 4 runs with orf_test = length; needs ALIGNMENTS_DIR,
#                     the folder of stage 3 alignment folders, and STAGE4_OPTIONS, the stage 4 method options of the run)
#     UNKNOWN_PRESENCE  optional TSV species,geneID,unknown_in: genomes with no usable call, which never add a death in
#                     step 5; with PRESENCE=orf it defaults to the table step 1a writes (set it to none to leave it out)
# Steps: 1 coordinates (1a presence with PRESENCE=orf), 2 loci, 3 overlap with annotated genes, 4 unclustered duplicates
# (BLASTp), 5 rates.
set -euo pipefail

[ $# -ge 1 ] || { sed -n '2,27p' "$0"; exit 1; }
CONF=$1; FIRST=${2:-1}; LAST=${3:-5}
[ -f "$CONF" ] || { echo "ERROR: config not found: $CONF" >&2; exit 1; }
case "$FIRST$LAST" in
    [1-5][1-5]) [ "$FIRST" -le "$LAST" ] || { echo "ERROR: first step $FIRST is after last step $LAST" >&2; exit 1; } ;;
    *) echo "ERROR: steps must be numbers 1-5 (got '$FIRST' '$LAST')" >&2; exit 1 ;;
esac
# shellcheck disable=SC1090
source "$CONF"
HERE=$(cd "$(dirname "$0")" && pwd)
S=$HERE/scripts
PYTHON=${PYTHON:-python3}; BLASTP=${BLASTP:-blastp}; MAKEBLASTDB=${MAKEBLASTDB:-makeblastdb}
THREADS=${THREADS:-4}; HOG_ID_MATCH=${HOG_ID_MATCH:-normalized}; GFF_TABLE=${GFF_TABLE:-}; NO_GFF_SPECIES=${NO_GFF_SPECIES:-}; HOG_SPAN=${HOG_SPAN:-split}; EXCLUDE_OVERLAP=${EXCLUDE_OVERLAP:-same_strand_cds}
PRESENCE=${PRESENCE:-tblastn}; UNKNOWN_PRESENCE=${UNKNOWN_PRESENCE:-}
MAKEBLASTDB=$(command -v "$MAKEBLASTDB" || echo "$MAKEBLASTDB"); BLASTP=$(command -v "$BLASTP" || echo "$BLASTP")
for v in TREE_CONFIG COORDS_DIR NODES_DIR TBLASTN_DIR HOGS PROTEINS_DIR BRANCH_LENGTHS GFF_TABLE OUTDIR; do
    [ -n "${!v:-}" ] || { echo "ERROR: $v is not set in $CONF" >&2; exit 1; }
done
case "$PRESENCE" in
    tblastn) ;;
    orf) for v in ALIGNMENTS_DIR STAGE4_OPTIONS; do
             [ -n "${!v:-}" ] || { echo "ERROR: PRESENCE=orf needs $v in $CONF" >&2; exit 1; }
         done
         [ -n "$UNKNOWN_PRESENCE" ] || UNKNOWN_PRESENCE=$OUTDIR/presence/unknown_presence.tsv ;;
    *) echo "ERROR: PRESENCE must be tblastn or orf (got '$PRESENCE')" >&2; exit 1 ;;
esac
[ "$UNKNOWN_PRESENCE" != none ] || UNKNOWN_PRESENCE=
if [ "$FIRST" -le 1 ] && [ -e "$OUTDIR" ]; then echo "ERROR: $OUTDIR exists; reruns go to a new folder" >&2; exit 1; fi
run() { [ "$FIRST" -le "$1" ] && [ "$LAST" -ge "$1" ]; }
case "$OUTDIR" in
    *" "*) if run 4; then echo "ERROR: OUTDIR contains a space; BLAST+ (step 4) cannot use such paths: $OUTDIR" >&2; exit 1; fi ;;
esac
# refuse before touching STATUS, so a finished run keeps its record
for pair in "2:loci.tsv" "3:overlap.tsv" "4:unclustered/unclustered_families.tsv" "5:rates_by_branch.tsv"; do
    s=${pair%%:*}; f=${pair#*:}
    if run "$s" && [ -e "$OUTDIR/$f" ]; then
        echo "ERROR: $OUTDIR/$f exists (step $s); reruns go to a new folder" >&2; exit 1
    fi
done
if run 4 && compgen -G "$OUTDIR/unclustered/*_unclustered.faa" > /dev/null; then
    echo "ERROR: $OUTDIR/unclustered already holds step 4 files; reruns go to a new folder" >&2; exit 1
fi
mkdir -p "$OUTDIR"
echo "RUNNING steps $FIRST-$LAST $(date)" > "$OUTDIR/STATUS"
fail() {
    echo "FAILED at step ${STEP:-?} $(date)" > "$OUTDIR/STATUS"
    echo "FAILED at step ${STEP:-?}; last lines of $OUTDIR/step${STEP:-?}.log:" >&2
    tail -n 15 "$OUTDIR/step${STEP:-?}.log" >&2 2>/dev/null || true
}
trap fail ERR

STEP=1
PRESENCE_DIR=$TBLASTN_DIR
[ "$PRESENCE" = orf ] && PRESENCE_DIR=$OUTDIR/presence
if run 1; then
    if [ "$PRESENCE" = orf ]; then
        "$PYTHON" "$S/step1a_orf_presence.py" --tree-config "$TREE_CONFIG" --tblastn "$TBLASTN_DIR" \
            --alignments "$ALIGNMENTS_DIR" --stage4-options "$STAGE4_OPTIONS" --outdir "$PRESENCE_DIR" 2> "$OUTDIR/step1a.log"
    fi
    "$PYTHON" "$S/step1_coordinates.py" --tree-config "$TREE_CONFIG" --coords "$COORDS_DIR" --nodes "$NODES_DIR" \
        --tblastn "$PRESENCE_DIR" --out "$OUTDIR/coords.tsv" 2> "$OUTDIR/step1.log"
fi
STEP=2
if run 2; then
    "$PYTHON" "$S/step2_loci.py" --tree-config "$TREE_CONFIG" --coords "$OUTDIR/coords.tsv" --outdir "$OUTDIR" \
        2> "$OUTDIR/step2.log"
fi
STEP=3
if run 3; then
    "$PYTHON" "$S/step3_overlap.py" --tree-config "$TREE_CONFIG" --coords "$OUTDIR/coords.tsv" --gff-table "$GFF_TABLE" \
        --no-gff "$NO_GFF_SPECIES" --out "$OUTDIR/overlap.tsv" 2> "$OUTDIR/step3.log"
fi
STEP=4
if run 4; then
    U=$OUTDIR/unclustered
    "$PYTHON" "$S/step4_unclustered_families.py" select --tree-config "$TREE_CONFIG" --loci "$OUTDIR/loci.tsv" \
        --hogs "$HOGS" --proteins "$PROTEINS_DIR" --outdir "$U" 2> "$OUTDIR/step4.log"
    for faa in "$U"/*_unclustered.faa; do
        [ -e "$faa" ] || continue
        sp=$(basename "$faa" _unclustered.faa)
        "$MAKEBLASTDB" -in "$faa" -dbtype prot -out "$U/${sp}_db" > /dev/null 2>> "$OUTDIR/step4.log"
        "$BLASTP" -query "$faa" -db "$U/${sp}_db" -evalue 1e-5 -max_hsps 1 -num_threads "$THREADS" \
            -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore" -out "$U/${sp}_unclustered.blastp.tsv" \
            2>> "$OUTDIR/step4.log"
        rm -f "$U/${sp}_db".p*
    done
    "$PYTHON" "$S/step4_unclustered_families.py" families --outdir "$U" 2>> "$OUTDIR/step4.log"
fi
STEP=5
if run 5; then
    "$PYTHON" "$S/step5_rates.py" --tree-config "$TREE_CONFIG" --loci "$OUTDIR/loci.tsv" --genes "$OUTDIR/genes.tsv" \
        --hogs "$HOGS" --branch-lengths "$BRANCH_LENGTHS" \
        --overlap "$OUTDIR/overlap.tsv" --exclude-overlap "$EXCLUDE_OVERLAP" \
        --families "$OUTDIR/unclustered/unclustered_families.tsv" --hog-id-match "$HOG_ID_MATCH" --hog-span "$HOG_SPAN" \
        ${UNKNOWN_PRESENCE:+--unknown "$UNKNOWN_PRESENCE"} --outdir "$OUTDIR" 2> "$OUTDIR/step5.log"
fi
echo "SUCCEEDED steps $FIRST-$LAST $(date)" > "$OUTDIR/STATUS"
cat "$OUTDIR/STATUS"
