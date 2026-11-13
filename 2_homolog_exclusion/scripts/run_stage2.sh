#!/usr/bin/env bash
# run_stage2.sh - LINGUA stage 2 end to end: the six screens, combine, the HOG rule and the paralog rule.
#
# Usage: bash 2_homolog_exclusion/scripts/run_stage2.sh stage2.conf [first_step] [last_step]
#   stage2.conf is a bash file of variables (see presets/stage2_run.conf.example):
#     CANDIDATES      folder of candidate protein FASTAs, one per species (stage 1)
#     OUTDIR          output folder (must not exist when starting at step 1)
#     BLAST_DBS       protein databases for BLASTp, space-separated (each searched as its own run, then pooled)
#     JACKHMMER_DB    protein FASTA for jackhmmer
#     CDD_DB          CDD database for RPS-BLAST
#     TE_DB, TE_CONFIG              TE protein database and its settings file
#     ORGANELLE_DB, ORGANELLE_CONFIG  organelle genome nucleotide database and its settings file
#     HOGS            OrthoFinder HOG table (N0.tsv, or the node holding all species)
#     PROTEOMES       full proteomes, one FASTA per species, named like the candidate files
#     CONFIG          settings file for the homology and domain searches and the two rules (optional; defaults otherwise)
#     THREADS         threads per search (default 4); SPECIES_MAP, IPRSCAN_BIN, PANTHER_FAMILY_SPECIES optional
#     SKIP_SCREENS    optional comma-separated screens left out (blast, jackhmmer, cdd, interpro, te, organelle);
#                     recorded in STATUS. The HOG rule and the paralog rule cannot be left out.
# Steps: 1 BLASTp, 2 jackhmmer, 3 CDD, 4 InterProScan, 5 TE, 6 organelle, 7 combine, 8 HOG rule, 9 paralog rule.
# The stage 2 output (input of stage 3) is OUTDIR/paralog/passed/ and OUTDIR/paralog/passed_ids.tsv.
# A module that leaves any gene unsearched exits with status 2; the run then stops as FAILED.
set -Eeuo pipefail

[ $# -ge 1 ] || { sed -n '2,23p' "$0"; exit 1; }
CONF=$1; FIRST=${2:-1}; LAST=${3:-9}
[ -f "$CONF" ] || { echo "ERROR: config not found: $CONF" >&2; exit 1; }
case "$FIRST$LAST" in
    [1-9][1-9]) [ "$FIRST" -le "$LAST" ] || { echo "ERROR: first step $FIRST is after last step $LAST" >&2; exit 1; } ;;
    *) echo "ERROR: steps must be numbers 1-9 (got '$FIRST' '$LAST')" >&2; exit 1 ;;
esac
# shellcheck disable=SC1090
source "$CONF"
HERE=$(cd "$(dirname "$0")" && pwd)
P() { bash "$HERE/homolog_exclusion_pipeline.sh" "$@"; }
THREADS=${THREADS:-4}; SKIP_SCREENS=${SKIP_SCREENS:-}; CONFIG=${CONFIG:-}
SPECIES_MAP=${SPECIES_MAP:-}; IPRSCAN_BIN=${IPRSCAN_BIN:-}; PANTHER_FAMILY_SPECIES=${PANTHER_FAMILY_SPECIES:-}

skipped() { case ",$SKIP_SCREENS," in *",$1,"*) return 0 ;; *) return 1 ;; esac; }
for s in ${SKIP_SCREENS//,/ }; do
    case "$s" in blast|jackhmmer|cdd|interpro|te|organelle) ;;
        *) echo "ERROR: SKIP_SCREENS: unknown screen '$s' (the HOG and paralog rules cannot be skipped)" >&2; exit 1 ;;
    esac
done
need=(CANDIDATES OUTDIR HOGS PROTEOMES)
skipped blast || need+=(BLAST_DBS)
skipped jackhmmer || need+=(JACKHMMER_DB)
skipped cdd || need+=(CDD_DB)
skipped te || need+=(TE_DB TE_CONFIG)
skipped organelle || need+=(ORGANELLE_DB ORGANELLE_CONFIG)
for v in "${need[@]}"; do
    [ -n "${!v:-}" ] || { echo "ERROR: $v is not set in $CONF" >&2; exit 1; }
done
for v in CANDIDATES HOGS PROTEOMES; do
    [ -e "${!v}" ] || { echo "ERROR: $v not found: ${!v}" >&2; exit 1; }
done
if [ "$FIRST" -le 1 ] && [ -e "$OUTDIR" ]; then echo "ERROR: $OUTDIR exists; reruns go to a new folder" >&2; exit 1; fi
mkdir -p "$OUTDIR"
if [ -f "$OUTDIR/STATUS" ] && grep -q '^RUNNING' "$OUTDIR/STATUS"; then
    echo "ERROR: $OUTDIR/STATUS says RUNNING; another run may be writing here" >&2; exit 1
fi
run() { [ "$FIRST" -le "$1" ] && [ "$LAST" -ge "$1" ]; }
SKIPNOTE=""; [ -n "$SKIP_SCREENS" ] && SKIPNOTE=" (screens left out: $SKIP_SCREENS)"
echo "RUNNING steps $FIRST-$LAST$SKIPNOTE $(date)" > "$OUTDIR/STATUS"
fail() {
    echo "FAILED at step ${STEP:-?}$SKIPNOTE $(date)" > "$OUTDIR/STATUS"
    echo "FAILED at step ${STEP:-?} (exit status 2 means some genes were not searched); last lines of $OUTDIR/step${STEP:-?}.log:" >&2
    tail -n 15 "$OUTDIR/step${STEP:-?}.log" >&2 2>/dev/null || true
}
trap fail ERR

COMMON=(--threads "$THREADS")
[ -n "$SPECIES_MAP" ] && COMMON+=(--species-map "$SPECIES_MAP")
CFG=(); [ -n "$CONFIG" ] && CFG=(--config "$CONFIG")

screen() {  # screen <step> <name> <command...>
    local n=$1 name=$2; shift 2
    STEP=$n
    if run "$n"; then
        if skipped "$name"; then
            echo "screen $name left out (SKIP_SCREENS)" > "$OUTDIR/step$n.log"
        else
            "$@" > "$OUTDIR/step$n.log" 2>&1
        fi
    fi
}

DBA=(); for d in ${BLAST_DBS:-}; do DBA+=(--db "$d"); done
screen 1 blast P homology-blast -i "$CANDIDATES" -o "$OUTDIR/blast" --label blast ${CFG[@]+"${CFG[@]}"} \
    ${DBA[@]+"${DBA[@]}"} "${COMMON[@]}"
screen 2 jackhmmer P homology-jackhmmer -i "$CANDIDATES" -o "$OUTDIR/jackhmmer" --label jackhmmer ${CFG[@]+"${CFG[@]}"} \
    --db "${JACKHMMER_DB:-}" "${COMMON[@]}"
screen 3 cdd P domain-cdd -i "$CANDIDATES" -o "$OUTDIR/cdd" --label cdd ${CFG[@]+"${CFG[@]}"} --db "${CDD_DB:-}" "${COMMON[@]}"
IPR=(); [ -n "$IPRSCAN_BIN" ] && IPR+=(--iprscan-bin "$IPRSCAN_BIN")
[ -n "$PANTHER_FAMILY_SPECIES" ] && IPR+=(--panther-family-species "$PANTHER_FAMILY_SPECIES")
screen 4 interpro P domain-interpro -i "$CANDIDATES" -o "$OUTDIR/interpro" --label interpro ${CFG[@]+"${CFG[@]}"} \
    ${IPR[@]+"${IPR[@]}"} "${COMMON[@]}"
screen 5 te P homology-blast -i "$CANDIDATES" -o "$OUTDIR/te" --label te --config "${TE_CONFIG:-}" --db "${TE_DB:-}" "${COMMON[@]}"
screen 6 organelle P homology-blast -i "$CANDIDATES" -o "$OUTDIR/organelle" --label organelle --config "${ORGANELLE_CONFIG:-}" \
    --db "${ORGANELLE_DB:-}" "${COMMON[@]}"

STEP=7
if run 7; then
    RUNS=()
    for s in blast jackhmmer cdd interpro te organelle; do skipped "$s" || RUNS+=("$OUTDIR/$s"); done
    P combine --runs "${RUNS[@]}" -o "$OUTDIR/combine" > "$OUTDIR/step7.log" 2>&1
fi
STEP=8
if run 8; then
    EV=()
    skipped blast || EV+=("blast*")
    skipped jackhmmer || EV+=("jackhmmer*")
    skipped cdd || EV+=("cdd")
    skipped interpro || EV+=("interpro")
    [ ${#EV[@]} -gt 0 ] || { echo "ERROR: no homology or domain screen left for the HOG rule" > "$OUTDIR/step8.log"; false; }
    EVS=$(IFS=,; echo "${EV[*]}")
    P hog-rule -i "$CANDIDATES" -o "$OUTDIR/hog" --combine "$OUTDIR/combine" --hogs "$HOGS" --evidence-runs "$EVS" \
        ${CFG[@]+"${CFG[@]}"} "${COMMON[@]}" > "$OUTDIR/step8.log" 2>&1
fi
STEP=9
if run 9; then
    P paralog-rule -i "$CANDIDATES" -o "$OUTDIR/paralog" --passed "$OUTDIR/hog" --proteomes "$PROTEOMES" \
        ${CFG[@]+"${CFG[@]}"} "${COMMON[@]}" > "$OUTDIR/step9.log" 2>&1
fi
echo "SUCCEEDED steps $FIRST-$LAST$SKIPNOTE $(date)" > "$OUTDIR/STATUS"
cat "$OUTDIR/STATUS"
