#!/usr/bin/env bash
# run_pipeline.sh
#
# Stage 4 of LINGUA: enabler detection and DNG/PDNG classification, steps 1-9, from the stage 3 results.
# Which rules apply is set by the method options file (data/method_options.json, or LINGUA_OPTIONS):
#   data/method_options.json                           the default rules (T2T mammal rules of 2026-11-13)
#   presets/method_options_brassicaceae_published.json the rules of the published Brassicaceae catalog (1,811 DNGs)
# The tree (species, node labels, excluded nodes) comes from data/tree_config.json (or LINGUA_CONFIG), written by
# scripts/make_tree_config.py. Final classification: WORKDIR/stage4_final/nodes_summary_by_species/.
#
# Usage:
#   bash run_pipeline.sh WORKDIR [FIRST_STEP] [LAST_STEP]
#
# WORKDIR must contain the upstream inputs described in README.md, namely
#   _all_species_data/                   stage 3 result archives, one .zip per species
#   clsgs_45k_1line/                     candidate de novo gene proteins, one FASTA per species
#   coordinates_alignments_all_species/  MAF block coordinates (only needed for step 9)
# All outputs are written inside WORKDIR. Steps can be rerun individually, for example
# `bash run_pipeline.sh work 4 4` reruns only the tBLASTn step. Set ONLY to a space-separated list of query
# prefixes from data/species.tsv (or SPECIES_TSV) to restrict steps 4, 7c, 7d and 7e. REMOVE_LIST (default
# data/lsgs_to_remove.txt, the loci of the original Brassicaceae stage 2 screens) is applied in step 7b; set
# REMOVE_LIST= (empty) to skip it, as with stage 2 of LINGUA v1.1 or later. REUSE_TBLASTN=1 keeps the step 4
# tBLASTn tables and reruns only step 4b (after a change of the ORF options). CONSERVED_NODES (step 7c) defaults to
# excluded_orf_nodes of the tree config.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$REPO/scripts"
TREE="${TREE:-$REPO/data/SpeciesTree_rooted_node_labels_oct_16.txt}"
SPECIES="${SPECIES_TSV:-$REPO/data/species.tsv}"
REMOVE_LIST="${REMOVE_LIST-$REPO/data/lsgs_to_remove.txt}"
export LINGUA_CONFIG="${LINGUA_CONFIG:-$REPO/data/tree_config.json}"
export LINGUA_OPTIONS="${LINGUA_OPTIONS:-$REPO/data/method_options.json}"
PY="${PYTHON:-python3}"
TBLASTN="${TBLASTN:-tblastn}"
MAKEBLASTDB="${MAKEBLASTDB:-makeblastdb}"
BLASTP="${BLASTP:-$(dirname "$(command -v "$TBLASTN" 2>/dev/null || echo ./blastp)")/blastp}"

# Thresholds used for the published analysis
FRAC_ORF=80      # minimum ORF conservation, % of the query length
PIDENT=40        # minimum tBLASTn identity, %
SYNTENY=40       # minimum aligned fraction of the focal CDS to call synteny, %

WORK="${1:?usage: bash run_pipeline.sh WORKDIR [FIRST_STEP] [LAST_STEP]}"
FIRST="${2:-1}"
LAST="${3:-9}"
WORK="$(cd "$WORK" && pwd)"
mkdir -p "$WORK/logs"

die() { echo "ERROR: $*" >&2; exit 1; }
need() { [ -e "$1" ] || die "missing $1"; }
run() { [ "$1" -ge "$FIRST" ] && [ "$1" -le "$LAST" ]; }
step() { echo; echo "=== step $1: $2 ==="; }
# value of a method option, or the given default; opt: 1 if the option is true, else 0
optval() { ( cd "$SCRIPTS" && "$PY" -c 'import sys, lingua_config as c; print(c.opt(sys.argv[1], sys.argv[2]))' "$1" "$2" ); }
opt() { ( cd "$SCRIPTS" && "$PY" -c 'import sys, lingua_config as c; print(1 if c.opt(sys.argv[1]) else 0)' "$1" ); }
only() { [ -z "${ONLY:-}" ] || [[ " $ONLY " == *" $1 "* ]]; }

for f in "$TREE" "$SPECIES" "$LINGUA_CONFIG"; do need "$f"; done
[ -f "$LINGUA_OPTIONS" ] || die "missing method options $LINGUA_OPTIONS"
[ -z "$REMOVE_LIST" ] || need "$REMOVE_LIST"
echo "method options: $LINGUA_OPTIONS"
echo "tree config:    $LINGUA_CONFIG"
cp "$LINGUA_OPTIONS" "$WORK/logs/method_options_used.json"

RESULTS="$WORK/brassicaceae_23species_LINGUA_results"
NOGAPS="$WORK/brassicaceae_23species_LINGUA_nogaps"
RENAMED="$WORK/45k_renamed"
TBLASTN_OUT="$WORK/cdngs_tblastn_results"
SYNTENY_OUT="$WORK/synteny_${SYNTENY}perc"
RULES="$WORK/stage4_rules"
OPTB="$WORK/stage4_optionB"
FINAL="$WORK/stage4_final"
NODES="$RULES/nodes_summary_by_species"
FINAL_NODES="$FINAL/nodes_summary_by_species"
INPUTS="$WORK/logs/stage4_inputs.tsv"

# per-species inputs (scripts/stage4_inputs.py): focal qpref opref nodes_prefix aln_dir indel_table stop_table ancseq
inputs() {
    ( cd "$SCRIPTS" && "$PY" stage4_inputs.py "$SPECIES" "$RESULTS" "$WORK/indels_dir" "$WORK/stops_dir" "$TBLASTN_OUT" ) > "$INPUTS"
}
ok() { [ "$1" != NA ] || die "$2"; }

if run 1; then
    step 1 "unpack LINGUA results into per-species alignment, indel and stop folders"
    need "$WORK/_all_species_data"
    "$PY" "$SCRIPTS/unzip_and_move_LINGUA_results.py" "$WORK/_all_species_data"
fi

if run 2; then
    step 2 "remove alignment gaps from every per-species sequence"
    need "$RESULTS"
    "$PY" "$SCRIPTS/remove_gaps_batch_recursive.py" -i "$RESULTS" -o "$NOGAPS"
fi

if run 3; then
    step 3 "rename candidate protein files and tabulate gene IDs by species"
    need "$WORK/clsgs_45k_1line"
    "$PY" "$SCRIPTS/rename_file_getIDs.py" \
        "$WORK/clsgs_45k_1line" "$RENAMED" "$WORK/brassicaceae_geneIDs_species_names.tsv"
fi

if run 4; then
    step 4 "tBLASTn of each candidate protein against its gap-free alignment (ORF node)"
    need "$RENAMED"; need "$NOGAPS"
    mkdir -p "$TBLASTN_OUT"
    tail -n +2 "$SPECIES" | while IFS=$'\t' read -r qpref nogaps opref focal; do
        only "$qpref" || continue
        faa=$(ls "$RENAMED"/"${qpref}"_cdngs_*.faa 2>/dev/null | head -1)
        [ -n "$faa" ] || die "no query FASTA for $qpref in $RENAMED"
        need "$NOGAPS/$nogaps"
        n=$(grep -c '^>' "$faa")
        base="$TBLASTN_OUT/${opref}_cdngs_tblastn_${n}"
        if [ -n "${REUSE_TBLASTN:-}" ]; then
            need "${base}_AncSeq.tsv"; echo "  $focal: tBLASTn table reused (REUSE_TBLASTN)"; continue
        fi
        echo "  $focal: $n proteins"
        rm -f "${base}_AncSeq.tblastn.tsv"   # a fresh step 4 table must not be shadowed by an old step 4b backup
        "$PY" "$SCRIPTS/blast_prts-vs-nogaps_AncSeq_v2.py" \
            "$faa" "$NOGAPS/$nogaps" "${base}.tsv" "${base}_summary.tsv" "${base}_AncSeq.tsv" \
            --tree "$TREE" --focal-species "$focal" \
            --frac_ORFss "$FRAC_ORF" --pident "$PIDENT" \
            --tblastn "$TBLASTN" --makeblastdb "$MAKEBLASTDB" --sanity-check
    done
    if [ "$(optval orf_test tblastn)" = length ]; then
        inputs
        while IFS=$'\t' read -r focal qpref opref npref alnd indt stpt anc; do
            only "$qpref" || continue
            [ "$anc" != NA ] || continue
            ok "$alnd" "step 4b: no alignment folder for $focal under $RESULTS"
            step 4b "length-based ORF test for $focal (>= $(optval orf_min_len 80)% of the CDS, >= $(optval orf_min_in_frame 0)% in frame)"
            "$PY" "$SCRIPTS/orf_length_node.py" "$anc" "$alnd" --focal-species "$focal" \
                --min-len "$(optval orf_min_len 80)" --min-cov "$(optval orf_min_cov 80)"
        done < "$INPUTS"
    fi
fi

if run 5; then
    step 5 "synteny conservation per species and the oldest syntenic node"
    need "$NOGAPS"
    "$PY" "$SCRIPTS/summary_synteny_tree_AncSeq.py" \
        --batch "$NOGAPS" --batch-outdir "$SYNTENY_OUT" --tree "$TREE" --synteny "$SYNTENY"
fi

if run 6; then
    need "$WORK/indels_dir"; need "$WORK/stops_dir"; need "$TBLASTN_OUT"
    inputs
    if [ "$(opt start_codon_events)" = 1 ]; then
        step 5b "add start-codon events (non-ATG at the focal ATG) to the stop tables"
        while IFS=$'\t' read -r focal qpref opref npref alnd indt stpt anc; do
            [ "$anc" != NA ] || continue
            ok "$alnd" "step 5b: no alignment folder for $focal"; ok "$stpt" "step 5b: no stop table for $focal"
            "$PY" "$SCRIPTS/add_start_codon_events.py" "$WORK/stops_dir" "$alnd" "$anc" --focal-species "$focal" --table "$stpt"
        done < "$INPUTS"
    elif ls "$WORK"/stops_dir/*.pre_start >/dev/null 2>&1; then
        # start_codon_events off after a run with it on: restore the original stop tables
        for f in "$WORK"/stops_dir/*.pre_start; do cp "$f" "${f%.pre_start}"; done
    fi
    step 6 "assign shared indels and premature stop codons to nodes"
    "$PY" "$SCRIPTS/indel_oldest_shared_node_AncSeq_v2.py" "$WORK/indels_dir"
    "$PY" "$SCRIPTS/stops_oldest_shared_node_AncSeq_v2.py" "$WORK/stops_dir"
    # step 6 just rewrote the event lists; a step 6b backup of an earlier run would otherwise be read instead of them
    rm -f "$WORK"/indels_dir/all_indels_nodes_*.unfiltered.tsv
    if [ "$(opt drop_terminal_indels)" = 1 ] || [ "$(opt homopolymer_needs_extant)" = 1 ]; then
        step 6b "remove indel events that are not evidence (terminal gaps; homopolymers without extant support)"
        while IFS=$'\t' read -r focal qpref opref npref alnd indt stpt anc; do
            [ "$anc" != NA ] || continue
            ok "$alnd" "step 6b: no alignment folder for $focal"; ok "$indt" "step 6b: no indel table for $focal"
            "$PY" "$SCRIPTS/filter_enabler_events.py" "$WORK/indels_dir" "$alnd" --focal-species "$focal" --detail "$indt"
        done < "$INPUTS"
    fi
    if [ "$(opt orf_bearing_nodes_filter)" = 1 ]; then
        # an ancestor that carries the ORF (AncSeq table) cannot be a disabler node; an event touching any
        # ORF-bearing node is dropped whole. The input of the last run is kept as *.pre_orf_bearing.tsv.
        step 6c "drop disabler events that touch an ORF-bearing node (whole event)"
        n=0
        while IFS=$'\t' read -r focal qpref opref npref alnd indt stpt anc; do
            [ "$anc" != NA ] || continue
            ind="$WORK/indels_dir/all_indels_nodes_${focal}.tsv"; stp="$WORK/stops_dir/all_stops_nodes_${focal}.tsv"
            [ -s "$ind" ] && [ -s "$stp" ] || die "step 6c: node tables missing for $focal ($opref)"
            tmp="$WORK/step6c_tmp"; rm -rf "$tmp"; mkdir -p "$tmp"
            "$PY" "$SCRIPTS/step6b_orf_bearing_nodes_filter.py" "$anc" "$ind" "$stp" "$tmp" \
                --log "$WORK/logs/step6c_orf_bearing_${focal}.tsv" --drop-whole-event 2>&1 | sed "s/^/  $focal: /"
            cp "$ind" "${ind%.tsv}.pre_orf_bearing.tsv"; cp "$stp" "${stp%.tsv}.pre_orf_bearing.tsv"
            mv "$tmp/$(basename "$ind")" "$ind"; mv "$tmp/$(basename "$stp")" "$stp"
            rmdir "$tmp"
            n=$((n + 1))
        done < "$INPUTS"
        echo "  step 6c: $n species filtered"
    fi
fi

if run 7; then
    step 7 "combine ORF, synteny and disabler nodes and classify DNGs and PDNGs"
    for d in "$RULES" "$OPTB" "$FINAL"; do [ ! -e "$d" ] || die "$d exists; remove it to rerun step 7"; done
    "$PY" "$SCRIPTS/node_summary_tree_v2.py" --batch \
        --tblastn-dir "$TBLASTN_OUT" --synteny-dir "$SYNTENY_OUT" \
        --indels-dir "$WORK/indels_dir" --stops-dir "$WORK/stops_dir" \
        --tree "$TREE" --outdir "$NODES"
    # node_summary_tree_v2.py names files after the canonical species; step 9 expects the
    # short prefixes defined in nodes_by_coordinates.py, so rename with that same mapping.
    ( cd "$SCRIPTS" && "$PY" -c '
import pathlib, sys
from nodes_by_coordinates import species_to_nodes_prefix
d = pathlib.Path(sys.argv[1])
for f in sorted(d.glob("*_nodes_summary.tsv")):
    sp = f.name[: -len("_nodes_summary.tsv")]
    new = d / f"{species_to_nodes_prefix(sp)}_nodes_summary.tsv"
    if new != f:
        if new.exists():
            sys.exit(f"ERROR: {new} already exists")
        f.rename(new)
        print(f"  renamed {f.name} -> {new.name}")
' "$NODES" )
    if [ -n "$REMOVE_LIST" ]; then
        step 7b "remove loci listed in $(basename "$REMOVE_LIST")"
        "$PY" "$SCRIPTS/filter_nodes_summary_by_geneids.py" "$REMOVE_LIST" "$NODES"
    else
        step 7b "no REMOVE_LIST: every locus kept"
        "$PY" "$SCRIPTS/filter_nodes_summary_by_geneids.py" /dev/null "$NODES"
    fi
    inputs
    if [ "$(optval conserved_orf_filter True)" != False ]; then
        CONSERVED_NODES="${CONSERVED_NODES:-$( cd "$SCRIPTS" && "$PY" -c 'import lingua_config as c; print(" ".join(sorted(c.EXCLUDED_ORF_NODES, key=lambda n: c.NODE_ORDER_YOUNG_TO_OLD.index(n))))' )}"
        step 7c "remove candidates with an ORF in an ancestor outside the ingroup ($CONSERVED_NODES)"
        if [ -n "$CONSERVED_NODES" ]; then
            while IFS=$'\t' read -r focal qpref opref npref alnd indt stpt anc; do
                only "$qpref" || continue
                [ -s "$NODES/${npref}_nodes_summary_filtered.tsv" ] || continue
                ok "$anc" "step 7c: no AncSeq table for $focal"
                "$PY" "$SCRIPTS/remove_conserved_orfs.py" "$NODES" "$anc" "$npref" $CONSERVED_NODES
            done < "$INPUTS"
        fi
    else
        echo "  step 7c skipped (conserved_orf_filter false)"
    fi
    if [ "$(opt remove_literature_etes)" = 1 ]; then
        ete_faa="$(optval ete_proteins '')"
        [ -n "$ete_faa" ] || die "remove_literature_etes needs ete_proteins (a protein FASTA) in the method options"
        case "$ete_faa" in /*) ;; *) ete_faa="$REPO/$ete_faa" ;; esac
        need "$ete_faa"
        step 7e "remove literature ETEs (BLASTp vs $(basename "$ete_faa"), >= 50% identity over >= 50% of the candidate or the ETE)"
        while IFS=$'\t' read -r focal qpref opref npref alnd indt stpt anc; do
            only "$qpref" || continue
            [ -s "$NODES/${npref}_nodes_summary_filtered.tsv" ] || continue
            faa=$(ls "$RENAMED"/"${qpref}"_cdngs_*.faa 2>/dev/null | head -1)
            [ -n "$faa" ] || die "no query FASTA for $qpref in $RENAMED"
            "$PY" "$SCRIPTS/remove_literature_etes.py" "$NODES" "$faa" "$npref" "$ete_faa" --blastp "$BLASTP" --makeblastdb "$MAKEBLASTDB"
        done < "$INPUTS"
    fi
    if [ "$(opt annotate_evidence)" = 1 ]; then
        step 7d "annotate the evidence behind every enabler node"
        while IFS=$'\t' read -r focal qpref opref npref alnd indt stpt anc; do
            only "$qpref" || continue
            [ -s "$NODES/${npref}_nodes_summary_filtered.tsv" ] || continue
            ok "$alnd" "step 7d: no alignment folder for $focal"; ok "$indt" "step 7d: no indel table for $focal"
            ok "$stpt" "step 7d: no stop table for $focal"
            allev=()
            if [ "$(opt evidence_any_valid_disabler)" = 1 ]; then
                allev=(--all-indel-events "$WORK/indels_dir/all_indels_nodes_${focal}.tsv"
                       --all-stop-events "$WORK/stops_dir/all_stops_nodes_${focal}.tsv")
            fi
            "$PY" "$SCRIPTS/annotate_evidence.py" "$NODES" "$npref" "$alnd" "$indt" \
                --focal-species "$focal" --stops-detail "$stpt" ${allev[@]+"${allev[@]}"}
        done < "$INPUTS"
    fi
    # Steps 7f and 7g write parallel folders with only Gene_type rewritten (DNG -> PDNG):
    # stage4_rules/ (before 7f), stage4_optionB/ (after 7f), stage4_final/ (after 7g, the final classification).
    if [ "$(opt option_b_extant_support)" = 1 ] || [ "$(opt consecutive_enablers)" = 1 ]; then
        [ "$(opt annotate_evidence)" = 1 ] || die "option_b_extant_support and consecutive_enablers need annotate_evidence"
    fi
    cur="$RULES"
    if [ "$(opt option_b_extant_support)" = 1 ]; then
        step 7f "option B: a DNG needs an extant carrier at each enabler node"
        "$PY" "$SCRIPTS/apply_option_b_reclassification.py" "$cur" "$OPTB"
        cur="$OPTB"
    fi
    if [ "$(opt consecutive_enablers)" = 1 ]; then
        step 7g "consecutive enablers: the parent and grandparent of the ORF node"
        "$PY" "$SCRIPTS/apply_consecutive_enablers.py" "$cur" "$FINAL"
    else
        mkdir -p "$FINAL"; cp -R "$cur/nodes_summary_by_species" "$FINAL/"
    fi
fi

if run 8; then
    step 8 "summaries of DNGs and PDNGs by node and species"
    "$PY" "$SCRIPTS/summarize_nodes_by_species_split.py" "$FINAL_NODES" "$WORK/brassicaceae" --tree "$TREE"
fi

if run 9; then
    step 9 "DNGs and PDNGs at internal nodes from alignment coordinate intersections"
    need "$WORK/coordinates_alignments_all_species"
    "$PY" "$SCRIPTS/nodes_by_coordinates.py" \
        "$WORK/coordinates_alignments_all_species" "$FINAL_NODES" "$TBLASTN_OUT" "$WORK/parsed_coordinates"
fi

echo; echo "Done (steps $FIRST-$LAST)."
