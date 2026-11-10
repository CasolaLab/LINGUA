#!/usr/bin/env bash
# run_pipeline.sh
#
# Runs the enabler-detection and node-assignment pipeline, steps 1-9, on the LINGUA
# outputs of the 23-species Progressive Cactus alignment. Steps 6 and 7 include the rules
# of the NEE revision (2026-10-28/29) that produced the published catalog of 1,811 DNGs,
# namely step 6b (ORF-bearing nodes cannot carry disablers; whole event dropped), the step 7
# classification with orf_clade_consistency and recorded enabler events (data/method_options.json),
# step 7c (annotate_evidence.py), step 7d (option B) and step 7e (consecutive enablers).
# Final classification: WORKDIR/stage4_final/nodes_summary_by_species/.
#
# Usage:
#   bash run_pipeline.sh WORKDIR [FIRST_STEP] [LAST_STEP]
#
# WORKDIR must contain the upstream inputs described in README.md, namely
#   _all_species_data/                   LINGUA result archives, one .zip per species
#   clsgs_45k_1line/                     candidate de novo gene proteins, one FASTA per species
#   coordinates_alignments_all_species/  MAF block coordinates (only needed for step 9)
# All outputs are written inside WORKDIR. Steps can be rerun individually, for
# example `bash run_pipeline.sh work 4 4` reruns only the tBLASTn step. Set ONLY to a
# space-separated list of query prefixes from data/species.tsv to restrict step 4, for
# example `ONLY=athaliana bash run_pipeline.sh work`; the other steps process whatever
# species are present in WORKDIR. Set REMOVE_LIST= (empty) to skip step 7b, the removal of
# loci flagged by the original stage 2 domain and TE screens (data/lsgs_to_remove.txt); with
# stage 2 of LINGUA v1.1 or later those loci never reach stage 3.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$REPO/scripts"
TREE="$REPO/data/SpeciesTree_rooted_node_labels_oct_16.txt"
SPECIES="$REPO/data/species.tsv"
REMOVE_LIST="${REMOVE_LIST-$REPO/data/lsgs_to_remove.txt}"
export LINGUA_CONFIG="${LINGUA_CONFIG:-$REPO/data/tree_config.json}"
export LINGUA_OPTIONS="${LINGUA_OPTIONS:-$REPO/data/method_options.json}"
PY="${PYTHON:-python3}"
TBLASTN="${TBLASTN:-tblastn}"
MAKEBLASTDB="${MAKEBLASTDB:-makeblastdb}"

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

for f in "$TREE" "$SPECIES" "$LINGUA_CONFIG" "$LINGUA_OPTIONS"; do need "$f"; done
[ -z "$REMOVE_LIST" ] || need "$REMOVE_LIST"

NOGAPS="$WORK/brassicaceae_23species_LINGUA_nogaps"
RENAMED="$WORK/45k_renamed"
TBLASTN_OUT="$WORK/cdngs_tblastn_results"
SYNTENY_OUT="$WORK/synteny_${SYNTENY}perc"
RULES="$WORK/stage4_rules"
OPTB="$WORK/stage4_optionB"
FINAL="$WORK/stage4_final"
NODES="$RULES/nodes_summary_by_species"
FINAL_NODES="$FINAL/nodes_summary_by_species"

if run 1; then
    step 1 "unpack LINGUA results into per-species alignment, indel and stop folders"
    need "$WORK/_all_species_data"
    "$PY" "$SCRIPTS/unzip_and_move_LINGUA_results.py" "$WORK/_all_species_data"
fi

if run 2; then
    step 2 "remove alignment gaps from every per-species sequence"
    need "$WORK/brassicaceae_23species_LINGUA_results"
    "$PY" "$SCRIPTS/remove_gaps_batch_recursive.py" \
        -i "$WORK/brassicaceae_23species_LINGUA_results" -o "$NOGAPS"
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
        if [ -n "${ONLY:-}" ] && [[ " $ONLY " != *" $qpref "* ]]; then continue; fi
        faa=$(ls "$RENAMED"/"${qpref}"_cdngs_*.faa 2>/dev/null | head -1)
        [ -n "$faa" ] || die "no query FASTA for $qpref in $RENAMED"
        need "$NOGAPS/$nogaps"
        n=$(grep -c '^>' "$faa")
        base="$TBLASTN_OUT/${opref}_cdngs_tblastn_${n}"
        echo "  $focal: $n proteins"
        "$PY" "$SCRIPTS/blast_prts-vs-nogaps_AncSeq_v2.py" \
            "$faa" "$NOGAPS/$nogaps" "${base}.tsv" "${base}_summary.tsv" "${base}_AncSeq.tsv" \
            --tree "$TREE" --focal-species "$focal" \
            --frac_ORFss "$FRAC_ORF" --pident "$PIDENT" \
            --tblastn "$TBLASTN" --makeblastdb "$MAKEBLASTDB" --sanity-check
    done
fi

if run 5; then
    step 5 "synteny conservation per species and the oldest syntenic node"
    need "$NOGAPS"
    "$PY" "$SCRIPTS/summary_synteny_tree_AncSeq.py" \
        --batch "$NOGAPS" --batch-outdir "$SYNTENY_OUT" --tree "$TREE" --synteny "$SYNTENY"
fi

if run 6; then
    step 6 "assign shared indels and premature stop codons to nodes"
    need "$WORK/indels_dir"; need "$WORK/stops_dir"
    "$PY" "$SCRIPTS/indel_oldest_shared_node_AncSeq_v2.py" "$WORK/indels_dir"
    "$PY" "$SCRIPTS/stops_oldest_shared_node_AncSeq_v2.py" "$WORK/stops_dir"
    step 6b "drop disabler events that touch an ORF-bearing node (pre-filter tables kept as *_unfiltered.tsv)"
    need "$TBLASTN_OUT"
    n=0
    for anc in "$TBLASTN_OUT"/*_cdngs_tblastn_*_AncSeq.tsv; do
        pre=$(basename "$anc" | sed 's/_cdngs_tblastn.*//')
        sp=$(cd "$SCRIPTS" && "$PY" -c 'import sys, lingua_config as c; print(c.canonicalize(sys.argv[1]))' "$pre")
        ind="$WORK/indels_dir/all_indels_nodes_${sp}.tsv"; stp="$WORK/stops_dir/all_stops_nodes_${sp}.tsv"
        [ -s "$ind" ] && [ -s "$stp" ] || die "step 6b: node tables missing for $sp ($pre)"
        tmp="$WORK/step6b_tmp"; mkdir -p "$tmp" "$WORK/logs"
        "$PY" "$SCRIPTS/step6b_orf_bearing_nodes_filter.py" "$anc" "$ind" "$stp" "$tmp" \
            --log "$WORK/logs/step6b_${sp}.tsv" --drop-whole-event | sed "s/^/  $sp: /"
        mv "$ind" "${ind%.tsv}_unfiltered.tsv"; mv "$stp" "${stp%.tsv}_unfiltered.tsv"
        mv "$tmp/$(basename "$ind")" "$ind"; mv "$tmp/$(basename "$stp")" "$stp"
        n=$((n + 1))
    done
    rmdir "$WORK/step6b_tmp"
    echo "  step 6b: $n species filtered"
fi

if run 7; then
    step 7 "combine ORF, synteny and disabler nodes and classify DNGs and PDNGs"
    [ ! -e "$RULES" ] || die "$RULES exists; remove it to rerun step 7"
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
        step 7b "remove loci flagged by the original stage 2 domain and TE screens"
        "$PY" "$SCRIPTS/filter_nodes_summary_by_geneids.py" "$REMOVE_LIST" "$NODES"
    fi
    step 7c "annotate the evidence behind every enabler node (any valid disabler at the node)"
    ( cd "$SCRIPTS" && "$PY" - "$WORK" "$NODES" ) > "$WORK/logs/step7c_inputs.tsv" <<'EOF'
import csv, glob, os, sys
import lingua_config as cfg
csv.field_size_limit(10 ** 9)
wd, nodes = sys.argv[1:]
aln = {}
for d in glob.glob(os.path.join(wd, "brassicaceae_23species_LINGUA_results", "*")):
    fs = glob.glob(os.path.join(d, "*.aln"))
    sps = {sp for sp in cfg.SPECIES for f in fs[:50] if os.path.basename(f).startswith(sp + "_")}
    if len(sps) == 1:
        aln[sps.pop()] = d
ind, stp = {}, {}
for sub, pat, store in (("indels_dir", "*_indels_summary.tsv", ind), ("stops_dir", "*_stop_codon_summary.tsv", stp)):
    for f in glob.glob(os.path.join(wd, sub, pat)):
        r = next(csv.DictReader(open(f), delimiter="\t"), None)
        hit = [sp for sp in cfg.SPECIES if r and r["file_name"].startswith(sp + "_")]
        if len(hit) == 1:
            store[hit[0]] = f
for f in sorted(glob.glob(os.path.join(nodes, "*_nodes_summary_filtered.tsv"))):
    pre = os.path.basename(f)[: -len("_nodes_summary_filtered.tsv")]
    sp = cfg.canonicalize(pre)
    if sp not in aln or sp not in ind or sp not in stp:
        sys.exit(f"ERROR: step 7c inputs missing for {pre} ({sp})")
    print("\t".join([pre, sp, aln[sp], ind[sp], stp[sp]]))
EOF
    while IFS=$'\t' read -r pre sp alnd indf stpf; do
        "$PY" "$SCRIPTS/annotate_evidence.py" "$NODES" "$pre" "$alnd" "$indf" \
            --focal-species "$sp" --stops-detail "$stpf" \
            --all-indel-events "$WORK/indels_dir/all_indels_nodes_${sp}.tsv" \
            --all-stop-events "$WORK/stops_dir/all_stops_nodes_${sp}.tsv"
    done < "$WORK/logs/step7c_inputs.tsv"
    step 7d "option B: a DNG needs extant support for its enabler nodes"
    "$PY" "$SCRIPTS/apply_option_b_reclassification.py" "$RULES" "$OPTB"
    step 7e "consecutive enablers: the parent and grandparent of the ORF node"
    "$PY" "$SCRIPTS/apply_consecutive_enablers.py" "$OPTB" "$FINAL"
fi

if run 8; then
    step 8 "family-wide summaries of DNGs and PDNGs by node and species"
    "$PY" "$SCRIPTS/summarize_nodes_by_species_split.py" "$FINAL_NODES" "$WORK/brassicaceae" --tree "$TREE"
fi

if run 9; then
    step 9 "DNGs and PDNGs at internal nodes from alignment coordinate intersections"
    need "$WORK/coordinates_alignments_all_species"
    "$PY" "$SCRIPTS/nodes_by_coordinates.py" \
        "$WORK/coordinates_alignments_all_species" "$FINAL_NODES" "$TBLASTN_OUT" "$WORK/parsed_coordinates"
fi

echo; echo "Done (steps $FIRST-$LAST)."
