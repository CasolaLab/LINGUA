#!/usr/bin/env bash
# run_stage3.sh
#
# LINGUA stage 3, steps 1-17 plus step 8c, for one reference species. These are the commands and
# the step order that produced the alignments behind the published catalog (NEE revision,
# 2026-10-28/29): step 4 fills the bases between MAF blocks from each genome, step 8c keeps one
# collinear row per species and writes the missing-data mask, and step 14 reads that mask so that
# missing exons are not called deletions (deletions overlapping the mask are trimmed to their
# unmasked part).
#
# Usage:
#   TREE_FILE=<newick> bash run_stage3.sh <workdir> <Species> <bed> <cds_fasta> <maf> <genome_dir>
#     <Species>    name used in the MAF and genome folder (e.g. Arabidopsis_thaliana)
#     <bed>        reference CDS exons (BED6 + transcript ID in column 7)
#     <cds_fasta>  reference CDS sequences, one per transcript ID in the BED
#     <maf>        Progressive Cactus MAF with <Species> as reference (hal2maf --onlyOrthologs)
#     <genome_dir> <species>.fa genome FASTA files (with .fai), named as in the MAF
# Optional environment:
#   THREADS=<n>     default 4
#   STEP1_DIR=<dir> reuse existing step 1 files (named <Species>.<chrom>_<start>_<end>.txt)
#   JOBS=<n>        parallel mafExtractor chunks when step 1 is run (default THREADS)
# Requires mafExtractor, samtools, GNU parallel and Python 3 with Biopython and pandas on PATH.
# Writes STATUS (RUNNING / SUCCEEDED / FAILED at <step>) in <workdir>; exits non-zero on failure.
set -uo pipefail

WD="${1:?usage: run_stage3.sh <workdir> <Species> <bed> <cds_fasta> <maf> <genome_dir>}"
SP="${2:?Species}"; BED_SRC="${3:?bed}"; CDS_SRC="${4:?cds_fasta}"; MAF="${5:?maf}"; G="${6:?genome_dir}"
S="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
T=${THREADS:-4}; JOBS=${JOBS:-$T}
R=$SP; B=$R.bed
: "${TREE_FILE:?set TREE_FILE to the Newick species tree (read by step 8)}"
export TREE_FILE
for f in "$BED_SRC" "$CDS_SRC" "$TREE_FILE" "$G"; do [ -e "$f" ] || { echo "missing input $f"; exit 2; }; done
[ -n "${STEP1_DIR:-}" ] || [ -e "$MAF" ] || { echo "missing MAF $MAF"; exit 2; }
[ -z "${STEP1_DIR:-}" ] || [ -d "$STEP1_DIR" ] || { echo "missing STEP1_DIR $STEP1_DIR"; exit 2; }
mkdir -p "$WD" && cd "$WD" || exit 2
[ -s STATUS ] && grep -q -E '^(RUNNING|SUCCEEDED)' STATUS && { echo "STATUS says $(cat STATUS); refusing"; exit 3; }
CURRENT=setup
finish() { rc=$?; if [ "$rc" -eq 0 ] && [ "$CURRENT" = done ]; then echo "SUCCEEDED $(date -Is)" > STATUS
           else echo "FAILED at $CURRENT (exit $rc) $(date -Is)" > STATUS; fi; cat STATUS; }
trap finish EXIT
step() { CURRENT="$1"; shift; echo "== $CURRENT $(date -Is)"; "$@" || { echo "step $CURRENT failed"; exit 1; }; }
date -Is > START; echo RUNNING > STATUS

prepare() {
    # CRLF and a trailing '+' on contig names (mafExtractor names lack it) are removed from the BED
    tr -d '\r' < "$BED_SRC" | awk 'BEGIN{FS=OFS="\t"} !/^#/{sub(/\+$/,"",$1)} {print}' > "$B"
    python3 - "$CDS_SRC" "$SP" <<'PY'
import sys
cds_src, sp = sys.argv[1:3]
exons, tx = [], []
for l in open(sp + ".bed"):
    l = l.rstrip("\n")
    if l.startswith("#") or not l:
        continue
    f = l.split("\t")
    if (f[0], f[1], f[2]) not in exons:
        exons.append((f[0], f[1], f[2]))
    if f[6] not in tx:
        tx.append(f[6])
with open("coords.tsv", "w") as out:
    out.write("seq\tstart\tstop\n")
    for e in exons:
        out.write("\t".join(e) + "\n")
open("N_EXONS", "w").write(str(len(exons)))
seqs, k = {}, None
for l in open(cds_src):
    l = l.strip()
    if l.startswith(">"):
        k = l[1:].split()[0]; seqs[k] = []
    elif k:
        seqs[k].append(l)
for t in tx:   # Phytozome CDS IDs carry a ".p" suffix that the BED IDs lack
    if t not in seqs and t + ".p" in seqs:
        seqs[t] = seqs[t + ".p"]
absent = [t for t in tx if t not in seqs]
assert not absent, f"{len(absent)} transcripts missing from the CDS FASTA, e.g. {absent[:3]}"
with open("ref_raw.fna", "w") as out:
    for t in tx:
        out.write(f">{t}\n{''.join(seqs[t])}\n")
print(f"prepared: {len(exons)} exons for {len(tx)} transcripts")
PY
}

step1() {
    mkdir -p step1
    if [ -n "${STEP1_DIR:-}" ]; then
        tail -n +2 coords.tsv | while IFS=$'\t' read -r c s e; do
            f="$SP.${c}_${s}_${e}.txt"
            [ -f "$STEP1_DIR/$f" ] && cp -p "$STEP1_DIR/$f" step1/ || echo "MISSING_STEP1 $f"
        done > step1_missing.txt
    else
        mkdir -p chunks
        tail -n +2 coords.tsv | split -n l/"$JOBS" -d -a 3 - chunks/c_
        pids=()
        for c in chunks/c_???; do
            { printf 'seq\tstart\tstop\n'; cat "$c"; } > "$c.tsv"
            bash "$S/step1_mafextractor.sh" "$SP" "$c.tsv" "$MAF" step1 > "$c.log" 2>&1 &
            pids+=($!)
        done
        fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
        [ "$fail" -eq 0 ] || { echo "a step 1 chunk failed (see chunks/*.log)"; return 1; }
        : > step1_missing.txt
    fi
    n=$(find step1 -name '*.txt' | wc -l); want=$(cat N_EXONS)
    echo "step 1: $n files for $want exons; $(wc -l < step1_missing.txt) missing"
    [ "$n" -eq "$want" ] || { echo "step 1 file count $n != $want (see step1_missing.txt)"; return 1; }
}

V=step8c/master_summary_final_validation.tsv
step prepare prepare
step step1 step1
step clean_headers python3 "$S/clean_nucleotide_fasta_headers.py" -f ref_raw.fna -b "$B" -o ref.fna
step step2 bash "$S/step2_missingbp_summaries.sh" step1 step2
step step3 python3 "$S/step3_filtering_species.py" -i step1 -o step3 -r "$R" -s step3_summary.csv -t "$T"
step step4 bash "$S/step4_reverse_complement_add_missingbp.sh" -m step3 -g "$G" -o step4 -b "$B" -c step2/combined_all_species.csv -t "$T"
step step5 python3 "$S/step5_validate_reverse_complement.py" step4 step1 "$B" "$R" "$T"
step step6 python3 "$S/step6_update_summary_files.py" step4/extraction_summary.tsv "$B" "$R"
step step7 python3 "$S/step7_first_stage_validation.py" step1 step3 step4 "$R"
step step8 bash "$S/step8_stitching_exons.sh" step4 "$B" step8 "$R"
step step8c python3 "$S/step8c_one_collinear_row.py" . step8c step8c_log.tsv
step step9 python3 "$S/step9_trimming_duplicate_bases.py" step8c step9 ref.fna "$B" "$T"
step step10 python3 "$S/step10_second_stage_validation.py" step9 ref.fna "$R"
step step11 python3 "$S/step11_update_stitched_exons_summary.py" step8c/master_summary.tsv step9 "$B" ref.fna "$R"
step step12 bash "$S/step12_identify_start_codons.sh" step9 step12 "$B" "$R" "$V"
step step13 bash "$S/step13_identify_stop_codons.sh" step9 step13 "$B" "$R" "$V"
step step14 bash "$S/step14_identify_indels.sh" step9 step14 "$B" "$R" "$V" step8c/step8c_missing_data.tsv
step step15 bash "$S/step15_identify_first_and_last_three_bases.sh" step9 step15 "$R" "$V"
step step16 python3 "$S/step16_extract_start_and_stops.py" step9 step16_starts_stops.tsv
step step17 python3 "$S/step17_combine_all_species_info.py" step9 step3 "$B" "$R" step17_all_species.tsv
CURRENT=done
echo "== done $(date -Is)"
