#!/usr/bin/env bash
# ==============================================================================
# Script: step1_mafextractor.sh
# Description: Wrapper for mafExtractor. Extracts localized blocks from a master MAF file based on provided coordinates.
# Usage: ./step1_mafextractor.sh <species_name> <reference_species_coordinates_file> <maf_file> <output_dir>
#
# Each exon is written to a temporary file and only renamed to <species>.<seq>_<start>_<stop>.txt when
# mafExtractor exits successfully and the output contains at least one reference ("s <species>.") line.
# Failed exons are listed on stderr and in <output_dir>/step1_failed.tsv, and the script exits with status 1,
# so a missing mafExtractor, an unreadable MAF or a wrong sequence name can no longer leave empty files that
# later steps silently process. The coordinate file is read, not modified.
# ==============================================================================
set -uo pipefail

# === MAIN EXECUTION ROUTINE ===
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <species_name> <reference_species_coordinates_file> <maf_file> <output_dir>"
    exit 1
fi

species_name="$1"
reference_species_coordinates_file="$2"
maf_file="$3"
output_dir="$4"

command -v mafExtractor >/dev/null || { echo "Error: mafExtractor not found on PATH" >&2; exit 1; }
[ -r "$maf_file" ] || { echo "Error: MAF not readable: $maf_file" >&2; exit 1; }
[ -r "$reference_species_coordinates_file" ] || { echo "Error: coordinate file not readable: $reference_species_coordinates_file" >&2; exit 1; }
mkdir -p "$output_dir"
failed_list="$output_dir/step1_failed.tsv"
n_ok=0; n_fail=0

# --- Execution Loop (Windows carriage returns removed while reading) ---
while IFS=$'\t' read -r seq start stop; do
    seq=$(echo "$seq" | xargs)
    start=$(echo "$start" | xargs)
    stop=$(echo "$stop" | xargs)

    # Skip bad lines
    if ! [[ "$start" =~ ^[0-9]+$ ]] || ! [[ "$stop" =~ ^[0-9]+$ ]]; then
        echo "Skipping invalid line: $seq $start $stop" >&2
        continue
    fi

    seq="${species_name}.$seq"
    base_name="${seq}_${start}_${stop}.txt"
    output_path="${output_dir}/${base_name}"
    tmp_path="${output_path}.partial"

    echo "Running mafExtractor for $seq:$start-$stop -> $output_path"

    if ! mafExtractor --maf "$maf_file" --seq "$seq" --start "$start" --stop "$stop" > "$tmp_path"; then
        reason="mafExtractor exited with an error"
    elif ! grep -q "^s ${species_name}\." "$tmp_path"; then
        reason="output has no reference (s ${species_name}.) lines"
    else
        mv "$tmp_path" "$output_path"
        n_ok=$((n_ok + 1))
        continue
    fi
    rm -f "$tmp_path"
    echo -e "${base_name}\t${reason}" >> "$failed_list"
    echo "FAILED: ${base_name}: ${reason}" >&2
    n_fail=$((n_fail + 1))
done < <(tail -n +2 "$reference_species_coordinates_file" | tr -d '\r')

echo "Extractions: $n_ok written, $n_fail failed."
if [ "$n_fail" -gt 0 ]; then
    echo "Error: $n_fail exon(s) failed; see $failed_list" >&2
    exit 1
fi
echo "All extractions completed."
