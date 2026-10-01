#!/usr/bin/env bash
# ==============================================================================
# Script: step1_mafextractor.sh
# Description: Wrapper for mafExtractor. Extracts localized blocks from a master MAF file based on provided coordinates.
# Usage: ./step1_mafextractor.sh <species_name> <reference_species_coordinates_file> <maf_file> <output_dir>
# ==============================================================================

# === MAIN EXECUTION ROUTINE ===
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <species_name> <reference_species_coordinates_file> <maf_file> <output_dir>"
    exit 1
fi

species_name="$1"
reference_species_coordinates_file="$2"
maf_file="$3"
output_dir="$4"

mkdir -p "$output_dir"

# Remove Windows carriage returns just in case
sed -i 's/\r$//' "$reference_species_coordinates_file"

# --- Execution Loop ---
tail -n +2 "$reference_species_coordinates_file" | while IFS=$'\t' read -r seq start stop; do
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

    echo "Running mafExtractor for $seq:$start-$stop -> $output_path"

    mafExtractor --maf "$maf_file" --seq "$seq" --start "$start" --stop "$stop" > "$output_path"
done

echo "All extractions completed."
