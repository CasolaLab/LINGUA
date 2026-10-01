#!/usr/bin/env bash
# ==============================================================================
# Script: step15_identify_first_and_last_three_bases.sh
# Description: Targeted sequence boundary extractor. Pulls the exact first and last 3 base pairs for a specific target species from the multiple sequence alignment to verify translation boundary integrity.
# Usage: ./step15_identify_first_and_last_three_bases.sh <step9_trimmed_dir> <output_dir> <reference_species> <step11_summary_tsv>
# ==============================================================================

set -euo pipefail
LC_ALL=C

if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <step9_trimmed_dir> <output_dir> <reference_species> <step11_summary_tsv>"
  exit 1
fi

step9_trimmed_dir="$1"
output_dir="$2"
reference_species="$3"
exons_summary="$4"

mkdir -p "$output_dir"
summary="$output_dir/master_first_last.tsv"

# Changed header to use actual Tabs
echo -e "file_name\tgff_transcriptID\tfirst_3\tlast_3" > "$summary"

# -------------------------------------------------
# Load TSV exons summary into an associative array
# -------------------------------------------------
declare -A file2tid

# Changed IFS to safely read Tab-Separated Values (TSV)
while IFS=$'\t' read -r file_name gff_tid _rest; do
    # Skip header
    [[ "$file_name" =~ ^file_name$ ]] && continue

    # Normalize file_name and gff_tid
    file_name="${file_name//\"/}"       # remove quotes
    file_name="${file_name//$'\r'/}"    # remove carriage returns
    file_name="${file_name// /}"        # remove spaces

    gff_tid="${gff_tid//\"/}"
    gff_tid="${gff_tid//$'\r'/}"
    gff_tid="${gff_tid// /}"

    if [[ -n "$file_name" && -n "$gff_tid" ]]; then
        file2tid["$file_name"]="$gff_tid"
    fi
done < "$exons_summary"

# Process ALN files
find "$step9_trimmed_dir" -maxdepth 1 -type f -name "*.aln" | while IFS= read -r aln_file; do
    filename="$(basename "$aln_file")"
    file_name="${filename%.aln}"

    gff_tid="${file2tid[$file_name]:-N/A}"
    if [[ "$gff_tid" == "N/A" ]]; then
        echo "WARNING: No transcript ID found in exons summary for $file_name" >&2
        continue
    fi

    # -------------------------------------------------
    # Sequence extraction (Updated to output TSV)
    # -------------------------------------------------
    awk -v OFS="\t" \
        -v file_name="$file_name" \
        -v gff_tid="$gff_tid" \
        -v reference_species="$reference_species" '
    function clean_and_extract(s,   L, first3, last3) {
      gsub(/-/, "", s)
      s = toupper(s)
      gsub(/[^ACGT]/, "", s)
      L = length(s)
      if (L >= 3) {
        first3 = substr(s, 1, 3)
        last3  = substr(s, L - 2, 3)
      } else if (L > 0) {
        first3 = s
        last3  = s
      } else {
        first3 = "N/A"
        last3  = "N/A"
      }
      # Return elements separated by tabs
      return first3 "\t" last3
    }

    BEGIN { collecting = 0; seq = "" }

    /^>/ {
      if (collecting && length(seq) > 0) {
        print file_name, gff_tid, clean_and_extract(seq)
      }
      collecting = 0
      seq = ""
      header = substr($0, 2)
      split(header, a, " ")
      # Match exact species name or substring
      if (a[1] == reference_species || index(a[1], reference_species) > 0) {
        collecting = 1
      }
      next
    }

    {
      if (collecting) {
        gsub(/\r/, "", $0)
        seq = seq $0
      }
    }

    END {
      if (collecting && length(seq) > 0) {
        print file_name, gff_tid, clean_and_extract(seq)
      }
    }
    ' "$aln_file" >> "$summary"
done

echo "Done! Wrote results to $summary"