#!/usr/bin/env bash
# ==============================================================================
# Script: step2_missingbp_summaries.sh
# Description: Maps species presence/absence and calculates missing base pairs between MAF alignment blocks to generate summary CSVs.
# Usage: ./step2_missingbp_summaries.sh <mafextractor_dir> <output_dir>
# ==============================================================================

set -euo pipefail

mafextractor_dir="$1"
output_dir="$2"

if [[ -z "${mafextractor_dir:-}" || -z "${output_dir:-}" ]]; then
  echo "Usage: $0 /path/to/mafextractor_dir /path/to/output_dir"
  exit 1
fi

mkdir -p "$output_dir"

combined_all_csv="$output_dir/combined_all_species.csv"
orf_summary_csv="$output_dir/orf_summary.csv"
species_summary_csv="$output_dir/species_summary.csv"

# Headers
echo "file_name,non_ref_ID,species_name,chr_scaff,prev_end,next_start,missing_bases,block_number,species_absent" > "$combined_all_csv"
echo "file_name,start_coord,end_coord,num_blocks,num_species_present" > "$orf_summary_csv"
echo "file_name,non_ref_ID,species_name,chr_scaff,missing_bases,species_absent,chrscaff" > "$species_summary_csv"

# Ensure the for-glob doesn’t iterate a literal pattern if no files match
shopt -s nullglob

# --- Execution Loop ---
for maf_file in "$mafextractor_dir"/*.txt; do
  [[ -e "$maf_file" ]] || continue
  base_name=$(basename "$maf_file" .txt)
  file_name=$(basename "$maf_file")
  out_csv="$output_dir/${base_name}_gaps.csv"
  echo "Processing $maf_file -> $out_csv"
  echo "file_name,non_ref_ID,species_name,chr_scaff,prev_end,next_start,missing_bases,block_number,species_absent" > "$out_csv"

  awk -v comb_all="$combined_all_csv" -v orf_sum="$orf_summary_csv" \
      -v species_sum="$species_summary_csv" -v out_csv="$out_csv" -v fname="$file_name" '
  BEGIN {
    OFS=","
    block=0
    total_blocks=0
  }
  {
    sub(/^[ \t]+/, "", $0)
  }
  /^a/ {
    block++
    total_blocks=block
    next
  }
  /^s/ {
    full_id=$2
    split(full_id, arr, ".")
    spname=arr[1]
    chr=(length(arr)>1 && arr[2] != "" ? arr[2] : 0)  # 0 if missing

    start=$3+0
    len=$4+0
    end=start+len

    species_seen[spname]=1
    block_species[block,spname]=full_id
    block_end[block,spname]=end
    block_start[block,spname]=start

    # Record canonical IDs on first sighting
    if (!(spname in species_full_id)) {
      species_full_id[spname]=full_id
      species_chr_first[spname]=chr
    }

    # Compute missing bases for this present occurrence
    missing=0
    if (spname in prev_end) {
      gap=start - prev_end[spname]
      if (gap>0) missing=gap
    }
    # Count blocks where missing>0 (not sum of lengths)
    if (missing > 0) species_missing_block_count[spname]++

    # Track scaffold/chromosome switch across present blocks
    if (spname in last_present_chr) {
      if (chr != last_present_chr[spname]) chr_switch_count[spname]++
    }
    last_present_chr[spname]=chr

    # Write present row
    print fname, full_id, spname, chr, (spname in prev_end ? prev_end[spname] : 0), start, missing, block, "No" >> comb_all
    print fname, full_id, spname, chr, (spname in prev_end ? prev_end[spname] : 0), start, missing, block, "No" >> out_csv

    prev_end[spname]=end
    next
  }
  END {
    # Species list
    spcount=0
    for (s in species_seen) sp_list[++spcount]=s

    # Absent rows and count absences
    for (i=1; i<=total_blocks; i++) {
      for (si=1; si<=spcount; si++) {
        sp=sp_list[si]
        key=i SUBSEP sp
        if (!(key in block_species)) {
          full_id_for_absent = (sp in species_full_id ? species_full_id[sp] : sp)
          print fname, full_id_for_absent, sp, 0, 0, 0, 0, i, "Yes" >> comb_all
          print fname, full_id_for_absent, sp, 0, 0, 0, 0, i, "Yes" >> out_csv
          species_absent_count[sp]++
        }
      }
    }

    # ORF summary
    n=split(fname, parts, "_")
    startc=parts[n-1]
    endc=parts[n]
    sub(/\.txt$/, "", endc)

    num_species_present=0
    for (si=1; si<=spcount; si++) {
      sp=sp_list[si]
      present_all=1
      for (i=1; i<=total_blocks; i++) {
        if (!(i SUBSEP sp in block_species)) { present_all=0; break }
      }
      if (present_all==1) num_species_present++
    }
    print fname, startc, endc, total_blocks, num_species_present >> orf_sum

    # Species summary per file
    for (si=1; si<=spcount; si++) {
      sp=sp_list[si]
      nonref = (sp in species_full_id ? species_full_id[sp] : sp)
      chrfirst = (sp in species_chr_first ? species_chr_first[sp] : 0)
      miss_blocks = (sp in species_missing_block_count ? species_missing_block_count[sp] : 0)
      absent_ct = (sp in species_absent_count ? species_absent_count[sp] : 0)
      switches = (sp in chr_switch_count ? chr_switch_count[sp] : 0)
      print fname, nonref, sp, chrfirst, miss_blocks, absent_ct, switches >> species_sum
    }
  }
  ' "$maf_file"
done

# Sorting outputs
tmp_sorted="$combined_all_csv.sorted"
(
  head -n 1 "$combined_all_csv"
  tail -n +2 "$combined_all_csv" | sort -t, -k1,1V -k5,5n -k6,6n
) > "$tmp_sorted"
mv "$tmp_sorted" "$combined_all_csv"

tmp_sorted2="$orf_summary_csv.sorted"
(
  head -n 1 "$orf_summary_csv"
  tail -n +2 "$orf_summary_csv" | sort -t, -k1,1V -k2,2n -k3,3n
) > "$tmp_sorted2"
mv "$tmp_sorted2" "$orf_summary_csv"

tmp_sorted3="$species_summary_csv.sorted"
(
  head -n 1 "$species_summary_csv"
  tail -n +2 "$species_summary_csv" | sort -t, -k1,1V -k3,3
) > "$tmp_sorted3"
mv "$tmp_sorted3" "$species_summary_csv"

echo "All done!"
echo "- Per-file CSVs saved in $output_dir"
echo "- Combined all-species CSV (sorted): $combined_all_csv"
echo "- ORF summary CSV (sorted): $orf_summary_csv"
echo "- Species summary CSV (sorted): $species_summary_csv"