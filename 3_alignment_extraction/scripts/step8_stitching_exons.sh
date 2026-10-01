#!/usr/bin/env bash
# ================================================================================
# SCRIPT: step8_stitching_exons.sh
# ================================================================================
# Description: 
#   Stitches individual exon alignment blocks into continuous transcript models
#   based on reference BED coordinates.
#
# Usage:
# export TREE_FILE=/path/to/tree.nwk
# export REF_SPECIES=reference_species_name
# export SPECIES_NAME="reference_species_name"

#   ./step8_stitching_exons.sh <step4_output_dir> <bed_file> <output_dir> <reference_species>
# ================================================================================

set -euo pipefail
LC_ALL=C

# USAGE:
# export TREE_FILE=/path/to/tree.nwk
# ./stitching_exons.sh <step4_output_dir> <bed_file> <output_dir> <reference_species>

if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <step4_output_dir> <bed_file> <output_dir> <reference_species>"
  exit 1
fi

mafft_root="$1"
combined_tsv="$2"
out_root="$3"
REF_SPECIES="$4"
SPECIES_NAME="$4"
SPECIES_NAME_CLEAN="${SPECIES_NAME// /_}"

mkdir -p "$out_root"

# --- CHANGE 1: TSV Header with Start/End Columns ---
master_summary="$out_root/master_summary.tsv"
echo -e "file_name\tgff_transcriptID\tchromosome\tstrand\tlength\tnumber_of_exons\taligned_species\tdistant_species\tdistant_node\tstart\tend" > "$master_summary"

ms_body_tmp="$out_root/.master_summary.body.tmp"
: > "$ms_body_tmp"

: "${TREE_FILE:?Please set TREE_FILE to a Newick file, e.g., export TREE_FILE=/path/to/tree.nwk}"

# [Helpers]
sanitize_tid() { local t="$1"; t="$(tr -d '\r\n' <<< "$t")"; t="$(sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//' <<< "$t")"; printf "%s" "$t"; }
sanitize_chr() { local c="$1"; c="$(tr -d '\r\n' <<< "$c")"; c="$(sed -E 's/^[Cc]hr/chr/' <<< "$c" | tr -d '[:space:]')"; while [[ "$c" == "-"* || "$c" == "_"* || "$c" == "."* ]]; do c="${c#-}"; c="${c#_}"; c="${c#.}"; done; printf "%s" "$c"; }
chr_to_num() { local c="$1"; c="$(sed -E 's/^[Cc]hr([0-9]+).*/\1/;t; s/.*/0/' <<< "$c")"; printf "%s" "$c"; }
clean_csv() { tr -d '\r\n' <<< "$1"; }
repeat_char() { local ch="$1" n="$2"; (( n > 0 )) && printf "%*s" "$n" "" | tr ' ' "$ch" || printf ""; }

# --- CRITICAL FIX: Match by Coordinates, Verify by Chromosome ---
declare -A COORD_DICT
while IFS= read -r -d '' p; do
  base="$(basename "$p" .aln)"
  
  # Read coordinates strictly from the right side of the filename
  start=$(echo "$base" | awk -F'_' '{print $(NF-1)}')
  end=$(echo "$base" | awk -F'_' '{print $NF}')
  
  # Store mapped by start and end (append | in case of multiple files sharing coordinates)
  COORD_DICT["${start}:${end}"]+="$p|"
done < <(find "$mafft_root" -type f -name '*.aln' -print0)

# Grouping logic
declare -A T_EXONS T_STRAND T_CHR VALIDATED_ALN_PATH
{
  read -r _header || true
  while IFS=$'\t' read -r chrom chromStart chromEnd name score strand proteinID || [[ -n "${chrom:-}" ]]; do
    # Strip invisible carriage returns from Windows/Mac files to prevent empty outputs
    chrom="$(echo "$chrom" | tr -d '\r\n ')"
    chromStart="$(echo "$chromStart" | tr -d '\r\n ')"
    chromEnd="$(echo "$chromEnd" | tr -d '\r\n ')"
    proteinID="$(echo "$proteinID" | tr -d '\r\n ')"
    
    [[ -z "${chrom:-}" || -z "${chromStart:-}" || -z "${proteinID:-}" ]] && continue
    
    # 1. Match by coordinates
    paths_str="${COORD_DICT[${chromStart}:${chromEnd}]:-}"
    [[ -z "$paths_str" ]] && continue
    
    # 2. Verify Chromosome ID is in the filename
    chrom_lower=$(echo "$chrom" | tr '[:upper:]' '[:lower:]')
    
    IFS='|' read -r -a path_array <<< "$paths_str"
    aln_path=""
    for p in "${path_array[@]}"; do
      [[ -z "$p" ]] && continue
      base="$(basename "$p" .aln)"
      base_lower=$(echo "$base" | tr '[:upper:]' '[:lower:]')
      
      if [[ "$base_lower" == *"${chrom_lower}"* ]]; then
        aln_path="$p"
        break
      fi
    done
    
    [[ -z "$aln_path" ]] && continue
    
    # Build the exact key expected by your original downstream logic
    key="${chrom_lower}:${chromStart}:${chromEnd}"
    VALIDATED_ALN_PATH["$key"]="$aln_path"
    
    T_EXONS["$proteinID"]+="${chromStart}:${key}"$'\n'
    T_CHR["$proteinID"]="$chrom"
    T_STRAND["$proteinID"]="$strand"
  done
} < "$combined_tsv"

# Tree distances 
dist_tsv="$out_root/_tree_distances.tsv"
python3 - << 'PY' > "${dist_tsv}"
from Bio import Phylo
import itertools, os
tree = Phylo.read(os.environ['TREE_FILE'], "newick")
taxa = [t.name for t in tree.get_terminals() if t.name]
for a,b in itertools.combinations(taxa, 2):
    try: d = tree.distance(a,b)
    except: d = None
    if d is not None: print(f"{a}\t{b}\t{float(d)}")
PY

declare -A DIST
while IFS=$'\t' read -r a b d; do DIST["$a|$b"]="$d"; DIST["$b|$a"]="$d"; done < "$dist_tsv"

mrca_info_for_set() {
  python3 - "$TREE_FILE" "$@" << 'PY'
import sys
from Bio import Phylo
t = Phylo.read(sys.argv[1], "newick")
tips = {x.name for x in t.get_terminals() if x.name}
present = [s for s in sys.argv[2:] if s in tips]
if not present: print("NA"); sys.exit(0)
mrca = t.common_ancestor(present)
print(getattr(mrca, "name", f"MRCA_unlabeled;clade_size={len(present)}"))
PY
}

# Process each transcript
for tid in "${!T_EXONS[@]}"; do
  exons="${T_EXONS[$tid]}"; chr="${T_CHR[$tid]}"; strand="${T_STRAND[$tid]}"
  if [[ "$strand" == "-" ]]; then mapfile -t ordered < <(printf "%s" "$exons" | sort -t: -k1,1nr)
  else mapfile -t ordered < <(printf "%s" "$exons" | sort -t: -k1,1n); fi

  exon_keys=(); for rec in "${ordered[@]}"; do [[ -n "$rec" ]] && exon_keys+=("${rec#*:}"); done
  exon_count="${#exon_keys[@]}"; tid_clean="$(sanitize_tid "$tid")"; chr_clean="$(sanitize_chr "$chr")"
  
  # Determine Coords
  min_start=""; max_end=""; total_len=0
  for key in "${exon_keys[@]}"; do
    IFS=':' read -r _c _s _e <<< "$key"
    [[ -z "$min_start" || $_s -lt $min_start ]] && min_start="$_s"
    [[ -z "$max_end" || $_e -gt $max_end ]] && max_end="$_e"
    (( total_len += (_e - _s + 1) ))
  done

  out_base="${SPECIES_NAME_CLEAN}_${tid_clean}_${chr_clean}_${min_start}_${max_end}"
  out_aln="$out_root/${out_base}.aln"

  # [Stitching Logic]
  declare -A stitched seen_species
  species_order=(); stitched_len=0; first_exon_done=0

  for key in "${exon_keys[@]}"; do
    aln_path="${VALIDATED_ALN_PATH[$key]:-}"
    declare -A exon_seqs
    exon_len=0; exon_order=()
    while IFS=$'\t' read -r hdr seq; do
      sp="${hdr#>}"; sp="${sp%%[[:space:]]*}"
      exon_seqs["$sp"]="$seq"; exon_order+=("$sp")
      l=${#seq}; (( l > exon_len )) && exon_len=$l
    done < <(awk 'BEGIN{hdr=""; seq=""} /^>/ { if(seq){print hdr "\t" seq} hdr=$0; seq=""; next } { gsub("\r",""); seq=seq $0 } END { if(seq) print hdr "\t" seq }' "$aln_path")
    
    for sp in "${species_order[@]}"; do [[ -z "${exon_seqs[$sp]:-}" ]] && stitched["$sp"]+="$(repeat_char "-" "$exon_len")"; done
    if (( first_exon_done == 0 )); then
      for sp in "${exon_order[@]}"; do
        stitched["$sp"]="${exon_seqs[$sp]}"; seen_species["$sp"]=1; species_order+=("$sp")
      done
      first_exon_done=1
    else
      for sp in "${exon_order[@]}"; do
        if [[ -z "${seen_species[$sp]:-}" ]]; then
          stitched["$sp"]="$(repeat_char "-" "$stitched_len")${exon_seqs[$sp]}"; seen_species["$sp"]=1; species_order+=("$sp")
        else stitched["$sp"]+="${exon_seqs[$sp]}"; fi
      done
    fi
    (( stitched_len += exon_len ))
    unset exon_seqs exon_order
  done

  # Write ALN
  { for sp in "${species_order[@]}"; do echo ">$sp"
    seq="${stitched[$sp]}"; i=0
    while (( i < ${#seq} )); do echo "${seq:i:60}"; (( i += 60 )); done
  done; } > "$out_aln"

  # Stats and Phylo
  species_count="${#species_order[@]}"
  mapfile -t sp_list < <(printf "%s\n" "${species_order[@]}" | sed 's/\..*$//')
  
  # [Distant Species / MRCA Logic]
  distant="NA"
  if [[ -n "$REF_SPECIES" ]] && printf '%s\n' "${sp_list[@]}" | grep -qx "$REF_SPECIES"; then
      maxd="-1"; far="NA"
      for s in "${sp_list[@]}"; do
        [[ "$s" == "$REF_SPECIES" ]] && continue
        d="${DIST["$REF_SPECIES|$s"]:-}"; [[ -z "$d" ]] && continue
        awk -v d="$d" -v m="$maxd" 'BEGIN{exit !(d+0>m+0)}' && { maxd="$d"; far="$s"; }
      done
      distant="$far"
  else
      maxd="-1"; pair="NA"
      for ((i=0;i<${#sp_list[@]};i++)); do
        for ((j=i+1;j<${#sp_list[@]};j++)); do
          a="${sp_list[i]}"; b="${sp_list[j]}"; d="${DIST["$a|$b"]:-}"; [[ -z "$d" ]] && continue
          awk -v d="$d" -v m="$maxd" 'BEGIN{exit !(d+0>m+0)}' && { maxd="$d"; pair="$a|$b"; }
        done
      done
      distant="$pair"
  fi
  dnode="$(mrca_info_for_set "${sp_list[@]}")"

  # --- CHANGE 2: Print as TSV with Start and End appended ---
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(clean_csv "$out_base")" \
    "$(sanitize_tid "$tid")" \
    "$(sanitize_chr "$chr")" \
    "$(clean_csv "$strand")" \
    "$total_len" \
    "$exon_count" \
    "$species_count" \
    "$(clean_csv "$distant")" \
    "$(clean_csv "$dnode")" \
    "$min_start" \
    "$max_end" >> "$ms_body_tmp"

  unset stitched seen_species species_order
done

# Sort and finish
tail -n +1 "$ms_body_tmp" | sort -t $'\t' -k2,2V -k3,3 -k10,10n >> "$master_summary"
rm -f "$dist_tsv" "$ms_body_tmp"
echo "Done. Summary at $master_summary"