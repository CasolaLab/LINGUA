#!/usr/bin/env bash
# ==============================================================================
# Script: step14_identify_indels.sh
# Description: Biological mutation identifier. Maps insertions and deletions across the alignments. Specifically flags catastrophic 3-bp insertions that themselves form a premature stop codon (TAG, TAA, TGA).
# Usage: ./step14_identify_indels.sh <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv>
# ==============================================================================

set -euo pipefail

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv>"
  exit 1
fi

python3 - "$1" "$2" "$3" "$4" "$5" << 'EOF'
import os, sys, glob

input_dir, output_dir, bed_file, reference_species, summary_file = sys.argv[1:6]
os.makedirs(output_dir, exist_ok=True)
out_master = os.path.join(output_dir, "master_indels.tsv")

STOP_CODONS = {"TAG", "TGA", "TAA"}

# We keep bed_map for legacy compatibility but prioritize summary_file data
bed_map = {}
with open(bed_file, "r") as f:
    for line in f:
        if line.strip() and not line.startswith("#"):
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                bed_map[parts[6].strip()] = {"chrom": parts[0], "strand": parts[5]}

# Load metadata strictly from the stitched summary file
# 1:tid, 2:chrom, 3:strand, 4:ref_length
meta_map = {}
with open(summary_file, "r") as f:
    next(f) # skip header
    for line in f:
        if line.strip():
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                meta_map[parts[0].strip()] = {
                    "tid": parts[1].strip(),
                    "chrom": parts[2].strip(),
                    "strand": parts[3].strip(),
                    "length": parts[4].strip()
                }

header = "file_name\tgff_transcriptID\tChromosome\tStart\tEnd\tSpecies\tIndelStart_RefPos\tIndelEnd_RefPos\tIndelLength\tgff_strand_ref\tcds_length\treference_length\n"
all_rows = []

for aln_file in glob.glob(os.path.join(input_dir, "*.aln")):
    file_name = os.path.basename(aln_file).replace(".aln", "")
    if not file_name.startswith(reference_species): continue
    
    # Pull attributes from summary map 
    metadata = meta_map.get(file_name, {})
    tid = metadata.get("tid", "NA")
    chrom = metadata.get("chrom", "NA")
    strand = metadata.get("strand", "NA")
    ref_length_val = metadata.get("length", "NA")
    
    # --- ROBUST ID EXTRACTION (For Coordinates) ---
    clean_fname = file_name
    if clean_fname.startswith(reference_species + "_"):
        clean_fname = clean_fname[len(reference_species)+1:]
        
    parts = clean_fname.split("_")
    start_coord = parts[-2] if len(parts) > 1 else "0"
    end_coord = parts[-1] if len(parts) > 1 else "0"

    species_seqs = {}
    with open(aln_file, "r") as f:
        curr_sp, curr_seq = None, []
        for line in f:
            if line.startswith(">"):
                if curr_sp: species_seqs[curr_sp] = "".join(curr_seq)
                curr_sp = line[1:].strip().split(".")[0]
                curr_seq = []
            else:
                curr_seq.append(line.strip())
        if curr_sp: species_seqs[curr_sp] = "".join(curr_seq)

    ref_seq = species_seqs.get(reference_species)
    if not ref_seq: continue
    
    cds_len = str(len(ref_seq.replace("-", "")))
    
    ref_pos_map = []
    curr_p = 0
    for char in ref_seq:
        if char != "-": curr_p += 1
        ref_pos_map.append(curr_p)

    file_specific_rows = []
    
    for sp, sp_seq in species_seqs.items():
        if sp == reference_species: continue
        in_indel, indel_type, start_col = False, None, 0
        
        for i in range(len(ref_seq)):
            r, s = ref_seq[i], sp_seq[i]
            cur = "Deletion" if (s == "-" and r != "-") else "Insertion" if (s != "-" and r == "-") else None
            
            if cur != indel_type:
                if in_indel:
                    length = (i - 1) - start_col + 1
                    is_insertion = (indel_type == "Insertion")
                    indel_seq = sp_seq[start_col:i].replace("-", "").upper()
                    
                    # --- STRICT BIOLOGICAL FILTER ---
                    keep = False
                    if length % 3 != 0: 
                        keep = True  # Frameshift
                    elif length == 3 and is_insertion and indel_seq in STOP_CODONS: 
                        keep = True  # 3bp Premature Stop
                        
                    if keep:
                        row_data = [file_name, tid, chrom, start_coord, end_coord, sp, str(ref_pos_map[start_col]), str(ref_pos_map[i-1]), str(length), strand, cds_len, str(ref_length_val)]
                        all_rows.append(row_data)
                        file_specific_rows.append(row_data)
                        
                if cur: 
                    in_indel, indel_type, start_col = True, cur, i
                else: 
                    in_indel, indel_type = False, None
        
        if in_indel:
            length = (len(ref_seq) - 1) - start_col + 1
            is_insertion = (indel_type == "Insertion")
            indel_seq = sp_seq[start_col:].replace("-", "").upper()
            
            # --- STRICT BIOLOGICAL FILTER (End of sequence) ---
            keep = False
            if length % 3 != 0: 
                keep = True
            elif length == 3 and is_insertion and indel_seq in STOP_CODONS: 
                keep = True
                
            if keep:
                row_data = [file_name, tid, chrom, start_coord, end_coord, sp, str(ref_pos_map[start_col]), str(ref_pos_map[-1]), str(length), strand, cds_len, str(ref_length_val)]
                all_rows.append(row_data)
                file_specific_rows.append(row_data)

    indiv_path = os.path.join(output_dir, f"{file_name}_indels.tsv")
    with open(indiv_path, "w") as ind:
        ind.write(header)
        for r in file_specific_rows: ind.write("\t".join(r) + "\n")

all_rows.sort(key=lambda x: (x[1], int(x[6]) if x[6].isdigit() else 0))
with open(out_master, "w") as f:
    f.write(header)
    for r in all_rows: f.write("\t".join(r) + "\n")

print(f"Indels done. Master summary at {out_master}")
EOF