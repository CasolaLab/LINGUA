#!/usr/bin/env bash
# ==============================================================================
# Script: step12_identify_start_codons.sh
# Description: Scans the stitched multi-species alignments to dynamically locate the biological start codon (ATG) in-frame.
# Usage: ./step12_identify_start_codons.sh <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv>
# ==============================================================================
set -euo pipefail

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv>"
  exit 1
fi

python3 - "$1" "$2" "$3" "$4" "$5" << 'EOF'
import os, sys, glob

input_dir, output_dir, bed_file, ref_species, summary_file = sys.argv[1:6]
os.makedirs(output_dir, exist_ok=True)
out_master = os.path.join(output_dir, "master_start.tsv")

bed_map = {}
with open(bed_file, "r") as f:
    for line in f:
        if line.strip() and not line.startswith("#"):
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                bed_map[parts[6].strip()] = {"chrom": parts[0], "strand": parts[5]}

# Load metadata strictly from the stitched summary file
# 1:tid, 2:chrom, 3:strand, 5:ref_length
meta_map = {}
with open(summary_file, "r") as f:
    next(f) # skip header
    for line in f:
        if line.strip():
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                meta_map[parts[0].strip()] = {
                    "tid": parts[1].strip(),
                    "chrom": parts[2].strip(),
                    "strand": parts[3].strip(),
                    "length": parts[5].strip()
                }

header = "file_name\tgff_transcriptID\tChromosome\tStart\tEnd\tSpecies\tStartCodonStart_RefPos\tStartCodonEnd_RefPos\tStartCodon\tgff_strand_ref\tcds_length\treference_length\n"
all_rows = []

for aln_file in glob.glob(os.path.join(input_dir, "*.aln")):
    file_name = os.path.basename(aln_file).replace(".aln", "")
    if not file_name.startswith(ref_species): continue
    
    # Pull attributes from summary map based on file_name
    metadata = meta_map.get(file_name, {})
    tid = metadata.get("tid", "NA")
    chrom = metadata.get("chrom", "NA")
    strand = metadata.get("strand", "NA")
    ref_length_val = metadata.get("length", "NA")

    # Keep original filename coordinate parsing
    parts = file_name.split("_")
    start_coord = parts[-2] if len(parts) > 1 else "0"
    end_coord = parts[-1] if len(parts) > 1 else "0"

    seqs, order = {}, []
    with open(aln_file, "r") as f:
        curr_sp, curr_seq = None, []
        for line in f:
            if line.startswith(">"):
                if curr_sp: seqs[curr_sp] = "".join(curr_seq)
                curr_sp = line[1:].strip().split(".")[0]
                order.append(curr_sp); curr_seq = []
            else: curr_seq.append(line.strip())
        if curr_sp: seqs[curr_sp] = "".join(curr_seq)

    ref_id = next((sp for sp in order if ref_species in sp or sp in tid or tid in sp), None)
    if not ref_id: continue
    
    ref_seq = seqs[ref_id].upper()
    ref_pos_map = []
    curr_p = 0
    for char in ref_seq:
        if char != "-": curr_p += 1
        ref_pos_map.append(curr_p)

    with open(os.path.join(output_dir, f"{file_name}_starts.tsv"), "w") as ind:
        ind.write(header)
        for sp in order:
            if sp == ref_id: continue
            sp_seq = seqs[sp].upper()
            sp_ungapped, sp_pos_map = "", []
            for i, char in enumerate(sp_seq):
                if char != "-":
                    sp_ungapped += char
                    sp_pos_map.append(i)
            
            # SCAN IN-FRAME ONLY (Steps of 3)
            for i in range(0, len(sp_ungapped) - 2, 3):
                codon = sp_ungapped[i:i+3]
                if codon == "ATG":
                    r_start, r_end = ref_pos_map[sp_pos_map[i]], ref_pos_map[sp_pos_map[i+2]]
                    row_data = [file_name, tid, chrom, start_coord, end_coord, sp, str(r_start), str(r_end), codon, strand, str(len(sp_ungapped)), str(ref_length_val)]
                    all_rows.append(row_data)
                    ind.write("\t".join(row_data) + "\n")

all_rows.sort(key=lambda x: (x[1], int(x[6]) if x[6].isdigit() else 0))
with open(out_master, "w") as f:
    f.write(header)
    for r in all_rows: f.write("\t".join(r) + "\n")

print(f"Start codons done. Master summary at {out_master}")
EOF