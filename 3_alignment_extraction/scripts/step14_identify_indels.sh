#!/usr/bin/env bash
# ==============================================================================
# Script: step14_identify_indels.sh
# Description: Biological mutation identifier. Maps insertions and deletions across the alignments. Specifically flags catastrophic 3-bp insertions that themselves form a premature stop codon (TAG, TAA, TGA).
# Usage: ./step14_identify_indels.sh <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv> [missing_data_mask_tsv]
#
# 2026-10-28 (LINGUA stage 3, NEE revision): optional 6th argument, the missing-data mask written by step 8c
# (step8c_missing_data.tsv: alignment, species, exon, ref_start, ref_end, reason). A deletion of a species whose
# whole span lies inside a masked reference range of that species is NOT reported: the species has no aligned
# sequence there (exon absent from the alignment, or exon piece removed by step 8c), which is missing data, not a
# deletion. Masked reference ranges (BED genomic coordinates, half-open) are converted to CDS positions through the
# exon order of the transcript. Without the argument the behaviour is unchanged.
# 2026-10-29 (his decision): a deletion that overlaps the mask is TRIMMED to its part(s) outside the mask: the masked
# exon is missing data, the bases missing beyond it in a neighbouring exon are a deletion. Each remainder is tested as
# a frameshift (length not a multiple of 3) on its own; remainders of 1 bp at a mask boundary are dropped (boundary
# ambiguity, the former 1-bp tolerance). A deletion with no overlap is unchanged.
# ==============================================================================

set -euo pipefail

if [ "$#" -lt 5 ] || [ "$#" -gt 6 ]; then
  echo "Usage: $0 <step9_trimmed_dir> <output_dir> <bed_file> <reference_species> <step11_summary_tsv> [missing_data_mask_tsv]"
  exit 1
fi

python3 - "$1" "$2" "$3" "$4" "$5" "${6:-}" << 'EOF'
import os, sys, glob

input_dir, output_dir, bed_file, reference_species, summary_file, mask_file = sys.argv[1:7]
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

# Load metadata strictly from the step 11 summary file
# 1:tid, 2:chrom, 3:strand; reference length = actual_nt (trimmed reference length),
# found by name because step 11 appends its columns after however many step 8 wrote
meta_map = {}
with open(summary_file, "r") as f:
    summary_header = f.readline().rstrip("\r\n").split("\t")
    if "actual_nt" not in summary_header:
        sys.exit(f"Error: no actual_nt column in {summary_file}; pass the step 11 summary")
    len_idx = summary_header.index("actual_nt")
    for line in f:
        if line.strip():
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) > len_idx:
                meta_map[parts[0].strip()] = {
                    "tid": parts[1].strip(),
                    "chrom": parts[2].strip(),
                    "strand": parts[3].strip(),
                    "length": parts[len_idx].strip()
                }

# missing-data mask -> per (alignment, species) list of CDS intervals [a, b] (1-based, inclusive)
exons_by_tid = {}
with open(bed_file, "r") as f:
    for line in f:
        if line.strip() and not line.startswith("#"):
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                exons_by_tid.setdefault(parts[6].strip(), []).append((int(parts[1]), int(parts[2]), parts[5].strip()))
def cds_interval(tid, r0, r1):
    """CDS positions covered by the genomic half-open range [r0, r1) of transcript tid, clipped to its exons."""
    exs = exons_by_tid.get(tid)
    if not exs:
        return None
    minus = exs[0][2] == "-"
    order = sorted(exs, key=lambda x: x[0], reverse=minus)
    cum, out = 0, []
    for s0, e0, _ in order:
        lo, hi = max(r0, s0), min(r1, e0)
        if lo < hi:
            if minus:
                a, b = cum + (e0 - hi) + 1, cum + (e0 - lo)
            else:
                a, b = cum + (lo - s0) + 1, cum + (hi - s0)
            out.append((a, b))
        cum += e0 - s0
    return out
masked = {}
n_mask_rows = 0
if mask_file:
    n_unmapped = 0
    with open(mask_file, "r") as f:
        next(f)
        for line in f:
            aln, sp, exon, r0, r1, reason = line.rstrip("\n").split("\t")
            # transcript of the alignment from the step 11 summary (contig names may contain underscores)
            tid = meta_map.get(aln, {}).get("tid")
            if tid is None:
                sys.exit(f"ERROR: mask alignment {aln} not in the step 11 summary")
            ivs = cds_interval(tid, int(r0), int(r1))
            if not ivs:
                # a block made only of the padding base step 4 adds beyond an exon end (trimmed by step 9) maps
                # to nothing and is ignored; anything else outside the exons is an error
                exs = exons_by_tid.get(tid, [])
                padding = int(r1) - int(r0) <= 1 and any(int(r0) in (s0 - 1, e0) for s0, e0, _ in exs)
                if padding:
                    continue
                n_unmapped += 1
            for a, b in ivs or []:
                masked.setdefault((aln, sp), []).append((a, b))
                n_mask_rows += 1
    if n_unmapped:
        sys.exit(f"ERROR: {n_unmapped} mask rows did not map to any exon of their transcript")
    # coalesce overlapping or adjacent masked intervals: a species absent from two consecutive exons shows one gap
    # run spanning both, and that run must count as masked
    for key, ivs in masked.items():
        ivs.sort()
        merged = []
        for a, b in ivs:
            if merged and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        masked[key] = merged
def unmasked_parts(aln, sp, a, b):
    """[a, b] minus the masked intervals of (aln, sp); remainders of <= 1 bp next to a mask boundary are dropped."""
    ivs = masked.get((aln, sp), [])
    if not any(x <= b and y >= a for x, y in ivs):
        return [(a, b)], False
    parts, cur = [], a
    for x, y in ivs:
        if y < cur or x > b:
            continue
        if x > cur:
            parts.append((cur, x - 1))
        cur = max(cur, y + 1)
    if cur <= b:
        parts.append((cur, b))
    return [(x, y) for x, y in parts if y - x + 1 > 1], True
n_masked_out = 0
n_trimmed = 0
n_trimmed_rows = 0
n_remainder_inframe = 0
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

    def indel_rows(is_insertion, length, indel_seq, sp, a, b):
        """Rows to report for one indel run: the run itself (insertion, or deletion outside the mask) or the
        unmasked remainders of a deletion overlapping the mask, each tested as a frameshift."""
        global n_masked_out, n_trimmed, n_trimmed_rows, n_remainder_inframe
        def row(x, y, L):
            return [file_name, tid, chrom, start_coord, end_coord, sp, str(x), str(y), str(L), strand, cds_len, str(ref_length_val)]
        if is_insertion:
            if length % 3 != 0 or (length == 3 and indel_seq in STOP_CODONS):
                return [row(a, b, length)]
            return []
        parts, overlapped = unmasked_parts(file_name, sp, a, b)
        if not overlapped:
            return [row(a, b, length)] if length % 3 != 0 else []
        if not parts:
            n_masked_out += 1
            return []
        out = []
        for x, y in parts:
            if (y - x + 1) % 3 != 0:
                out.append(row(x, y, y - x + 1))
            else:
                n_remainder_inframe += 1
        n_trimmed += 1
        n_trimmed_rows += len(out)
        return out

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
    
    ref_pos_map = []
    curr_p = 0
    for char in ref_seq:
        if char != "-": curr_p += 1
        ref_pos_map.append(curr_p)

    file_specific_rows = []
    
    for sp, sp_seq in species_seqs.items():
        if sp == reference_species: continue
        cds_len = str(len(sp_seq.replace("-", "")))
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
                    for row_data in indel_rows(is_insertion, length, indel_seq, sp, ref_pos_map[start_col], ref_pos_map[i-1]):
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
            for row_data in indel_rows(is_insertion, length, indel_seq, sp, ref_pos_map[start_col], ref_pos_map[-1]):
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
if mask_file:
    print(f"Missing-data mask: {n_mask_rows} masked CDS intervals; {n_masked_out} deletions inside them not reported; "
          f"{n_trimmed} deletions overlapping them trimmed to {n_trimmed_rows} unmasked frameshift remainder(s) "
          f"({n_remainder_inframe} in-frame remainders dropped)")
EOF