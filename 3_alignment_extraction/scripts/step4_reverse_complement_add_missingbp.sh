#!/usr/bin/env bash
# ==============================================================================
# Script: step4_reverse_complement_add_missingbp.sh
# Description: Converts MAF to ALN (FASTA). Concatenates MAF blocks per species, pads species absent from a block with dashes, and ensures correct BED orientation.
# Usage: ./step4_reverse_complement_add_missingbp.sh -m <step3_filtered_dir> -o <output_dir> -b <bed_file> [-t <threads>]
#
# Bases a species skips between consecutive MAF blocks are unaligned sequence and
# are deliberately NOT added back: filling them would force gap columns into the
# reference and shift the frame of the other species relative to it.
# ==============================================================================

# Default values
THREADS=4

# Usage information
usage() {
    echo "Usage: $0 -m <step3_filtered_dir> -o <output_dir> -b <bed_file> [-t <threads>]"
    echo "  -m : Directory containing MAF .txt files"
    echo "  -o : Output directory for .aln files"
    echo "  -b : reference BED file (for strand lookup)"
    echo "  -t : Number of parallel threads (default: 4)"
    exit 1
}

# -----------------------------------------------------------------------------
# WORKER FUNCTION: Block Concatenation
# -----------------------------------------------------------------------------
process_file() {
    local maf_file="$1"
    local output_dir="$2"
    local bed_file="$3"

    maf_name="${maf_file##*/}"
    base="${maf_name%.txt}"
    file_outdir="$output_dir/$base"
    mkdir -p "$file_outdir"
    fna_out="$file_outdir/${base}.aln"
    > "$fna_out"

    # 1. Get Strand from BED
    coord_part=$(echo "$base" | grep -oE '[0-9]+_[0-9]+$')
    m_start="${coord_part%_*}"; m_end="${coord_part#*_}"
    ref_strand=$(grep -m 1 "_${m_start}_${m_end}" "$bed_file" | cut -f6 || echo "+")

    # 2. Embedded Python logic to concatenate blocks
    python3 -c '
import sys

maf_file = sys.argv[1]
out_file = sys.argv[2]
strand = sys.argv[3]

# --- 1. PARSE MAF BLOCKS ---
sp_order = []
with open(maf_file) as f:
    for line in f:
        if line.startswith("s "):
            sp_full = line.split()[1]
            sp = sp_full.split(".")[0] if "." in sp_full else sp_full
            if sp not in sp_order: sp_order.append(sp)

blocks = []
current_block = {}
with open(maf_file) as f:
    for line in f:
        # Matches all versions of mafExtractor headers
        if line.startswith("a ") or line.startswith("ascore=") or line.startswith("a\t") or line.strip() == "a":
            if current_block:
                blocks.append(current_block)
                current_block = {}
        elif line.startswith("s "):
            parts = line.split()
            sp_full = parts[1]
            sp = sp_full.split(".")[0] if "." in sp_full else sp_full
            chr_name = sp_full.split(".", 1)[1] if "." in sp_full else "unknown"
            seq = parts[6].upper()
            current_block[sp] = {"chr": chr_name, "seq": seq}

    if current_block: blocks.append(current_block)

final_seqs = {sp: "" for sp in sp_order}
chr_maps = {}

# --- 2. CONCATENATE BLOCKS ---
for b in blocks:
    block_len = len(list(b.values())[0]["seq"])
    for sp in sp_order:
        if sp in b:
            chr_maps[sp] = b[sp]["chr"]
            final_seqs[sp] += b[sp]["seq"]
        else:
            final_seqs[sp] += "-" * block_len

# Reverse complement if needed
if strand == "-":
    comp = str.maketrans("ATCG", "TAGC")
    for sp in sp_order:
        final_seqs[sp] = final_seqs[sp][::-1].translate(comp)

# Write output
with open(out_file, "w") as out:
    for sp in sp_order:
        chr_val = chr_maps.get(sp, "unknown")
        out.write(f">{sp}.{chr_val}\n")
        seq = final_seqs[sp]
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")
' "$maf_file" "$fna_out" "$ref_strand"
}

export -f process_file

# -----------------------------------------------------------------------------
# MAIN ORCHESTRATOR
# -----------------------------------------------------------------------------

# Parse command line arguments
while getopts "m:o:b:t:" opt; do
    case ${opt} in
        m ) step3_filtered_dir=$OPTARG ;;
        o ) OUTPUT_DIR=$OPTARG ;;
        b ) BED_FILE=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        * ) usage ;;
    esac
done

if [[ -z "$step3_filtered_dir" || -z "$OUTPUT_DIR" || -z "$BED_FILE" ]]; then
    usage
fi

echo "--- Starting MAF Block Concatenation Pipeline ---"
echo "Threads: $THREADS"
echo "Input:   $step3_filtered_dir"

find "$step3_filtered_dir" -name "*.txt" | parallel --progress -j "$THREADS" process_file {} "$OUTPUT_DIR" "$BED_FILE"

echo "--- Generating Summary File ---"
SUMMARY_FILE="$OUTPUT_DIR/extraction_summary.tsv"
echo -e "file_name\ttranscript_id\tchromosome\tstrand\tlength\tnum_species\tstart\tend" > "$SUMMARY_FILE"

find "$OUTPUT_DIR" -name "*.aln" | while read -r aln_path; do
    base_name=$(basename "$aln_path" .aln)
    transcript_id=$(echo "$base_name" | awk -F'_' '{print $(NF-3)}')
    chrom=$(echo "$base_name" | awk -F'_' '{print $(NF-2)}')
    start=$(echo "$base_name" | awk -F'_' '{print $(NF-1)}')
    end=$(echo "$base_name" | awk -F'_' '{print $NF}')
    seq_len=$(grep -v ">" "$aln_path" | tr -d '\n' | wc -c)
    num_sp=$(grep -c ">" "$aln_path")
    strand=$(grep -m 1 "_${start}_${end}" "$BED_FILE" | cut -f6 || echo "+")
    echo -e "${base_name}\t${transcript_id}\t${chrom}\t${strand}\t${seq_len}\t${num_sp}\t${start}\t${end}" >> "$SUMMARY_FILE"
done

echo "Summary saved to: $SUMMARY_FILE"
echo "--- All files processed successfully ---"
