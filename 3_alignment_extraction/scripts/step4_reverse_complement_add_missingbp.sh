#!/usr/bin/env bash
# ==============================================================================
# Script: step4_reverse_complement_add_missingbp.sh
# Description: Converts MAF to ALN (FASTA). Analyzes missing BP CSV to dynamically inject genome sequence or padding for authorized gaps, and ensures correct BED orientation.
# Usage: ./step4_reverse_complement_add_missingbp.sh -m <step3_filtered_dir> -g <genome_fasta_dir> -o <output_dir> -b <bed_file> -c <step2_combined_summary_file> [-t <threads>]
# ==============================================================================

# Default values
THREADS=4

# Usage information
usage() {
    echo "Usage: $0 -m <step3_filtered_dir> -g <genome_fasta_dir> -o <output_dir> -b <bed_file> -c <missing_bp_csv> [-t <threads>]"
    echo "  -m : Directory containing MAF .txt files"
    echo "  -g : Directory containing genome .fasta files"
    echo "  -o : Output directory for .aln files"
    echo "  -b : reference BED file (for strand lookup)"
    echo "  -c : Combined missing bases summary CSV file (to validate gap extraction)"
    echo "  -t : Number of parallel threads (default: 4)"
    exit 1
}

# -----------------------------------------------------------------------------
# WORKER FUNCTION: CSV-Guarded Block-Anchored Gap Injection
# -----------------------------------------------------------------------------
process_file() {
    local maf_file="$1"
    local fasta_genome_dir="$2"
    local output_dir="$3"
    local bed_file="$4"
    local step2_combined_summary_file="$5"

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

    # 2. Embedded Python logic to dynamically inject gaps based on CSV
    python3 -c '
import sys, os, subprocess, csv

maf_file = sys.argv[1]
genome_fasta_dir = sys.argv[2]
out_file = sys.argv[3]
strand = sys.argv[4]
csv_path = sys.argv[5]

base_name = os.path.basename(maf_file)

# --- 1. PARSE CSV WHITELIST ---
# Only species listed in the CSV for THIS specific file are allowed to have gaps filled
authorized_species = set()
try:
    with open(csv_path, "r") as cf:
        # Check if comma or tab separated
        sample_line = cf.readline()
        delim = "," if "," in sample_line else "\t"
        cf.seek(0)
        
        reader = csv.reader(cf, delimiter=delim)
        for row in reader:
            if len(row) >= 2:
                csv_fname = os.path.basename(row[0])
                csv_sp = row[1]
                if csv_fname == base_name:
                    authorized_species.add(csv_sp)
except Exception as e:
    pass

# --- 2. PARSE MAF BLOCKS ---
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
            start = int(parts[2])
            size = int(parts[3])
            seq = parts[6].upper()
            current_block[sp] = {"chr": chr_name, "start": start, "end": start + size, "seq": seq}
            
    if current_block: blocks.append(current_block)

final_seqs = {sp: "" for sp in sp_order}
last_ends = {sp: -1 for sp in sp_order}
chr_maps = {}

# --- 3. CONCATENATE BLOCKS & DYNAMICALLY INJECT GAPS ---
for b in blocks:
    block_len = len(list(b.values())[0]["seq"])
    
    inter_gaps = {}
    max_gap = 0
    
    for sp in sp_order:
        if sp in b:
            chr_maps[sp] = b[sp]["chr"]
            if last_ends[sp] != -1:
                gap = b[sp]["start"] - last_ends[sp]
                
                # CSV GUARD: Only acknowledge the gap if the species is in the whitelist!
                # This prevents huge introns from being pulled into the alignment.
                if gap > 0 and sp in authorized_species:
                    inter_gaps[sp] = gap
                    if gap > max_gap: 
                        max_gap = gap
                else:
                    inter_gaps[sp] = 0
            else:
                inter_gaps[sp] = 0
        else:
            inter_gaps[sp] = 0
            
    # Inject padding/samtools seqs for the entire matrix
    if max_gap > 0:
        for sp in sp_order:
            # If this species has an authorized gap, fill it!
            if inter_gaps[sp] > 0:
                s_start = last_ends[sp] + 1
                s_end = b[sp]["start"]
                g_files = [f for f in os.listdir(genome_fasta_dir) if f.startswith(sp)]
                inter_seq = ""
                if g_files:
                    g_path = os.path.join(genome_fasta_dir, g_files[0])
                    cmd = ["samtools", "faidx", g_path, f"{chr_maps[sp]}:{s_start}-{s_end}"]
                    try:
                        res = subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode()
                        inter_seq = "".join(res.strip().split("\n")[1:]).upper()
                    except: pass
                
                pad_len = max(0, max_gap - len(inter_seq))
                final_seqs[sp] += inter_seq + ("-" * pad_len)
            
            # If this species does NOT have a gap, insert dashes to keep alignment perfect
            else:
                final_seqs[sp] += "-" * max_gap

    # Add the aligned block
    for sp in sp_order:
        if sp in b:
            final_seqs[sp] += b[sp]["seq"]
            last_ends[sp] = b[sp]["end"]
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
' "$maf_file" "$fasta_genome_dir" "$fna_out" "$ref_strand" "$step2_combined_summary_file"
}

export -f process_file

# -----------------------------------------------------------------------------
# MAIN ORCHESTRATOR
# -----------------------------------------------------------------------------

# Parse command line arguments
while getopts "m:g:o:b:c:t:" opt; do
    case ${opt} in
        m ) step3_filtered_dir=$OPTARG ;;
        g ) GENOME_FASTA_DIR=$OPTARG ;;
        o ) OUTPUT_DIR=$OPTARG ;;
        b ) BED_FILE=$OPTARG ;;
        c ) step2_combined_summary_file=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        * ) usage ;;
    esac
done

if [[ -z "$step3_filtered_dir" || -z "$GENOME_FASTA_DIR" || -z "$OUTPUT_DIR" || -z "$BED_FILE" || -z "$step2_combined_summary_file" ]]; then
    usage
fi

echo "--- Starting CSV-Guarded Extraction Pipeline ---"
echo "Threads: $THREADS"
echo "Input:   $step3_filtered_dir"
echo "CSV Map: $step2_combined_summary_file"

find "$step3_filtered_dir" -name "*.txt" | parallel --progress -j "$THREADS" process_file {} "$GENOME_FASTA_DIR" "$OUTPUT_DIR" "$BED_FILE" "$step2_combined_summary_file"

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