#!/usr/bin/env python3
"""
==============================================================================
Script: step5_validate_reverse_complement.py
Description: Validates and corrects FASTA sequence orientation against a reference BED file using local dynamic programming alignments.
Usage: python3 step5_validate_reverse_complement.py <step4_output_dir> <mafextractor_dir> <bed_file> <reference_species> <threads>
==============================================================================
"""

import os
import sys
import concurrent.futures
from Bio import SeqIO, Align
from Bio.Seq import Seq

# --- Core Function --- 
def get_raw_maf_seq(maf_path, reference_species):
    """Extracts the raw, un-reversed sequence directly from the source MAF."""
    raw_seq = ""
    with open(maf_path, 'r') as f:
        for line in f:
            if line.startswith("s "):
                parts = line.split()
                if parts[1].startswith(reference_species):
                    raw_seq += parts[6].upper().replace("-", "")
    return raw_seq

def check_and_fix_orientation(aln_path, mafextractor_dir, bed_map, valid_chroms, reference_species):
    try:
        base_name = os.path.basename(aln_path).replace(".aln", "")
        maf_path = os.path.join(mafextractor_dir, base_name + ".txt")
        
        if not os.path.exists(maf_path):
            return "ERROR: RAW MAF FILE MISSING"
            
        # 1. Strip the reference species prefix if it exists in the filename
        base_name_no_sp = base_name
        if base_name.startswith(f"{reference_species}."):
            base_name_no_sp = base_name[len(f"{reference_species}."):]
            
        # 2. Find the exact chromosome match using the BED dictionary
        matched_chrom = None
        # Sort chromosomes by length (longest first) to prevent 'Chr1' from accidentally matching 'Chr10'
        for chrom in sorted(list(valid_chroms), key=len, reverse=True):
            if base_name_no_sp.startswith(f"{chrom}_"):
                matched_chrom = chrom
                break
                
        if not matched_chrom:
            return "ERROR: UNKNOWN CHROMOSOME FORMAT"
            
        # 3. Extract the exact coordinates 
        coord_str = base_name_no_sp[len(f"{matched_chrom}_"):]
        try:
            m_start, m_end = coord_str.split('_')
        except ValueError:
            return "ERROR: INVALID COORDINATE FORMAT"
            
        # 4. Lookup the intended strand in the BED map
        lookup_key = f"{matched_chrom}_{m_start}_{m_end}"
        target_strand = bed_map.get(lookup_key, "+") # Default to forward if missing from BED
        
        # 5. Get the original raw sequence from the MAF
        raw_seq = get_raw_maf_seq(maf_path, reference_species)
        if not raw_seq:
            return "ERROR: NO RAW REF SEQUENCE"

        # 6. Read the extracted ALN sequence
        records = list(SeqIO.parse(aln_path, "fasta"))
        if not records:
            return "NO_RECORDS"
            
        ref_rec = None
        for r in records:
            if r.id.lower().startswith(reference_species.lower()):
                ref_rec = r
                break
        if not ref_rec:
            return "NO_reference_species"

        aln_seq = str(ref_rec.seq).upper().replace("-", "")
        
        # 7. Check Current Orientation using Local DP Alignment
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 2.0
        aligner.mismatch_score = -3.0
        aligner.open_gap_score = -5.0
        aligner.extend_gap_score = -1.0
        
        score_fwd = aligner.score(raw_seq, aln_seq)
        score_rev = aligner.score(raw_seq, str(Seq(aln_seq).reverse_complement()))
        
        current_is_fwd = score_fwd >= score_rev
        
        # 8. Enforce BED Rules
        needs_flip = False
        if target_strand == "-" and current_is_fwd:
            needs_flip = True  # BED says -, but file is currently forward
        elif target_strand == "+" and not current_is_fwd:
            needs_flip = True  # BED says +, but file is currently reverse-complemented
            
        if needs_flip:
            # File orientation is wrong! Fix it.
            for i in range(len(records)):
                records[i].seq = records[i].seq.reverse_complement()
            # Overwrite the file with the corrected orientation
            SeqIO.write(records, aln_path, "fasta")
            return f"FIXED: FLIPPED TO MATCH BED ({target_strand})"
        else:
            return "PASS"
            
    except Exception as e:
        return f"ERROR: {str(e)}"

# === MAIN EXECUTION ===
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 reverse_complement_validation.py <step4_output_dir> <mafextractor_dir> <bed_file> <reference_species> <threads>")
        sys.exit(1)

    step4_output_dir = sys.argv[1]
    mafextractor_dir = sys.argv[2]
    bed_file = sys.argv[3]
    reference_species = sys.argv[4]
    
    try:
        threads = int(sys.argv[5])
    except ValueError:
        print("Error: threads must be an integer")
        sys.exit(1)
        
    print(f"Loading BED file to learn chromosome structure for {reference_species}...")
    bed_map = {}
    valid_chroms = set()
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split()
                # Ensure it's a valid BED line with at least 6 columns
                if len(parts) >= 6:
                    chrom = parts[0]
                    start = parts[1]
                    end = parts[2]
                    strand = parts[5]
                    valid_chroms.add(chrom)
                    bed_map[f"{chrom}_{start}_{end}"] = strand
    except Exception as e:
        print(f"Failed to load BED: {str(e)}")
        sys.exit(1)
        
    # Recursively search all subdirectories to find the .aln files
    aln_paths = []
    for root, dirs, files in os.walk(step4_output_dir):
        for f in files:
            if f.endswith(".aln"):
                aln_paths.append(os.path.join(root, f))
                
    total_files = len(aln_paths)
    print(f"Learned {len(valid_chroms)} unique chromosomes. Loaded {len(bed_map)} target coordinates.")
    print(f"Auditing orientation for {total_files} alignment files using {threads} threads...")
    
    completed = 0
    fixed_count = 0
    errors = 0
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(check_and_fix_orientation, path, mafextractor_dir, bed_map, valid_chroms, reference_species): path for path in aln_paths}
        
        for future in concurrent.futures.as_completed(futures):
            path = futures[future]
            base = os.path.basename(path).replace(".aln", "")
            status = future.result()
            completed += 1
            
            if "FIXED" in status:
                print(f"[{base}.aln]: {status}")
                fixed_count += 1
            elif "ERROR" in status or "NO_" in status:
                errors += 1
                
            if completed % 1000 == 0:
                print(f"Progress: {completed} / {total_files}")
                
    print("-" * 50)
    print("BED Orientation Audit Complete!")
    print(f"Total Files Checked: {total_files}")
    print(f"Files Fixed (Flipped to match BED): {fixed_count}")
    if errors > 0:
        print(f"Note: {errors} files could not be evaluated.")