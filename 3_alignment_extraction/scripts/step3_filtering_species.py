#!/usr/bin/env python3
"""
==============================================================================
Script: step3_filtering_species.py
Description: Universal Reference-Protected MAF Filter. Removes species from alignment blocks if they jump chromosomes or jump further than a specified threshold (default 50kb).
Usage: python3 step3_filtering_species.py -i <mafextractor_dir> -o <output_dir> -r <reference_species> -s <output_summary_csv_name> [-j 50000] [-t 8]
==============================================================================
"""

import os
import sys
import argparse
from multiprocessing import Pool

# === ARGUMENT PARSING ===
def parse_args():
    parser = argparse.ArgumentParser(description="Universal Reference-Protected MAF Filter")
    parser.add_argument("-i", "--mafextractor_dir", required=True, help="mafextractor_dir directory containing raw MAF files")
    parser.add_argument("-o", "--output_dir", required=True, help="output_dir directory for filtered MAF files")
    parser.add_argument("-r", "--ref", required=True, help="Reference species name (e.g., Arabidopsis_thaliana)")
    parser.add_argument("-s", "--summary", required=True, help="Path for the output_dir purge summary CSV")
    parser.add_argument("-j", "--jump", type=int, default=50000, help="Jump threshold in bp (default: 50000)")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads to use")
    return parser.parse_args()

# --- Core Function --- 
def process_maf(args_tuple):
    in_path, out_dir, ref_name, jump_thresh = args_tuple
    base_name = os.path.basename(in_path)
    out_path = os.path.join(out_dir, base_name)
    
    to_purge = set()
    
    try:
        with open(in_path, 'r') as f:
            lines = f.readlines()
            
        blocks = []
        current_block = []
        
        # Parse all blocks into memory
        for line in lines:
            is_header = line.startswith('a ') or line.startswith('a\t') or line.startswith('ascore=') or line.strip() == 'a'
            
            if is_header:
                if current_block:
                    blocks.append(current_block)
                current_block = [line]
            elif line.strip():
                if current_block:
                    current_block.append(line)
                    
        # Flush the final block
        if current_block:
            blocks.append(current_block)
            
        # ==========================================
        # PASS 1: AUDIT FOR CHROMOSOME & 50KB JUMPS
        # ==========================================
        last_coords = {}
        for block in blocks:
            for line in block[1:]:
                if line.startswith('s '):
                    parts = line.split()
                    sp_full = parts[1]
                    sp_name = sp_full.split('.')[0] if '.' in sp_full else sp_full
                    chr_name = sp_full.split('.', 1)[1] if '.' in sp_full else "unknown"
                    
                    start = int(parts[2])
                    size = int(parts[3])
                    end = start + size # Calculate the exact end coordinate
                    
                    if sp_name in last_coords:
                        last_chr, last_end = last_coords[sp_name]
                        
                        if sp_name != ref_name:
                            # 1. Did it jump to a completely different chromosome?
                            if last_chr != chr_name:
                                to_purge.add(sp_name)
                            # 2. Did it stay on the same chromosome but jump too far?
                            else:
                                jump = abs(start - last_end)
                                if jump > jump_thresh:
                                    to_purge.add(sp_name)
                                    
                    # Save the end coordinate for the next loop
                    last_coords[sp_name] = (chr_name, end)

        # ==========================================
        # PASS 2: WRITE CLEANSED BLOCKS
        # ==========================================
        final_blocks = []
        for block in blocks:
            filtered_lines = []
            has_ref_seq = False
            
            header = block[0].strip()
            filtered_lines.append(header)
            
            for line in block[1:]:
                if line.startswith('s '):
                    parts = line.split()
                    sp_full = parts[1]
                    sp_name = sp_full.split('.')[0] if '.' in sp_full else sp_full
                    
                    if sp_name == ref_name:
                        filtered_lines.append(line.strip())
                        has_ref_seq = True
                    elif sp_name not in to_purge:
                        filtered_lines.append(line.strip())
            
            # Only write the block if the Reference Sequence is still in it
            if has_ref_seq:
                final_blocks.append("\n".join(filtered_lines) + "\n\n")
                
        # Write to disk
        with open(out_path, 'w') as f:
            f.write("##maf version=1\n\n")
            for fb in final_blocks:
                f.write(fb)
        
        # Log results
        if to_purge:
            return [f"{base_name},{p},purged_due_to_jump" for p in to_purge]
        else:
            return [f"{base_name},None,clean"]
            
    except Exception as e:
        return [f"{base_name},ERROR,{str(e)}"]

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    files = [f for f in os.listdir(args.mafextractor_dir) if f.lower().endswith(('.maf', '.txt'))]
    task_args = [(os.path.join(args.mafextractor_dir, f), args.output_dir, args.ref, args.jump) for f in files]

    print(f"Processing {len(files)} files using {args.threads} threads...")
    with Pool(args.threads) as pool:
        all_results = pool.map(process_maf, task_args)
    
    with open(args.summary, 'w') as out:
        out.write("file_name,species_purged,status\n")
        for res_list in all_results:
            for res in res_list:
                out.write(res + "\n")
                
    print(f"Done! Summary saved to {args.summary}")

if __name__ == "__main__":
    main()