#!/usr/bin/env python3
"""
==============================================================================
Script: step9_trimming_duplicate_bases.py
Description: Open Reading Frame Protector. Trims extraneous padding and intron bleed-over from stitched sequences by tightly anchoring the alignment against the true Reference FNA sequence. Guarantees the 3-bp codon frame is preserved.
Usage: python3 step9_trimming_duplicate_bases.py <step8_stitched_dir> <output_dir> <reference_nucleotide_fna_file> <bed_file> <threads>
==============================================================================
"""

import os
import sys
import concurrent.futures
from Bio import SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def get_projection_instructions(fna_seq, aln_ref_seq_gapped):
    """
    Aligns the untrimmed reference sequence against the biological truth (FNA).
    Generates a list of instructions to project the exact length/frame of the FNA 
    onto the multiple sequence alignment.
    """
    ungapped_ref = aln_ref_seq_gapped.replace("-", "").upper()
    
    if not ungapped_ref:
        return []
        
    # Map the ungapped sequence indices back to their original columns in the ALN
    ungapped_to_col = {}
    curr_u = 0
    for col_idx, char in enumerate(aln_ref_seq_gapped):
        if char != '-':
            ungapped_to_col[curr_u] = col_idx
            curr_u += 1
            
    # --- THE FRAME-SAFE SCORING MATRIX ---
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 5.0          
    aligner.mismatch_score = -1.0      
    aligner.open_gap_score = -5.0      
    aligner.extend_gap_score = -0.1    
    aligner.target_end_gap_score = 0.0 
    aligner.query_end_gap_score = 0.0  
    
    alignments = aligner.align(fna_seq, ungapped_ref)
    best_alignment = alignments[0]
    
    # Extract the gapped alignment strings (Biopython securely stores them as a tuple)
    fna_aln = best_alignment[0]
    ref_aln = best_alignment[1]
    
    instructions = []
    u_idx = 0 # Tracks our position in the ungapped_ref
    
    # Step through the alignment column by column to enforce the FNA frame
    for f_char, r_char in zip(fna_aln, ref_aln):
        if f_char != '-':
            # RULE A: The FNA expects a base here. We MUST output exactly one column.
            if r_char != '-':
                orig_col = ungapped_to_col[u_idx]
                instructions.append(('keep', orig_col))
                u_idx += 1
            else:
                instructions.append(('insert', None))
        else:
            # RULE B: The FNA does NOT have a base here (it's a MAF insertion artifact).
            if r_char != '-':
                u_idx += 1
                
    return instructions

def process_single_file(fname, input_dir, output_dir, truth_map, ref_species):
    try:
        aln_path = os.path.join(input_dir, fname)
        records = list(SeqIO.parse(aln_path, "fasta"))
        if not records:
            return "NO_RECORDS"
            
        ref_rec = None
        for r in records:
            if r.id.lower().startswith(ref_species.lower()):
                ref_rec = r
                break
        if not ref_rec:
            return f"NO_REF_SPECIES (Looking for: {ref_species})"
            
        fna_seq = None
        for key in truth_map:
            if key in fname:
                fna_seq = truth_map[key]
                break
                
        if not fna_seq:
            return "NO_TRUTH_MATCH"
            
        # Get the Frame-Safe mapping instructions
        instructions = get_projection_instructions(fna_seq, str(ref_rec.seq).upper())
        
        if not instructions:
            return "FAILED_ALIGNMENT"
            
        trimmed_records = []
        for rec in records:
            new_seq_chars = []
            is_ref = rec.id.lower().startswith(ref_species.lower())
            original_seq = str(rec.seq)
            
            # Apply the instructions to every species in the ALN file
            for action, val in instructions:
                if action == 'keep':
                    new_seq_chars.append(original_seq[val])
                elif action == 'insert':
                    if is_ref:
                        new_seq_chars.append('N') 
                    else:
                        new_seq_chars.append('-') 
                        
            trimmed_seq = "".join(new_seq_chars)
            
            # Only keep species that actually have data remaining after trimming
            if len(trimmed_seq.replace("-", "").replace("N", "")) > 0:
                trimmed_records.append(SeqRecord(Seq(trimmed_seq), id=rec.id, description=""))
                
        # --- BUILT IN VALIDATION ---
        status_msg = "SUCCESS"
        if trimmed_records:
            final_ref_seq = ""
            for r in trimmed_records:
                if r.id.lower().startswith(ref_species.lower()):
                    final_ref_seq = str(r.seq).replace("-", "").upper() 
                    break
                    
            if len(final_ref_seq) != len(fna_seq):
                status_msg = f"ERROR: Trimming failed. Output ({len(final_ref_seq)}) != FNA ({len(fna_seq)})"
            elif final_ref_seq != fna_seq:
                mismatches = sum(1 for a, b in zip(final_ref_seq, fna_seq) if a != b)
                if 'N' in final_ref_seq:
                    status_msg = f"SUCCESS (Frame Preserved: Injected Ns. {mismatches} chars differ)"
                else:
                    status_msg = f"SUCCESS (Frame Preserved: {mismatches} ambiguous bases absorbed)"

            out_path = os.path.join(output_dir, fname)
            SeqIO.write(trimmed_records, out_path, "fasta")
            return status_msg
            
    except Exception as e:
        return f"ERROR: {str(e)}"

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 step9_trimming_duplicate_bases.py <step8_stitched_dir> <output_dir> <reference_nucleotide_fna_file> <bed_file> <threads>")
        sys.exit(1)

    # RE-ORDERED INPUTS TO MATCH YOUR COMMAND
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    reference_nucleotide_fna_file = sys.argv[3]
    bed_file = sys.argv[4]
    
    # Extract the species prefix from the BED file name (e.g., "Homo_sapiens.bed" -> "Homo_sapiens")
    ref_species = os.path.basename(bed_file).replace(".bed", "")
    
    try:
        threads = int(sys.argv[5])
    except ValueError:
        print("Error: threads must be an integer")
        sys.exit(1)
        
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Targeting reference species: {ref_species}")
    print("Loading Truth Sequences from FNA...")
    truth_map = {}
    try:
        for rec in SeqIO.parse(reference_nucleotide_fna_file, "fasta"):
            clean_id = rec.id.split()[0].strip()
            truth_map[clean_id] = str(rec.seq).upper()
    except Exception as e:
        print(f"Failed to load FNA file: {str(e)}")
        sys.exit(1)
        
    aln_files = [f for f in os.listdir(input_dir) if f.endswith(".aln")]
    total_files = len(aln_files)
    print(f"Loaded {len(truth_map)} relevant truth sequences.")
    print(f"Processing {total_files} alignment files using {threads} threads...")
    
    completed = 0
    warnings = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(process_single_file, fname, input_dir, output_dir, truth_map, ref_species): fname for fname in aln_files}
        
        for future in concurrent.futures.as_completed(futures):
            fname = futures[future]
            status = future.result()
            completed += 1
            
            # FIXED: Catch and print all skipped files and silent errors
            if status != "SUCCESS":
                print(f"[{fname}]: {status}")
                if "ERROR" in status or "NO_" in status or "FAILED_" in status:
                    warnings += 1
                
            if completed % 1000 == 0:
                print(f"Progress: {completed} / {total_files}")
                
    print("-" * 50)
    print(f"Trimming complete! {total_files} files processed.")
    if warnings > 0:
        print(f"Note: {warnings} files encountered warnings, errors, or were skipped.")