#!/usr/bin/env python3
"""
==============================================================================
Script: step17_combine_all_species_info.py
Description: Master Coordinate Aggregator. Compiles the mapped genomic coordinates of the functional exons for ALL species relative to the reference sequence. Outputs the final coordinate-space dictionary for downstream genomic plotting and analysis.
Usage: python3 step17_combine_all_species_info.py <step9_trimmed_dir> <step3_filtered_dir> <bed_file> <reference_species> <output.tsv>
==============================================================================
"""

import os
import sys
import glob

def main():
    if len(sys.argv) != 6:
        print("Usage: python3 step17_combine_all_species_info.py <step9_trimmed_dir> <step3_filtered_dir> <bed_file> <reference_species> <output.tsv>")
        sys.exit(1)

    step9_trimmed_dir = sys.argv[1]
    step3_filtered_dir = sys.argv[2]
    bed_file = sys.argv[3]
    reference_species = sys.argv[4]
    out_tsv = sys.argv[5]

    print("1. Parsing Reference BED file to map exons...")
    genes = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): 
                continue
            cols = line.strip().split('\t')
            if len(cols) >= 7:
                chrom = cols[0].strip()
                start = int(cols[1].strip())
                end = int(cols[2].strip())
                
                # Clean the Transcript ID (e.g., Ah119U00020.1__CDS_... -> Ah119U00020.1)
                tid_raw = cols[6].strip()
                tid = tid_raw.split("__")[0]
                
                if tid not in genes:
                    genes[tid] = {'chrom': chrom, 'exons': []}
                genes[tid]['exons'].append({'start': start, 'end': end})

    # Precompute the overall start/end for the reference species
    for tid, data in genes.items():
        starts = [e['start'] for e in data['exons']]
        ends = [e['end'] for e in data['exons']]
        data['min_s'] = min(starts)
        data['max_e'] = max(ends)

    out_rows = []
    out_rows.append(["file_name", "species_name", "chromosome", "start", "end"])
    
    aln_files = glob.glob(os.path.join(step9_trimmed_dir, "*.aln"))
    print(f"2. Found {len(aln_files)} final alignment files. Scanning for multi-species coordinates...")

    for count, aln_path in enumerate(aln_files, 1):
        if count % 1000 == 0:
            print(f"   Processed {count} / {len(aln_files)} files...")
            
        fname = os.path.basename(aln_path)
        clean_fname = fname.replace(".aln", "")
        
        # Strip off the reference species prefix if it exists
        if clean_fname.startswith(reference_species + "_"):
            clean_fname = clean_fname[len(reference_species)+1:]
            
        # --- BULLETPROOF ID EXTRACTION ---
        # Instead of guessing where the chromosome starts by counting underscores,
        # we check the filename against our verified list of BED transcript IDs.
        matched_tid = None
        for tid in genes.keys():
            if clean_fname.startswith(tid + "_"):
                matched_tid = tid
                break
                    
        if not matched_tid:
            print(f"Warning: Could not link {fname} to BED file records. Skipping.")
            continue
            
        gene_data = genes[matched_tid]
        ref_chrom = gene_data['chrom']
        ref_overall_start = gene_data['min_s']
        ref_overall_end = gene_data['max_e']

        # Determine which species actually survived in this final ALN file
        species_in_aln = set()
        with open(aln_path, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    sp = line[1:].strip().split('.')[0]
                    species_in_aln.add(sp)

        # Initialize coordinate tracking for non-reference species
        sp_coords = {sp: {'chroms': set(), 'starts': [], 'ends': []} for sp in species_in_aln if sp != reference_species}

        # Hunt down the raw MAF .txt files for every exon in this gene
        for exon in gene_data['exons']:
            maf_fname = f"{reference_species}.{ref_chrom}_{exon['start']}_{exon['end']}.txt"
            maf_path = os.path.join(step3_filtered_dir, maf_fname)
                    
            if not os.path.exists(maf_path):
                continue # If the exon was completely dropped during extraction, skip it
                    
            with open(maf_path, 'r') as f:
                for line in f:
                    if line.startswith('s '):
                        pts = line.strip().split()
                        if len(pts) >= 7:
                            sp_chrom = pts[1]
                            
                            # Parse species.chromosome safely
                            if '.' in sp_chrom:
                                sp_name, chrom = sp_chrom.split('.', 1)
                            else:
                                sp_name, chrom = sp_chrom, "NA"
                            
                            # If this is one of our target species, extract true coordinates
                            if sp_name in sp_coords:
                                try:
                                    m_start = int(pts[2])
                                    m_length = int(pts[3])
                                    m_strand = pts[4]
                                    m_srcSize = int(pts[5])
                                    
                                    # Convert MAF negative strand backwards-coordinates into standard forward genomic coordinates
                                    if m_strand == '+':
                                        true_s = m_start
                                        true_e = m_start + m_length
                                    else:
                                        true_s = m_srcSize - m_start - m_length
                                        true_e = m_srcSize - m_start
                                        
                                    sp_coords[sp_name]['chroms'].add(chrom)
                                    sp_coords[sp_name]['starts'].append(true_s)
                                    sp_coords[sp_name]['ends'].append(true_e)
                                except ValueError:
                                    pass

        # Generate the final master row entries for every species in this file
        for sp in sorted(list(species_in_aln)):
            if sp == reference_species:
                # We already know the reference species info perfectly from the BED mapping
                out_rows.append([fname, sp, ref_chrom, str(ref_overall_start), str(ref_overall_end)])
            else:
                # Aggregate the non-reference spanning coordinates across all valid exons
                if sp_coords[sp]['starts']:
                    # Join chroms by comma if a species bizarrely jumped scaffolds mid-gene
                    chroms_sorted = sorted(list(sp_coords[sp]['chroms']))
                    final_chrom = ",".join(chroms_sorted) if chroms_sorted else "NA"
                    
                    # The gene span is from the start of the first mapped exon to the end of the last mapped exon
                    final_start = str(min(sp_coords[sp]['starts']))
                    final_end = str(max(sp_coords[sp]['ends']))
                    
                    out_rows.append([fname, sp, final_chrom, final_start, final_end])
                else:
                    # Fallback if the species was in the ALN but missing coordinate data
                    out_rows.append([fname, sp, "NA", "NA", "NA"])

    print("3. Writing Master Coordinates Summary...")
    with open(out_tsv, 'w') as f:
        for row in out_rows:
            f.write("\t".join(row) + "\n")

    print("-" * 50)
    print(f"Success! Processed {len(aln_files)} files and logged {len(out_rows) - 1} species coordinate entries.")
    print(f"Saved to: {out_tsv}")

if __name__ == "__main__":
    main()