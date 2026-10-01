#!/usr/bin/env python3
"""
==============================================================================
Script: step16_extract_start_and_stops.py
Description: Fast directory-wide boundary parser. Pulls the first and last 3 bases of sequences, strips gaps, and validates if they match known canonical translation start/stop motifs.
Usage: python3 step16_extract_start_and_stops.py <step9_trimmed_dir> <output_tsv>
==============================================================================
"""

import os
import sys

STOP_CODONS = {"tag", "tga", "taa"}
# Added .aln and .fasta_aln to supported extensions
FASTA_EXTS  = {".fa", ".fasta", ".fna", ".aln", ".fasta_aln"}

def parse_first_sequence(filepath):
    """Return (header, sequence) for the first FASTA record in a file."""
    header = None
    seq_parts = []
    try:
        with open(filepath, "r") as fh:
            for line in fh:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        break
                    header = line[1:]
                else:
                    if header is not None:
                        seq_parts.append(line)
    except Exception as e:
        print(f"[error] could not read {filepath}: {e}", file=sys.stderr)
        return None, None
        
    if header is None:
        return None, None
    return header, "".join(seq_parts)

def process_dir(step9_trimmed_dir, out_fh):
    total = 0
    for dirpath, dirnames, filenames in os.walk(step9_trimmed_dir):
        dirnames.sort()   # Reproducible order
        for fname in sorted(filenames):
            ext = os.path.splitext(fname)[1].lower()
            if ext not in FASTA_EXTS:
                continue
            
            filepath = os.path.join(dirpath, fname)
            header, seq = parse_first_sequence(filepath)
            
            if header is None:
                print(f"[skip] {filepath}: no valid sequence", file=sys.stderr)
                continue

            # Clean sequence: lowercase and remove gaps/Ns
            seq_clean = seq.lower().replace("-", "").replace("n", "")
            
            if len(seq_clean) < 3:
                print(f"[skip] {filepath}: sequence too short after cleaning", file=sys.stderr)
                continue

            first3    = seq_clean[:3]
            last3     = seq_clean[-3:]
            no_start  = 0 if first3 == "atg" else 1
            no_stop   = 0 if last3 in STOP_CODONS else 1

            # Extract just the filename from the path
            just_filename = os.path.basename(filepath)

            out_fh.write(f"{just_filename}\t{header}\t{first3}\t{last3}\t{no_start}\t{no_stop}\n")
            
            total += 1
            if total % 10000 == 0:
                print(f"[progress] {total} files processed", file=sys.stderr)

    print(f"[done] {total} files processed", file=sys.stderr)

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 extract_starts_stops.py <step9_trimmed_dir> <output_tsv>", file=sys.stderr)
        sys.exit(1)

    step9_trimmed_dir  = sys.argv[1]
    output_tsv = sys.argv[2]

    if not os.path.isdir(step9_trimmed_dir):
        print(f"[error] not a directory: {step9_trimmed_dir}", file=sys.stderr)
        sys.exit(1)

    with open(output_tsv, "w") as out_fh:
        # Header for the TSV
        out_fh.write("Filename\tgeneID\t1stCodon\tLastCodon\tNoSTART\tNoSTOP\n")
        process_dir(step9_trimmed_dir, out_fh)

if __name__ == "__main__":
    main()