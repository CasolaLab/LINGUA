#!/usr/bin/env python3
"""
Usage:
python3 rename_file_getIDs.py clsgs_45k_1line 45k_renamed brassicaceae_geneIDs_species_names.tsv

python3 rename_file_getIDs.py input_dir output_dir geneIDs_species_names.tsv

"""


import os
import sys


def parse_fasta_headers(filepath):
    """
    Parse FASTA headers and return a list of (geneID, species_name).
    Header format expected:
    >geneID__SpeciesName
    Example:
    >Aa_G1000330.h1.t1__Aalpina
    """
    entries = []

    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]

                if "__" not in header:
                    print(f"WARNING: header without '__' separator in file {filepath}: {header}", file=sys.stderr)
                    continue

                gene_id, species_name = header.split("__", 1)
                entries.append((gene_id, species_name))

    return entries


def copy_file(infile, outfile):
    with open(infile, "r") as fin, open(outfile, "w") as fout:
        for line in fin:
            fout.write(line)


def main():
    if len(sys.argv) != 4:
        print("Usage: python3 rename_file_getIDs.py input_dir output_dir geneIDs_species_names.tsv")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    tsv_path = sys.argv[3]

    if not os.path.isdir(input_dir):
        print(f"ERROR: input_dir does not exist: {input_dir}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    fasta_exts = {".fa", ".faa", ".fasta", ".fas"}

    all_rows = []

    for filename in sorted(os.listdir(input_dir)):
        infile = os.path.join(input_dir, filename)

        if not os.path.isfile(infile):
            continue

        ext = os.path.splitext(filename)[1].lower()
        if ext not in fasta_exts:
            continue

        entries = parse_fasta_headers(infile)

        if not entries:
            print(f"WARNING: no valid FASTA headers found in {filename}", file=sys.stderr)
            continue

        seq_count = len(entries)

        # Use species name from headers
        species_names = {species for _, species in entries}
        if len(species_names) != 1:
            print(f"WARNING: multiple species names found in {filename}: {species_names}", file=sys.stderr)

        # Use the first species name
        species_name = entries[0][1]
        new_filename = f"{species_name.lower()}_cdngs_{seq_count}.faa"
        outfile = os.path.join(output_dir, new_filename)

        copy_file(infile, outfile)

        for gene_id, species in entries:
            all_rows.append((new_filename, gene_id, species))

        print(f"Processed: {filename} -> {new_filename}")

    with open(tsv_path, "w") as outtsv:
        outtsv.write("New_Filename\tgeneID\tSpecies_name\n")
        for new_filename, gene_id, species in all_rows:
            outtsv.write(f"{new_filename}\t{gene_id}\t{species}\n")

    print(f"\nDone.")
    print(f"Renamed FASTA files written to: {output_dir}")
    print(f"TSV table written to: {tsv_path}")


if __name__ == "__main__":
    main()