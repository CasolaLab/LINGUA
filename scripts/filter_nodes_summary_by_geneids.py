#!/usr/bin/env python3

"""
python3 filter_nodes_summary_by_geneids.py lsgs_to_remove.txt hyper_summary/nodes_summary_by_species

Filter <species>_nodes_summary.tsv files by removing geneIDs listed in an input file.

Usage:
python3 filter_nodes_summary_by_geneids.py \
    geneIDs_to_remove.txt \
    input_folder

Example:
python3 filter_nodes_summary_by_geneids.py bad_genes.txt nodes_summary_dir

Output:
For each file:
    <species>_nodes_summary.tsv
Produces:
    <species>_nodes_summary_filtered.tsv
"""

import sys
import csv
from pathlib import Path


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def load_geneids(path):
    geneids = set()
    with open(path) as f:
        for line in f:
            gid = line.strip()
            if gid:
                geneids.add(gid)
    return geneids


def normalize_geneid(gid):
    """Handle .p vs no .p consistency"""
    gid = gid.strip()
    if gid.endswith(".p"):
        return gid[:-2]
    return gid


def main():
    if len(sys.argv) != 3:
        die("Usage: python3 filter_nodes_summary_by_geneids.py geneIDs.txt input_folder")

    gene_file = Path(sys.argv[1])
    folder = Path(sys.argv[2])

    if not gene_file.exists():
        die(f"Gene list file not found: {gene_file}")
    if not folder.exists():
        die(f"Folder not found: {folder}")

    # load geneIDs (both exact + normalized)
    raw_ids = load_geneids(gene_file)
    remove_ids = set()
    for gid in raw_ids:
        remove_ids.add(gid)
        remove_ids.add(normalize_geneid(gid))

    files = sorted(folder.glob("*_nodes_summary.tsv"))

    if not files:
        die(f"No *_nodes_summary.tsv files found in {folder}")

    for infile in files:
        outfile = infile.with_name(infile.name.replace("_nodes_summary.tsv", "_nodes_summary_filtered.tsv"))

        with open(infile, newline="") as fin, open(outfile, "w", newline="") as fout:
            reader = csv.DictReader(fin, delimiter="\t")
            writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t")

            writer.writeheader()

            kept = 0
            removed = 0

            for row in reader:
                gid = row.get("geneID", "").strip()
                gid_norm = normalize_geneid(gid)

                if gid in remove_ids or gid_norm in remove_ids:
                    removed += 1
                    continue

                writer.writerow(row)
                kept += 1

        print(f"{infile.name}: kept={kept}, removed={removed}", file=sys.stderr)


if __name__ == "__main__":
    main()