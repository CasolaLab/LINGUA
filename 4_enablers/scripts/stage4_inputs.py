#!/usr/bin/env python3
"""
stage4_inputs.py

One row per focal species of data/species.tsv with the stage 4 files that belong to it, for the per-species
steps of run_pipeline.sh (4b, 5b, 6b, 6c, 7c, 7d, 7e). Columns, tab-separated, no header:
  focal  query_prefix  out_prefix  nodes_prefix  aln_dir  indel_table  stop_table  ancseq_table
A file that is not there is NA; each step checks the columns it needs.

Files are found by content, not by name, because stage 3 names them per species in more than one way:
  aln_dir       the folder under RESULTS_DIR whose .aln files all start with "<focal>_"
  indel_table   the *_indels_summary.tsv in INDELS_DIR whose first row's file_name starts with "<focal>_"
  stop_table    the *_stop_codon_summary.tsv in STOPS_DIR found the same way
  ancseq_table  <out_prefix>_cdngs_tblastn_*_AncSeq.tsv in TBLASTN_DIR
nodes_prefix is the prefix of the step 7 files in nodes_summary_by_species/ (nodes_by_coordinates.py names them).

Usage: python3 stage4_inputs.py SPECIES_TSV RESULTS_DIR INDELS_DIR STOPS_DIR TBLASTN_DIR
"""
import csv
import glob
import os
import sys

import lingua_config as cfg
from nodes_by_coordinates import species_to_nodes_prefix

csv.field_size_limit(10 ** 9)


def first_file_name(path):
    with open(path, newline="") as fh:
        row = next(csv.DictReader(fh, delimiter="\t"), None)
    return (row or {}).get("file_name") or ""


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    species_tsv, results, indels, stops, tblastn = sys.argv[1:]
    with open(species_tsv, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for col in ("query_prefix", "out_prefix", "focal_species"):
        if rows and col not in rows[0]:
            sys.exit(f"ERROR: {species_tsv} has no {col} column")
    focal = [cfg.canonicalize(r["focal_species"]) for r in rows]

    aln = {}
    for d in sorted(glob.glob(os.path.join(results, "*"))):
        names = [os.path.basename(f) for f in glob.glob(os.path.join(d, "*.aln"))[:50]]
        hit = {sp for sp in focal for n in names if n.startswith(sp + "_")}
        if len(hit) == 1:
            sp = hit.pop()
            if sp in aln:
                sys.exit(f"ERROR: two alignment folders for {sp}: {aln[sp]} and {d}")
            aln[sp] = d

    def by_first_row(folder, pattern):
        out = {}
        for f in sorted(glob.glob(os.path.join(folder, pattern))):
            name = first_file_name(f)
            hit = [sp for sp in focal if name.startswith(sp + "_")]
            if len(hit) == 1:
                if hit[0] in out:
                    sys.exit(f"ERROR: two tables for {hit[0]}: {out[hit[0]]} and {f}")
                out[hit[0]] = f
        return out

    ind = by_first_row(indels, "*_indels_summary.tsv")
    stp = by_first_row(stops, "*_stop_codon_summary.tsv")
    for r, sp in zip(rows, focal):
        anc = sorted(glob.glob(os.path.join(tblastn, f"{r['out_prefix']}_cdngs_tblastn_*_AncSeq.tsv")))
        if len(anc) > 1:
            sys.exit(f"ERROR: more than one AncSeq table for {r['out_prefix']}: {anc}")
        print("\t".join([sp, r["query_prefix"], r["out_prefix"], species_to_nodes_prefix(sp),
                         aln.get(sp, "NA"), ind.get(sp, "NA"), stp.get(sp, "NA"), anc[0] if anc else "NA"]))


if __name__ == "__main__":
    main()
