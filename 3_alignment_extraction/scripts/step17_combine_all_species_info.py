#!/usr/bin/env python3
"""
==============================================================================
Script: step17_combine_all_species_info.py
Description: Master Coordinate Aggregator. Compiles the mapped genomic coordinates of the functional exons for ALL
species relative to the reference sequence. Outputs the final coordinate-space dictionary for downstream genomic
plotting and analysis.

Usage: python3 step17_combine_all_species_info.py <step9_trimmed_dir> <step3_filtered_dir> <bed_file>
                                                  <reference_species> <output.tsv> [--mask <step8c_missing_data.tsv>]
                                                  [--tree <newick>]

For every non-reference species (living or ancestral) in each step 9 alignment, the coordinates are the span of the
step 3 MAF blocks that supplied the bases the species kept in that alignment:
  * only blocks on the species' kept contig, read from its header in the step 9 alignment (species.contig);
  * with --mask (step8c/step8c_missing_data.tsv), not blocks whose reference range lies inside a range masked for
    that species (exons step 8c gapped as other contig, opposite strand, beyond the span limit or not collinear);
  * not blocks that cover only the extra base step 1 extracts past the exon end (no overlap with the BED exon).
Before 2026-11-10 every block of every exon was used, so a species whose exons step 8 / 8c had gapped on another
contig got a comma-joined chromosome list and a span reaching across contigs.

With --mask, the run stops if a block outside the mask lies on a contig other than the species' kept contig (step 8c
and step 9 disagree). Without --mask, blocks gapped by step 8c on the kept contig cannot be recognized, and a species
with several rows (step 8c not run) keeps the old comma-joined chromosome list of its rows.

Output columns: file_name, species_name, chromosome, start, end, species_type. start/end are 0-based half-open
plus-strand coordinates. species_type is living (a tip of --tree, default $TREE_FILE), ancestor (a labeled internal
node of the tree, whatever its naming scheme), not_in_tree, or NA without a tree. The reference row is the BED span.
==============================================================================
"""

import argparse
import collections
import glob
import os
import sys


def merge(ivs):
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def inside(ivs, a, b):
    return any(x <= a and b <= y for x, y in ivs)


def tree_node_types(path):
    if not path:
        return None
    from Bio import Phylo
    tree = Phylo.read(path, "newick")
    types = {c.name: "ancestor" for c in tree.get_nonterminals() if c.name}
    types.update({c.name: "living" for c in tree.get_terminals() if c.name})
    if not any(v == "ancestor" for v in types.values()):
        print(f"Warning: {path} has no internal node labels; ancestors will be typed not_in_tree", file=sys.stderr)
    return types


def read_maf_blocks(path, ref):
    """[(ref_start, ref_end, [(species, contig, start, end), ...])] with species coordinates on the plus strand."""
    blocks, cur, ref_range = [], [], None

    def flush():
        nonlocal cur, ref_range
        if cur or ref_range:
            blocks.append((ref_range, cur))
        cur, ref_range = [], None

    for line in open(path):
        if line.startswith("a"):
            flush()
            continue
        if not line.startswith("s "):
            continue
        pts = line.split()
        if len(pts) < 7:
            continue
        name = pts[1]
        sp, contig = name.split(".", 1) if "." in name else (name, "NA")
        start, size, strand, src = int(pts[2]), int(pts[3]), pts[4], int(pts[5])
        if sp == ref:
            if strand != "+":
                sys.exit(f"Error: {path}: reference block on the minus strand")
            ref_range = (start, start + size)
            continue
        if strand == "-":
            start = src - start - size
        cur.append((sp, contig, start, start + size))
    flush()
    return blocks


def main():
    ap = argparse.ArgumentParser(description="Step 17: per-species genomic span of the exons kept in each alignment")
    ap.add_argument("step9_trimmed_dir")
    ap.add_argument("step3_filtered_dir")
    ap.add_argument("bed_file")
    ap.add_argument("reference_species")
    ap.add_argument("out_tsv")
    ap.add_argument("--mask", help="step8c/step8c_missing_data.tsv")
    ap.add_argument("--tree", default=os.environ.get("TREE_FILE") or None, help="Newick tree (default $TREE_FILE)")
    a = ap.parse_args()
    for p in (a.step9_trimmed_dir, a.step3_filtered_dir):
        if not os.path.isdir(p):
            sys.exit(f"Error: missing folder {p}")
    for p in [a.bed_file] + [x for x in (a.mask, a.tree) if x]:
        if not os.path.isfile(p):
            sys.exit(f"Error: missing file {p}")
    ref = a.reference_species
    types = tree_node_types(a.tree)

    print("1. Parsing Reference BED file to map exons...")
    genes = {}
    with open(a.bed_file) as f:
        for line in f:
            cols = line.rstrip("\r\n").split("\t")
            # '#' comments and a header line without '#' (non-numeric start) are skipped
            if line.startswith("#") or len(cols) < 7 or not cols[1].strip().isdigit():
                continue
            tid = cols[6].strip().split("__")[0]
            genes.setdefault(tid, {"chrom": cols[0].strip(), "exons": []})["exons"].append(
                (int(cols[1]), int(cols[2])))
    for data in genes.values():
        data["min_s"] = min(s for s, _ in data["exons"])
        data["max_e"] = max(e for _, e in data["exons"])

    masked = collections.defaultdict(list)
    if a.mask:
        with open(a.mask) as f:
            next(f)
            for line in f:
                aln, sp, _exon, r0, r1, _reason = line.rstrip("\n").split("\t")
                masked[(aln, sp)].append((int(r0), int(r1)))
        masked = {k: merge(v) for k, v in masked.items()}
        print(f"   missing-data mask: {sum(len(v) for v in masked.values())} merged ranges for {len(masked)} "
              "alignment-species pairs")
    else:
        print("Warning: no --mask; exons step 8c gapped on the kept contig are not excluded", file=sys.stderr)

    rows = [["file_name", "species_name", "chromosome", "start", "end", "species_type"]]
    stats = collections.Counter()
    aln_files = sorted(glob.glob(os.path.join(a.step9_trimmed_dir, "*.aln")))
    print(f"2. Found {len(aln_files)} final alignment files. Scanning for multi-species coordinates...")
    for aln_path in aln_files:
        fname = os.path.basename(aln_path)
        stem = fname[:-4]
        clean = stem[len(ref) + 1:] if stem.startswith(ref + "_") else stem
        # longest BED transcript ID that prefixes the file name (IDs may share prefixes)
        tid = max((t for t in genes if clean.startswith(t + "_")), key=len, default=None)
        if tid is None:
            print(f"Warning: Could not link {fname} to BED file records. Skipping.")
            stats["alignments not linked to the BED"] += 1
            continue
        gene = genes[tid]

        kept_contigs = collections.defaultdict(set)
        for line in open(aln_path):
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                sp, contig = name.split(".", 1) if "." in name else (name, "NA")
                kept_contigs[sp].add(contig)
        coords = {sp: [] for sp in kept_contigs if sp != ref}

        for s0, e0 in gene["exons"]:
            maf_path = os.path.join(a.step3_filtered_dir, f"{ref}.{gene['chrom']}_{s0}_{e0}.txt")
            if not os.path.exists(maf_path):
                continue
            for ref_range, block in read_maf_blocks(maf_path, ref):
                for sp, contig, b0, b1 in block:
                    if sp not in coords:
                        continue
                    if ref_range is None:
                        stats["blocks skipped: no reference line"] += 1
                        continue
                    if ref_range[1] <= s0 or ref_range[0] >= e0:
                        stats["blocks skipped: only past the exon end"] += 1
                        continue
                    # clipped to the exon: a block running into the extra base past the exon end must still count
                    # as inside a whole-exon mask range
                    clipped = (max(ref_range[0], s0), min(ref_range[1], e0))
                    if a.mask and inside(masked.get((stem, sp), []), *clipped):
                        stats["blocks skipped: masked by step 8c"] += 1
                        continue
                    if contig not in kept_contigs[sp]:
                        if a.mask:
                            sys.exit(f"Error: {fname}: {sp} block on {contig} is outside the step 8c mask but the "
                                     f"species' row is on {','.join(sorted(kept_contigs[sp]))}")
                        stats["blocks skipped: not on the kept contig"] += 1
                        continue
                    coords[sp].append((contig, b0, b1))
                    stats["blocks used"] += 1

        for sp in sorted(kept_contigs):
            sp_type = "NA" if types is None else types.get(sp, "not_in_tree")
            if sp == ref:
                rows.append([fname, sp, gene["chrom"], str(gene["min_s"]), str(gene["max_e"]), sp_type])
                continue
            stats[f"species type: {sp_type}"] += 1
            if coords[sp]:
                chroms = sorted({c for c, _, _ in coords[sp]})
                if a.mask and len(chroms) > 1:
                    sys.exit(f"Error: {fname}: {sp} has kept blocks on {len(chroms)} contigs")
                rows.append([fname, sp, ",".join(chroms), str(min(b for _, b, _ in coords[sp])),
                             str(max(e for _, _, e in coords[sp])), sp_type])
            else:
                stats["species rows without kept blocks (NA)"] += 1
                rows.append([fname, sp, "NA", "NA", "NA", sp_type])

    print("3. Writing Master Coordinates Summary...")
    with open(a.out_tsv, "w") as f:
        for r in rows:
            f.write("\t".join(r) + "\n")
    for k, v in sorted(stats.items()):
        print(f"   {k}: {v:,}")
    print("-" * 50)
    print(f"Success! Processed {len(aln_files)} files and logged {len(rows) - 1} species coordinate entries.")
    print(f"Saved to: {a.out_tsv}")


if __name__ == "__main__":
    main()
