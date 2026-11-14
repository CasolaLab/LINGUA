#!/usr/bin/env python3
"""
remove_conserved_orfs.py  (stage 4, step 7c)

Remove candidate de novo genes, DNG and PDNG alike, whose ORF is present in any of the listed
ancestral genomes: a candidate with an ORF in an ancestor outside (or at the base of) the ingroup
cannot be shown to be de novo, so it is dropped regardless of its Gene_type. ORF presence is the
per-node column of the step 4 *_AncSeq.tsv table (1 = a tBLASTn HSP in that ancestral sequence
passing the step 4 thresholds, frac_ORFss and pident). Every node is checked, on or off the focal
lineage.

Reads and rewrites <out_prefix>_nodes_summary_filtered.tsv (the step 7b output that step 8 reads)
and writes <out_prefix>_conserved_orf_removed.tsv with one row per removed gene and the nodes that
caused it. run_pipeline.sh runs it right after step 7b, which regenerates the filtered file, so the
log always describes one pass.

Nodes are given by their original (Cactus) names or N labels and are resolved through
data/tree_config.json.

Usage:
    python3 remove_conserved_orfs.py NODES_DIR ANCSEQ_TSV OUT_PREFIX NODE [NODE ...]
    e.g. ... nodes_summary_by_species cdngs_tblastn_results/hsapiens_cdngs_tblastn_150_AncSeq.tsv \
             hsapiens fullTreeAnc238 fullTreeAnc115 fullTreeAnc235 fullTreeAnc178 fullTreeAnc38
"""
import csv
import os
import sys

import lingua_config as cfg


def die(msg):
    sys.exit(f"ERROR: {msg}")


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    nodes_dir, ancseq, prefix = sys.argv[1:4]
    requested = sys.argv[4:]

    orig = cfg.CFG["node_original_names"]
    back = {v: k for k, v in orig.items()}
    labels = []
    for n in requested:
        lab = n if n in orig else back.get(n)
        if lab is None:
            die(f"node {n!r} is neither an N label nor an original node name in the tree config")
        labels.append(lab)
    name = {lab: f"{orig[lab]}({lab})" for lab in labels}

    summ = os.path.join(nodes_dir, f"{prefix}_nodes_summary_filtered.tsv")
    for p in (summ, ancseq):
        if not os.path.isfile(p):
            die(f"not found: {p}")

    with open(ancseq, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        missing_cols = [lab for lab in labels if lab not in rd.fieldnames]
        if missing_cols:
            die(f"{ancseq} has no column for {missing_cols}")
        orf = {r["geneID"].strip(): r for r in rd}

    with open(summ, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        fields = rd.fieldnames
        rows = list(rd)

    kept, removed = [], []
    for r in rows:
        g = r["geneID"].strip()
        a = orf.get(g)
        if a is None:
            die(f"{g} is in {summ} but not in {ancseq}")
        hits = [lab for lab in labels if a[lab].strip() == "1"]
        if hits:
            removed.append((g, r.get("Gene_type", ""), r.get("ORF_Node", ""),
                            ",".join(name[h] for h in hits)))
        else:
            kept.append(r)
    if len(kept) + len(removed) != len(rows):
        die("kept + removed does not equal the input rows")

    tmp = summ + ".tmp"
    with open(tmp, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        wr.writeheader()
        wr.writerows(kept)
    os.replace(tmp, summ)

    log = os.path.join(nodes_dir, f"{prefix}_conserved_orf_removed.tsv")
    with open(log, "w") as fh:
        fh.write("geneID\tGene_type\tORF_Node\tORF_present_in\n")
        for r in removed:
            fh.write("\t".join(r) + "\n")

    by_type = {}
    for _, t, _, _ in removed:
        by_type[t] = by_type.get(t, 0) + 1
    print(f"{os.path.basename(summ)}: in={len(rows)}, kept={len(kept)}, removed={len(removed)} "
          f"{by_type} (ORF in any of {', '.join(name[l] for l in labels)})", file=sys.stderr)


if __name__ == "__main__":
    main()
