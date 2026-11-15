#!/usr/bin/env python3
"""
older_orf_pdng.py  (stage 4, step 7h; method option older_orf_makes_pdng)

A DNG whose focal lineage has an ancestor OLDER than its ORF node that also has the ORF (step 4b table) becomes a
PDNG (2026-11-15). The ORF node is the oldest ancestor of the unbroken run of ORF-bearing ancestors above the focal
species, so such an older ancestor is separated from the run by ancestors without the ORF: the ORF was present, lost
and present again, and the disablers placed below the ORF node then disable an ORF that already existed, which is not
the birth of a new one (Brassicaceae, default rules with N2 and N3 allowed as ORF nodes: 34 DNGs, 22 of them with the
older ORF at N2 or N3, e.g. VVA98552.1 with a complete ORF at N3, lost at N5 and N7, present in Arabis nemorensis).
With orf_clade_consistency, an ORF node moved younger has ORF-bearing ancestors above it as well, which this rule
treats the same way.

Rewrites Gene_type (DNG -> PDNG) in place in every <prefix>_nodes_summary_filtered.tsv and <prefix>_nodes_summary.tsv
of RESULTS/nodes_summary_by_species (the stage 4 final folder, after step 7g), and writes
<prefix>_older_orf_pdng.log (geneID, ORF_Node, the older ancestors with the ORF) and one line in RESULTS/PROVENANCE.txt.

Usage: python3 older_orf_pdng.py RESULTS STAGE4_INPUTS_TSV
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lingua_config as cfg  # noqa: E402

csv.field_size_limit(10 ** 9)


def die(m):
    sys.exit(f"ERROR: {m}")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    res, inputs = sys.argv[1:3]
    nd = os.path.join(res, "nodes_summary_by_species")
    for p in (nd, inputs):
        if not os.path.exists(p):
            die(f"not found: {p}")
    tot_in = tot_re = 0
    for line in open(inputs):
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[7] == "NA":
            continue
        F, npref, anc = f[0], f[3], f[7]
        lineage = sorted([n for n, c in cfg.NODE_CLADES.items() if F in c], key=lambda n: -len(cfg.NODE_CLADES[n]))
        cells = {}
        for r in csv.DictReader(open(anc), delimiter="\t"):
            g = r["geneID"].strip()
            on = {n for n in lineage if (r.get(n) or "").strip() in ("1", "1.0")}
            cells[g] = on
            cells.setdefault(g[:-2] if g.endswith(".p") else g, on)
        log = []
        for suffix in ("_nodes_summary_filtered.tsv", "_nodes_summary.tsv"):
            path = os.path.join(nd, npref + suffix)
            if not os.path.exists(path):
                continue
            rd = csv.DictReader(open(path), delimiter="\t")
            fields, rows = rd.fieldnames, list(rd)
            n_in = n_re = 0
            for r in rows:
                if r["Gene_type"] != "DNG":
                    continue
                n_in += 1
                g = r["geneID"].strip()
                if g not in cells:
                    die(f"{g} of {path} is not in {anc}")
                orf = r["ORF_Node"]
                k = lineage.index(orf) if orf in lineage else len(lineage)
                if orf not in lineage and orf != f"{F}-specific":
                    die(f"{g}: ORF node {orf} is not on the lineage of {F}")
                older = [n for n in lineage[:k] if n in cells[g]]
                if older:
                    r["Gene_type"] = "PDNG"
                    n_re += 1
                    if suffix == "_nodes_summary_filtered.tsv":
                        log.append((g, orf, ",".join(older)))
            tmp = path + ".tmp"
            with open(tmp, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
            os.replace(tmp, path)
            if suffix == "_nodes_summary_filtered.tsv":
                tot_in += n_in
                tot_re += n_re
                print(f"  {npref + suffix:<44}{n_in:>8,} DNG in, {n_re:>5,} -> PDNG (older ancestor with the ORF)", file=sys.stderr)
        with open(os.path.join(nd, f"{npref}_older_orf_pdng.log"), "w") as fh:
            fh.write("geneID\tORF_Node\tolder_ancestors_with_ORF\n")
            fh.writelines("\t".join(x) + "\n" for x in log)
    with open(os.path.join(res, "PROVENANCE.txt"), "a") as fh:
        fh.write(f"older_orf_pdng.py (step 7h): {tot_in} DNG in, {tot_re} -> PDNG because an ancestor older than the ORF "
                 f"node has the ORF (filtered tables)\n")
    print(f"  TOTAL: {tot_in} DNG in, {tot_re} -> PDNG", file=sys.stderr)


if __name__ == "__main__":
    main()
