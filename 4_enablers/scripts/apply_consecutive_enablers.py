"""
apply_consecutive_enablers.py  (LINGUA stage 4, final filter; his rule of 2026-10-29)

A DNG must have its two enabler nodes equal to the two lineage nodes immediately older than its ORF node (the parent
and the grandparent of ORF_Node on the focal lineage). Step 7 picks the two nearest older nodes carrying the chosen
disablers, which can leave a gap (e.g. N2|N14 with N3, N5, N8, N11 lacking the disabler); ORF-Enablers_Node_dist <= 1
already forces the younger enabler to be the parent, this rule adds that the older one is the grandparent.
Writes a parallel results directory with only Gene_type rewritten (DNG -> PDNG when the pair is not consecutive), as
apply_option_b_reclassification.py does, so the set is strictly nested in the option B set.

Usage:
  python3 apply_consecutive_enablers.py <src_results> <dst_results>
  (tree configuration from ../data/tree_config.json, or the file named by LINGUA_CONFIG)
  (src_results holds nodes_summary_by_species/ with *_nodes_summary_filtered.tsv, *_nodes_summary.tsv, *_evidence.tsv)
"""
import csv
import glob
import json
import os
import shutil
import sys
from collections import Counter

csv.field_size_limit(10 ** 9)
CFG = os.environ.get("LINGUA_CONFIG", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                                  "data", "tree_config.json"))


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1:3]
    sd, dd = os.path.join(src, "nodes_summary_by_species"), os.path.join(dst, "nodes_summary_by_species")
    if not os.path.isdir(sd):
        sys.exit(f"ERROR: missing {sd}")
    if os.path.exists(dd):
        sys.exit(f"ERROR: {dd} exists; refusing to overwrite")
    cfg = json.load(open(CFG))
    clades = {n: set(v) for n, v in cfg["node_clades"].items()}
    aliases = cfg["aliases"]
    os.makedirs(dd)
    tot = Counter()
    print(f"  {'file':<44}{'DNG in':>8}{'DNG out':>9}{'->PDNG':>8}", file=sys.stderr)
    for f in sorted(glob.glob(os.path.join(sd, "*"))):
        base = os.path.basename(f)
        if not base.endswith("_nodes_summary_filtered.tsv") and not base.endswith("_nodes_summary.tsv"):
            shutil.copy2(f, os.path.join(dd, base))
            continue
        pre = base.split("_nodes_summary")[0]
        F = aliases.get(pre, aliases.get(pre.lower(), pre))
        lineage = sorted([n for n in clades if F in clades[n]], key=lambda n: -len(clades[n])) + [f"{F}-specific"]
        rows = list(csv.DictReader(open(f), delimiter="\t"))
        n_in = n_out = n_re = 0
        with open(os.path.join(dd, base), "w", newline="") as oh:
            wr = csv.DictWriter(oh, delimiter="\t", fieldnames=list(rows[0].keys()) if rows else [], lineterminator="\n")
            wr.writeheader()
            for r in rows:
                if r["Gene_type"] == "DNG":
                    n_in += 1
                    if r["ORF_Node"] not in lineage:
                        sys.exit(f"ERROR: {base} {r['geneID']}: ORF node {r['ORF_Node']} not on the lineage of {F}")
                    o = lineage.index(r["ORF_Node"])
                    required = {lineage[o - 1], lineage[o - 2]} if o >= 2 else set()
                    if set(r["Enablers_Nodes"].split("|")) == required:
                        n_out += 1
                    else:
                        r["Gene_type"] = "PDNG"
                        n_re += 1
                wr.writerow(r)
        assert n_in == n_out + n_re
        if base.endswith("_filtered.tsv"):
            print(f"  {base:<44}{n_in:>8,}{n_out:>9,}{n_re:>8,}", file=sys.stderr)
            tot["in"] += n_in
            tot["out"] += n_out
            tot["re"] += n_re
    print(f"  {'TOTAL (filtered tables, all species)':<44}{tot['in']:>8,}{tot['out']:>9,}{tot['re']:>8,}", file=sys.stderr)
    with open(os.path.join(dst, "PROVENANCE.txt"), "w") as fh:
        fh.write(f"apply_consecutive_enablers.py on {src}: DNG kept only when Enablers_Nodes == parent + grandparent of ORF_Node; "
                 f"{tot['in']} DNG in, {tot['out']} out, {tot['re']} -> PDNG (filtered tables)\n")


if __name__ == "__main__":
    main()
