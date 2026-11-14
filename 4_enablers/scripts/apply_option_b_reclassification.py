"""
Option B. Keep a locus as a DNG only when BOTH of its enabler nodes are corroborated by a
LIVING genome, and reclassify the rest as PDNGs.

THE RULE, which is a conjunction and strictly narrows the existing one:
  condition 1, unchanged, from the ancestral sequences -- the same disabler is reconstructed
    at the two ancestral nodes immediately older than the ORF node, with
    ORF-Enablers_Node_dist <= 1. All current DNGs already satisfy this.
  condition 2, what Option B adds, from living genomes -- for EACH of those two nodes, at
    least one extant species on that node's side branch (descendants of the node that do not
    descend from the next younger lineage node) carries the SAME event, same position and for
    indels the same length.
Condition 2 is read from `enabler_nodes_without_extant_support` in the step-7d evidence table,
which annotate_evidence.py builds with exactly that geometry and that position/length match.
Every Option B DNG is therefore also a DNG under the old rule; the set is strictly nested.

WHY. Under the competing history -- the ORF is ancestral and was lost in the other lineages --
the same mutation at the same position would have to arise independently in two lineages that
diverge at DIFFERENT depths on the focal path. Observing it in both makes that implausible, so
parsimony fixes the disabler as ancestral and the enabler as derived on the focal branch.
Condition 1 alone settles polarity by reconstruction; condition 2 settles it by observation.
annotate_evidence.py puts it in six words -- "Extant genomes are observations, ancestors
reconstructions."

WHAT THIS SCRIPT DOES, and deliberately does not do. It writes a PARALLEL results directory
with only the Gene_type column rewritten, so every downstream script runs unchanged against
it and stage 4 itself stays untouched and reproducible. It does not re-run any alignment, any
node call, or any search. TDNG per locus is preserved, so the classified loci are only
reapportioned between the two classes.

Usage:
    python3 apply_option_b_reclassification.py <src_results> <dst_results>
(src_results must hold nodes_summary_by_species/ with the step-7d <prefix>_evidence.tsv files
next to the <prefix>_nodes_summary_filtered.tsv files)
"""

import csv
import glob
import os
import shutil
import sys
from collections import Counter

csv.field_size_limit(10 ** 9)

def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1:3]
    sd = os.path.join(src, "nodes_summary_by_species")
    dd = os.path.join(dst, "nodes_summary_by_species")
    if not os.path.isdir(sd):
        sys.exit(f"ERROR: missing {sd}")
    os.makedirs(dd, exist_ok=True)
    for extra in ("indels_dir", "stops_dir", "tree_config.json"):
        s = os.path.join(src, extra)
        d = os.path.join(dst, extra)
        if os.path.exists(s) and not os.path.exists(d):
            (shutil.copytree if os.path.isdir(s) else shutil.copy2)(s, d)

    tot = Counter()
    print(f"  {'species':<30}{'DNG in':>8}{'DNG out':>9}{'->PDNG':>8}", file=sys.stderr)
    # T2T: species are the *_nodes_summary_filtered.tsv files present (the Brassicaceae list was hardcoded);
    # every other file of the folder (ETE, conserved-ORF, BLASTp tables...) is copied unchanged
    species = sorted(os.path.basename(f)[: -len("_nodes_summary_filtered.tsv")]
                     for f in glob.glob(os.path.join(sd, "*_nodes_summary_filtered.tsv")))
    if not species:
        sys.exit(f"ERROR: no *_nodes_summary_filtered.tsv in {sd}")
    for f in glob.glob(os.path.join(sd, "*")):
        b = os.path.basename(f)
        if not (b.endswith("_nodes_summary_filtered.tsv") or b.endswith("_evidence.tsv")):
            shutil.copy2(f, os.path.join(dd, b))
    for sp in species:
        f_sum = os.path.join(sd, f"{sp}_nodes_summary_filtered.tsv")
        f_ev = os.path.join(sd, f"{sp}_evidence.tsv")
        for p in (f_sum, f_ev):
            if not os.path.isfile(p):
                sys.exit(f"ERROR: missing {p}")
        ev = {r["geneID"]: r for r in csv.DictReader(open(f_ev), delimiter="\t")}
        rows = list(csv.DictReader(open(f_sum), delimiter="\t"))
        miss = [r["geneID"] for r in rows if r["geneID"] not in ev]
        if miss:
            sys.exit(f"ABORT: {sp}: {len(miss)} loci have no evidence row, e.g. {miss[:3]}; "
                     f"the reclassification would silently keep them")
        n_in = n_out = 0
        for r in rows:
            if r["Gene_type"] != "DNG":
                continue
            n_in += 1
            without = str(ev[r["geneID"]]["enabler_nodes_without_extant_support"]).strip()
            if without.lower() == "none":
                n_out += 1
            else:
                r["Gene_type"] = "PDNG"
        with open(os.path.join(dd, f"{sp}_nodes_summary_filtered.tsv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        shutil.copy2(f_ev, os.path.join(dd, f"{sp}_evidence.tsv"))
        c = Counter(r["Gene_type"] for r in rows)
        tot.update(c)
        tot["dng_in"] += n_in; tot["dng_out"] += n_out
        tot["rows_in"] += len(rows)
        print(f"  {sp:<30}{n_in:>8,}{n_out:>9,}{n_in-n_out:>8,}", file=sys.stderr)

    print(f"  {'TOTAL':<30}{tot['dng_in']:>8,}{tot['dng_out']:>9,}"
          f"{tot['dng_in']-tot['dng_out']:>8,}", file=sys.stderr)
    print(f"\n  reclassified table: DNG {tot['DNG']:,}  PDNG {tot['PDNG']:,}  "
          f"TDNG {tot['DNG']+tot['PDNG']:,}", file=sys.stderr)
    if tot["DNG"] + tot["PDNG"] != tot["rows_in"]:
        sys.exit(f"ABORT: TDNG changed, {tot['DNG'] + tot['PDNG']:,} classified vs {tot['rows_in']:,} rows in")
    if tot["DNG"] != tot["dng_out"]:
        sys.exit("ABORT: DNG count does not match the retained count")
    print(f"  wrote {dst}", file=sys.stderr)


if __name__ == "__main__":
    main()
