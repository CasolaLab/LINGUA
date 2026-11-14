#!/usr/bin/env python3
"""
filter_enabler_events.py  (stage 4, step 6b; method options drop_terminal_indels, homopolymer_needs_extant)

Remove indel events that are not evidence of a mutation from the step 6 event list before step 7 chooses enablers
(T2T mammals, his decisions 2026-11-12):
  drop_terminal_indels      an indel whose gap run, in any carrier ancestor, reaches the first or the last base of the
                            focal CDS: the ancestor has no aligned sequence there, so a deletion cannot be told from
                            unaligned flanking sequence (hand check: 4 of 4 such enabler events were unaligned starts;
                            enabler indels 1-8 nt from an end showed no excess over the background)
  homopolymer_needs_extant  an indel in a homopolymer (>= 4 identical bases within 4 nt of the event in the focal CDS, the
                            definition of annotate_evidence.py) counts only when an extant species on the side branch at one
                            of its carrier nodes (descendants of that node not descending from the next younger lineage
                            node) carries the same indel (same start, same length)
  homopolymer_support_any_outgroup  (with homopolymer_needs_extant) the supporting extant species may be any species outside
                            the clade of the lineage node just younger than the youngest carrier (i.e. outside the clade
                            that lacks the event), not only the carriers' side branches (his decision 2026-11-12; e.g. a gap
                            shared by goat and sheep supports a primate N2/N4 event when the root is not a carrier)
Stop events are not touched. The step 6 list is kept as all_indels_nodes_<F>.unfiltered.tsv (read from there on reruns) and
every removed event is written to indels_dir/removed_enabler_events_<F>.tsv with its reason.

Usage: python3 filter_enabler_events.py INDELS_DIR ALN_DIR --focal-species Homo_sapiens
"""
import argparse
import collections
import csv
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lingua_config as cfg  # noqa: E402


def die(m):
    sys.exit(f"ERROR: {m}")


def read_aln(path):
    rows, k = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split(".")[0]
            rows.setdefault(k, "")
        elif k:
            rows[k] += line.upper()
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("indels_dir")
    ap.add_argument("aln_dir")
    ap.add_argument("--focal-species", required=True)
    ap.add_argument("--detail", default=None, help="indel detail table (default <indels_dir>/<F>_DNG_indels_summary.tsv)")
    a = ap.parse_args()
    F = a.focal_species
    terminal_on, homo_on = cfg.opt("drop_terminal_indels"), cfg.opt("homopolymer_needs_extant")
    any_outgroup = cfg.opt("homopolymer_support_any_outgroup")
    events_path = os.path.join(a.indels_dir, f"all_indels_nodes_{F}.tsv")
    backup = events_path[:-4] + ".unfiltered.tsv"
    detail_path = a.detail or os.path.join(a.indels_dir, f"{F}_DNG_indels_summary.tsv")
    for p in (detail_path, a.aln_dir):
        if not os.path.exists(p):
            die(f"not found: {p}")
    if not os.path.exists(backup):
        if not os.path.exists(events_path):
            die(f"not found: {events_path}")
        shutil.copy2(events_path, backup)
    events = list(csv.DictReader(open(backup), delimiter="\t"))
    fields = list(events[0].keys()) if events else ["geneID", "IndelStart_RefPos", "Nodes_indels"]

    detail = collections.defaultdict(list)
    for r in csv.DictReader(open(detail_path), delimiter="\t"):
        detail[(r["gff_transcriptID"], r["IndelStart_RefPos"])].append(r)
    lineage = sorted([n for n, c in cfg.NODE_CLADES.items() if F in c], key=lambda n: -len(cfg.NODE_CLADES[n]))
    species = set(cfg.SPECIES)

    def side_branch(node):
        i = lineage.index(node)
        younger = set(cfg.NODE_CLADES[lineage[i + 1]]) if i + 1 < len(lineage) else {F}
        return (set(cfg.NODE_CLADES[node]) - younger) & species

    def outside_non_carrier_clade(nodes):
        youngest = max((n for n in nodes if n in lineage), key=lineage.index)
        i = lineage.index(youngest)
        younger = set(cfg.NODE_CLADES[lineage[i + 1]]) if i + 1 < len(lineage) else {F}
        return species - younger

    cds_cache = {}

    def focal_cds(g):
        if g not in cds_cache:
            gk = g[:-2] if g.endswith(".p") else g  # Phytozome IDs carry ".p" in the protein tables, not in alignments
            p = glob.glob(os.path.join(a.aln_dir, f"{F}_{gk}_*.aln"))
            if len(p) != 1:
                die(f"{len(p)} alignments for {g} in {a.aln_dir}")
            cds_cache[g] = read_aln(p[0])[F].replace("-", "")
        return cds_cache[g]

    kept, removed = [], []
    reasons = collections.Counter()
    for r in events:
        g, pos = r["geneID"], r["IndelStart_RefPos"]
        nodes = r["Nodes_indels"].split("|")
        rows = [x for x in detail.get((g, pos), []) if x["Species"] in nodes]
        if not rows:
            die(f"no step 14 detail rows for {g} at {pos} in {nodes}")
        why = []
        if terminal_on and any(int(x["IndelStart_RefPos"]) == 1 or int(x["IndelEnd_RefPos"]) == int(x["reference_length"])
                               for x in rows):
            why.append("terminal_gap")
        if homo_on:
            p = int(pos)
            w = focal_cds(g)[max(0, p - 5):p + 4]
            if any(w[j:j + 4] == b * 4 for b in "ACGT" for j in range(len(w) - 3)):
                lens = {x["IndelLength"] for x in rows}
                carriers = {x["Species"] for x in detail[(g, pos)] if x["IndelLength"] in lens}
                if any_outgroup:
                    side = outside_non_carrier_clade(nodes)
                else:
                    side = set().union(*[side_branch(n) for n in nodes if n in lineage]) if nodes else set()
                if not carriers & side:
                    why.append("homopolymer_no_extant_support")
        if why:
            removed.append({**r, "reason": ",".join(why)})
            reasons[",".join(why)] += 1
        else:
            kept.append(r)
    with open(events_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(kept)
    with open(os.path.join(a.indels_dir, f"removed_enabler_events_{F}.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields + ["reason"], delimiter="\t")
        w.writeheader()
        w.writerows(removed)
    assert len(kept) + len(removed) == len(events)
    print(f"{F}: {len(events)} indel events, {len(removed)} removed {dict(reasons)}, {len(kept)} kept", file=sys.stderr)


if __name__ == "__main__":
    main()
