#!/usr/bin/env python3
"""
parsimony_pdng.py  (stage 4, step 7i; method option orf_clade_parsimony_pdng; T2T 2026-11-15, his rule)

A DNG becomes a PDNG when parsimony on the living species does not place the ORF-compatible state at its ORF node for every chosen enabler
event. Meant to run with orf_clade_living_parent_orf "living_never" (the ORF node is never moved by the ORF-clade rule): a living species
inside the ORF clade that carries the disabler is then judged by parsimony instead of moving the node or being ignored.

Per DNG of RESULTS/nodes_summary_by_species/<npref>_nodes_summary_filtered.tsv and per chosen event (Indels_event_RefPos, Stops_event_RefPos):
  living states  1  carries the event: a row of the step 13/14 detail table at the same position (indels: the length carried by the event's
                    nearest lineage node older than the ORF node among Indels_Node / Stops_Node);
                 0  not a carrier and the site is aligned: stops and start codons, the three bases at the focal codon of the event are all
                    present (not '-' or N); indels, the row has bases at every focal position of the event interval (pos .. pos + length - 1);
                 ?  otherwise (missing data). The focal species is 0.
  parsimony      Sankoff, unit cost, two states, on the tree of the tree config (living species as tips, '?' = either state); the set of
                 states at the ORF node over all most-parsimonious assignments. The gene stays DNG only if that set is {0} for every event.
Living species only: the Cactus ancestral sequences are not used. Rewrites Gene_type DNG -> PDNG in place in <npref>_nodes_summary_filtered.tsv
and <npref>_nodes_summary.tsv, writes <npref>_parsimony_pdng.log (one row per DNG and event: tip states, ORF-node set, verdict) and a line in
RESULTS/PROVENANCE.txt.

Usage: python3 parsimony_pdng.py RESULTS STAGE4_INPUTS_TSV
  STAGE4_INPUTS_TSV columns: focal, qpref, opref, npref, aln_dir, indel_table, stop_table, ancseq (as run_pipeline.sh writes it)
"""
import collections
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lingua_config as cfg  # noqa: E402

csv.field_size_limit(10 ** 9)
INF = float("inf")


def die(m):
    sys.exit(f"ERROR: {m}")


def read_aln(p):
    rows, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split(".")[0]
            rows[k] = []
        elif k:
            rows[k].append(line.upper())
    return {a: "".join(b) for a, b in rows.items()}


def build_tree():
    clades = {n: frozenset(c) for n, c in cfg.NODE_CLADES.items()}
    children = {}
    for n, c in clades.items():
        inner = [m for m, d in clades.items() if d < c]
        direct = [m for m in inner if not any(clades[m] < clades[o] for o in inner)]
        covered = set().union(*(clades[m] for m in direct)) if direct else set()
        children[n] = direct + sorted(c - covered)
    root = max(clades, key=lambda n: len(clades[n]))
    return root, children


def mp_sets(root, children, tips):
    """Sankoff: node -> frozenset of states present in some most-parsimonious assignment. tips: species -> allowed states."""
    order, stack = [], [root]
    while stack:
        n = stack.pop()
        order.append(n)
        stack += children.get(n, [])
    down = {}
    for n in reversed(order):
        if n not in children:
            ok = tips.get(n, (0, 1))
            down[n] = [0 if s in ok else INF for s in (0, 1)]
        else:
            down[n] = [sum(min(down[c][t] + (s != t) for t in (0, 1)) for c in children[n]) for s in (0, 1)]
    up = {root: [0, 0]}
    for n in order:
        for c in children.get(n, []):
            rest = [[min(down[d][u] + (s != u) for u in (0, 1)) for d in children[n] if d != c] for s in (0, 1)]
            up[c] = [min(up[n][s] + sum(rest[s]) + (s != t) for s in (0, 1)) for t in (0, 1)]
    best = min(down[root])
    return {n: frozenset(s for s in (0, 1) if down[n][s] + up[n][s] == best) for n in down}


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    res, inputs = sys.argv[1:3]
    nd = os.path.join(res, "nodes_summary_by_species")
    for p in (nd, inputs):
        if not os.path.exists(p):
            die(f"not found: {p}")
    root, children = build_tree()
    living = set(cfg.SPECIES)
    tot_in = tot_re = 0
    for line in open(inputs):
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[7] == "NA":
            continue
        F, npref, alnd, indt, stpt = f[0], f[3], f[4], f[5], f[6]
        for p in (alnd, indt, stpt):
            if not os.path.exists(p):
                die(f"not found: {p}")
        lineage_old = sorted([n for n, c in cfg.NODE_CLADES.items() if F in c], key=lambda n: -len(cfg.NODE_CLADES[n]))
        young = [f"{F}-specific"] + lineage_old[::-1]
        age = {n: i for i, n in enumerate(young)}
        det = collections.defaultdict(lambda: collections.defaultdict(set))
        for r in csv.DictReader(open(indt), delimiter="\t"):
            det[(r["gff_transcriptID"], "indel", r["IndelStart_RefPos"])][r["Species"].split(".")[0]].add(r["IndelLength"])
        for r in csv.DictReader(open(stpt), delimiter="\t"):
            det[(r["gff_transcriptID"], "stop", r["StopCodonStart_RefPos"])][r["Species"].split(".")[0]].add("stop")
        verdict, log = {}, []
        path = os.path.join(nd, npref + "_nodes_summary_filtered.tsv")
        if not os.path.exists(path):
            continue
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["Gene_type"] != "DNG":
                continue
            g, orf = r["geneID"].strip(), r["ORF_Node"]
            gk = g[:-2] if g.endswith(".p") else g
            ap = glob.glob(os.path.join(alnd, f"{F}_{gk}_*.aln"))
            if len(ap) != 1:
                die(f"{len(ap)} alignments for {g} in {alnd}")
            aln = read_aln(ap[0])
            fc = [i for i, ch in enumerate(aln[F]) if ch != "-"]
            keep = True
            for kind, pc, nc in (("indel", "Indels_event_RefPos", "Indels_Node"), ("stop", "Stops_event_RefPos", "Stops_Node")):
                pos = r.get(pc, "NA")
                if pos in ("NA", ""):
                    continue
                car = det.get((gk, kind, pos), {})  # Phytozome IDs carry ".p" in the summaries, not in the detail tables
                if not car:
                    die(f"{g}: no detail rows for the chosen {kind} at {pos}")
                older = sorted((age[n], n) for n in (r.get(nc) or "").split("|") if n in age and age[n] > age.get(orf, -1) and n in car)
                lengths = car[older[0][1]] if older else set()
                if not lengths:
                    die(f"{g}: the chosen {kind} at {pos} has no lineage carrier older than {orf}")
                p0 = int(pos) - 1
                span = 3 if kind == "stop" else max(int(x) for x in lengths)
                site = [fc[j] for j in range(p0, min(p0 + span, len(fc)))]
                tips = {}
                for sp in living:
                    if sp == F:
                        tips[sp] = (0,)
                    elif car.get(sp, set()) & lengths:
                        tips[sp] = (1,)
                    elif sp in aln and site and all(aln[sp][j] not in "-N" for j in site):
                        tips[sp] = (0,)
                    else:
                        tips[sp] = (0, 1)
                sets = mp_sets(root, children, tips)
                at = sets[orf] if orf in sets else frozenset({0})
                ok = at == frozenset({0})
                keep = keep and ok
                log.append((g, orf, kind, pos, ",".join(sorted(sp for sp, v in tips.items() if v == (1,))),
                            ",".join(sorted(sp for sp, v in tips.items() if v == (0, 1))), "".join(map(str, sorted(at))),
                            "intact" if ok else "not_intact"))
            verdict[g] = keep
        for suffix in ("_nodes_summary_filtered.tsv", "_nodes_summary.tsv"):
            p = os.path.join(nd, npref + suffix)
            if not os.path.exists(p):
                continue
            rd = csv.DictReader(open(p), delimiter="\t")
            fields, rows = rd.fieldnames, list(rd)
            n_in = n_re = 0
            for x in rows:
                g = x["geneID"].strip()
                if x["Gene_type"] == "DNG" and g in verdict:
                    n_in += 1
                    if not verdict[g]:
                        x["Gene_type"] = "PDNG"
                        n_re += 1
            tmp = p + ".tmp"
            with open(tmp, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
            os.replace(tmp, p)
            if suffix == "_nodes_summary_filtered.tsv":
                tot_in += n_in
                tot_re += n_re
                print(f"  {npref + suffix:<44}{n_in:>8,} DNG in, {n_re:>5,} -> PDNG (ORF node not intact by parsimony)", file=sys.stderr)
        with open(os.path.join(nd, f"{npref}_parsimony_pdng.log"), "w") as fh:
            fh.write("geneID\tORF_Node\tkind\tposition\tliving_carriers\tliving_missing\torf_node_states\tverdict\n")
            fh.writelines("\t".join(x) + "\n" for x in log)
    with open(os.path.join(res, "PROVENANCE.txt"), "a") as fh:
        fh.write(f"parsimony_pdng.py (step 7i): {tot_in} DNG in, {tot_re} -> PDNG because parsimony on the living species does not make the "
                 f"ORF node intact for every chosen enabler event (filtered tables)\n")
    print(f"  TOTAL: {tot_in} DNG in, {tot_re} -> PDNG", file=sys.stderr)


if __name__ == "__main__":
    main()
