#!/usr/bin/env python3
"""
annotate_evidence.py  (stage 4, step 7d)

Annotate, without filtering, the quality of the enabler evidence behind every classified locus.
Needs the event positions written by step 7 (method option record_enabler_events). For each event
used (the indel and the stop rows whose nodes entered Enablers_Nodes), it records:
    homopolymer       (indels) >= 4 identical bases within 4 nt of the event in the focal CDS
    cds_edge          event within 9 nt of either end of the focal CDS
    N_nearby          an N within 5 alignment columns in the row of any carrier ancestor
    multi_length      (indels) the carrier ancestors do not all share one indel length
and per locus whether its enabler nodes rest only on flagged events (evidence_flagged).
NEE revision 2026-10-28: an indel carrier must match the length of the enabler node's own copy; with
--all-indel-events/--all-stop-events an enabler node is supported by ANY valid disabler at that node (absent from
every node and living species of the ORF clade, carried at the same position by a living side-branch species),
including nodes step 7 added from living evidence (column Enablers_living_inferred).
With --stops-detail it also records extant support (item 4): for each enabler node X, the extant
species on the side branch at X (descendants of X that do not descend from the next younger
lineage node) that carry the same event (same position; for indels also the same length), and
those that contradict it (aligned bases within 10 columns of the event, but not the event); species
with no aligned base there are no data. Extant genomes are observations, ancestors reconstructions.

Output: <out_prefix>_evidence.tsv in NODES_DIR, one row per locus of the filtered summary.

Usage: python3 annotate_evidence.py NODES_DIR OUT_PREFIX RESULTS_ALN_DIR INDELS_DETAIL_TSV
                                    --focal-species Homo_sapiens [--stops-detail STOPS_DETAIL_TSV]
"""
import argparse
import collections
import csv
import glob
import os
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
    ap.add_argument("nodes_dir")
    ap.add_argument("out_prefix")
    ap.add_argument("aln_dir")
    ap.add_argument("indels_detail")
    ap.add_argument("--focal-species", required=True)
    ap.add_argument("--stops-detail", default=None)
    ap.add_argument("--all-indel-events", default=None,
                    help="step-6 all_indels_nodes table; with --all-stop-events, an enabler node is supported by ANY "
                         "valid disabler, not only the event step 7 chose (NEE rule 5)")
    ap.add_argument("--all-stop-events", default=None)
    a = ap.parse_args()
    if bool(a.all_indel_events) != bool(a.all_stop_events) or (a.all_indel_events and not a.stops_detail):
        die("--all-indel-events and --all-stop-events go together, and need --stops-detail")
    F = a.focal_species
    summ_path = os.path.join(a.nodes_dir, f"{a.out_prefix}_nodes_summary_filtered.tsv")
    for p in (summ_path, a.aln_dir, a.indels_detail):
        if not os.path.exists(p):
            die(f"not found: {p}")
    rows = list(csv.DictReader(open(summ_path), delimiter="\t"))
    if rows and "Indels_event_RefPos" not in rows[0]:
        die("event positions missing; switch on record_enabler_events (step 7)")

    lengths = collections.defaultdict(lambda: collections.defaultdict(set))
    for r in csv.DictReader(open(a.indels_detail), delimiter="\t"):
        lengths[(r["gff_transcriptID"], r["IndelStart_RefPos"])][r["IndelLength"]].add(r["Species"])

    stop_carriers = collections.defaultdict(set)
    if a.stops_detail:
        for r in csv.DictReader(open(a.stops_detail), delimiter="\t"):
            stop_carriers[(r["gff_transcriptID"], r["StopCodonStart_RefPos"])].add(r["Species"])
    all_events = collections.defaultdict(list)
    if a.all_indel_events:
        for kind, path, pcol, ncol in (("indel", a.all_indel_events, "IndelStart_RefPos", "Nodes_indels"),
                                       ("stop", a.all_stop_events, "StopCodonStart_RefPos", "Nodes_stops")):
            for r in csv.DictReader(open(path), delimiter="\t"):
                g0 = r["geneID"][:-2] if r["geneID"].endswith(".p") else r["geneID"]
                all_events[g0].append((kind, r[pcol], set(r[ncol].split("|")) - {"NA", ""}))
    lineage = sorted([n for n, c in cfg.NODE_CLADES.items() if F in c], key=lambda n: -len(cfg.NODE_CLADES[n]))
    species = set(cfg.SPECIES)

    def side_branch(node):
        i = lineage.index(node)
        younger = set(cfg.NODE_CLADES[lineage[i + 1]]) if i + 1 < len(lineage) else {F}
        return set(cfg.NODE_CLADES[node]) - younger

    def node_lengths(kind, gk, pos, node, nodes, inferred_nodes=()):
        """Lengths of the event at `node` (indels), else at the nearest older lineage node of the event that
        carries it (living-inferred nodes have no reconstructed copy of the event: their own reconstructed indel at
        the same start, when it has another length, is a different event and is not used; fix 2026-11-13); 'stop'
        for stops."""
        if kind == "stop":
            return {"stop"}
        own = {L for L, c in lengths[(gk, pos)].items() if node in c} if node not in inferred_nodes else set()
        if own:
            return own
        older = [x for x in lineage if x in nodes and x not in inferred_nodes
                 and lineage.index(x) < lineage.index(node)] if node in lineage else []
        for x in reversed(older):
            got = {L for L, c in lengths[(gk, pos)].items() if x in c}
            if got:
                return got
        return set()

    def event_carriers(kind, gk, pos, L):
        if kind == "stop":
            return set(stop_carriers[(gk, pos)])
        return set().union(*[c for k, c in lengths[(gk, pos)].items() if k in L]) if L else set()

    out = []
    for r in rows:
        g = r["geneID"]
        # Phytozome IDs carry ".p" in the summaries but not in alignment names or detail tables (Brassicaceae)
        gk = g[:-2] if g.endswith(".p") else g
        paths = glob.glob(os.path.join(a.aln_dir, f"{F}_{gk}_*.aln"))
        enablers0 = set((r["Enablers_Nodes"] or "NA").split("|")) - {"NA"}
        if not paths and not enablers0:
            # no alignment and no enabler nodes, so no evidence to annotate (2 Brassicaceae PDNGs)
            rec = {"geneID": g, "Gene_type": r["Gene_type"], "Enablers_Nodes": r["Enablers_Nodes"],
                   "indel_event_flags": "NA", "stop_event_flags": "NA", "evidence_flagged": "NA"}
            if a.stops_detail:
                rec.update({k: "NA" for k in ("extant_support", "extant_contradicting",
                                              "enabler_nodes_without_extant_support",
                                              "enabler_nodes_contradicted_no_support")})
            out.append(rec)
            print(f"  {g}: no alignment, no enabler nodes; written as NA", file=sys.stderr)
            continue
        if len(paths) != 1:
            die(f"{len(paths)} alignments for {g} in {a.aln_dir}")
        aln = read_aln(paths[0])
        ref = aln[F]
        cds = ref.replace("-", "")
        cols = [i for i, c in enumerate(ref) if c != "-"]
        enablers = set((r["Enablers_Nodes"] or "NA").split("|")) - {"NA"}
        rec = {"geneID": g, "Gene_type": r["Gene_type"], "Enablers_Nodes": r["Enablers_Nodes"]}
        support = collections.defaultdict(list)
        via = {}
        extant = collections.defaultdict(set)
        contra = collections.defaultdict(set)
        for kind, pcol, ncol in (("indel", "Indels_event_RefPos", "Indels_Node"),
                                 ("stop", "Stops_event_RefPos", "Stops_Node")):
            pos = r.get(pcol, "NA")
            nodes = [n for n in (r[ncol] or "NA").split("|") if n != "NA"]
            if pos in ("", "NA") or not nodes:
                rec[f"{kind}_event_flags"] = "NA"
                continue
            p = int(pos)
            fl = []
            w = cds[max(0, p - 5):p + 4]
            # homopolymers make indel calls unreliable, not single-base substitutions (stops)
            if kind == "indel" and any(w[j:j + 4] == b * 4 for b in "ACGT" for j in range(len(w) - 3)):
                fl.append("homopolymer")
            if p <= 9 or p >= len(cds) - 9:
                fl.append("cds_edge")
            c0 = cols[min(p - 1, len(cols) - 1)]
            if any(n in aln and "N" in aln[n][max(0, c0 - 5):c0 + 6] for n in nodes):
                fl.append("N_nearby")
            if kind == "indel":
                anc_len = {L for L, carriers in lengths[(gk, pos)].items() if any(c in nodes for c in carriers)}
                if len(anc_len) > 1:
                    fl.append("multi_length")
            rec[f"{kind}_event_flags"] = ",".join(fl) or "ok"
            for n in nodes:
                if n in enablers:
                    support[n].append(not fl)
            if a.stops_detail:
                lo_c, hi_c = max(0, c0 - 10), c0 + 11
                for n in nodes:
                    if n in enablers and n in lineage:
                        # exact match: same position and, for indels, the length of THIS node's copy (NEE fix 1)
                        carriers = event_carriers(kind, gk, pos, node_lengths(kind, gk, pos, n, nodes))
                        side = side_branch(n) & species
                        extant[n] |= carriers & side
                        with_data = {s for s in side if s in aln and any(ch not in "-N" for ch in aln[s][lo_c:hi_c])}
                        contra[n] |= with_data - carriers
        if a.all_indel_events and enablers:
            orf = r["ORF_Node"]
            oclade = set(cfg.NODE_CLADES[orf]) if orf in cfg.NODE_CLADES else {F}
            inside_nodes = {n for n, c in cfg.NODE_CLADES.items() if set(c) <= oclade}
            inferred = collections.defaultdict(set)
            for item in (r.get("Enablers_living_inferred") or "NA").split(";"):
                if item != "NA":
                    x, k, q = item.split(":")
                    inferred[(k, q)].add(x)
            extant = collections.defaultdict(set)
            proxy = collections.defaultdict(dict)  # option_b_outer_side_branch: node -> {older branch node: carriers}
            outer = cfg.opt("option_b_outer_side_branch")
            for kind, pos, ev_nodes in all_events[gk]:
                nodes_e = ev_nodes | inferred.get((kind, pos), set())
                if nodes_e & inside_nodes:
                    continue
                for n in enablers & nodes_e:
                    if n not in lineage:
                        continue
                    L = node_lengths(kind, gk, pos, n, nodes_e, inferred.get((kind, pos), set()))
                    car = event_carriers(kind, gk, pos, L)
                    if car & (oclade - {F}):
                        continue
                    hit = car & side_branch(n) & species
                    if hit:
                        extant[n] |= hit
                    elif outer:
                        # the node's own side branch has no aligned data at the event: a living carrier on an older
                        # side branch of the focal lineage may stand in (T2T 2026-11-13, rule 2)
                        c0 = cols[min(int(pos) - 1, len(cols) - 1)]
                        lo_c, hi_c = max(0, c0 - 10), c0 + 11
                        side_n = side_branch(n) & species
                        if any(sp in aln and any(ch not in "-N" for ch in aln[sp][lo_c:hi_c]) for sp in side_n):
                            continue
                        for x in lineage[:lineage.index(n)]:
                            h = car & side_branch(x) & species
                            if h:
                                proxy[n].setdefault(x, set()).update(h)
            via = {}
            if outer:
                # each enabler node needs its own observed side branch, or a distinct older one standing in; one
                # outgroup branch cannot support both nodes (a lone outgroup cannot polarize the root)
                ens = [n for n in enablers if n in lineage]
                choices = [[n] if extant[n] else sorted(proxy[n]) for n in ens]
                best = None
                import itertools
                for combo in itertools.product(*choices) if all(choices) else ():
                    if len(set(combo)) == len(combo):
                        best = combo
                        break
                if best is not None:
                    for n, x in zip(ens, best):
                        if not extant[n]:
                            extant[n] |= proxy[n][x]
                            via[n] = x
            for n in enablers:
                if n in lineage and inferred and any(n in v for v in inferred.values()) and not extant[n]:
                    die(f"{g}: living-inferred enabler {n} has no living carrier; step 7 and step 7d disagree")
        rec["evidence_flagged"] = ("NA" if not enablers else
                                   "yes" if any(not any(v) for n, v in support.items() if n in enablers)
                                   or any(n not in support for n in enablers) else "no")
        if a.stops_detail:
            ordered = sorted(enablers, key=lambda x: (lineage.index(x) if x in lineage else -1, x))
            rec["extant_support"] = ";".join(f"{n}:{'via ' + via[n] + ':' if n in via else ''}{','.join(sorted(extant[n])) or '-'}"
                                             for n in sorted(enablers, key=lambda x: lineage.index(x) if x in lineage else -1)) or "NA"
            rec["extant_contradicting"] = ";".join(f"{n}:{','.join(sorted(contra[n] - extant[n])) or '-'}"
                                                   for n in sorted(enablers, key=lambda x: lineage.index(x) if x in lineage else -1)) or "NA"
            rec["enabler_nodes_without_extant_support"] = (
                "NA" if not enablers else ",".join(n for n in ordered if not extant[n]) or "none")
            rec["enabler_nodes_contradicted_no_support"] = (
                "NA" if not enablers else ",".join(n for n in ordered if not extant[n] and contra[n]) or "none")
        out.append(rec)

    path = os.path.join(a.nodes_dir, f"{a.out_prefix}_evidence.tsv")
    cols_out = ["geneID", "Gene_type", "Enablers_Nodes", "indel_event_flags", "stop_event_flags", "evidence_flagged"]
    if a.stops_detail:
        cols_out += ["extant_support", "extant_contradicting", "enabler_nodes_without_extant_support",
                     "enabler_nodes_contradicted_no_support"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols_out, delimiter="\t")
        w.writeheader()
        w.writerows(out)
    c = collections.Counter((o["Gene_type"], o["evidence_flagged"]) for o in out)
    print(f"{os.path.basename(path)}: {len(out)} loci; evidence_flagged by type {dict(c)}", file=sys.stderr)


if __name__ == "__main__":
    main()
