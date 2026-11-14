#!/usr/bin/env python3
"""
make_tree_config.py

Build the tree configuration that every LINGUA script reads, from a rooted Newick tree
with named tips and named internal nodes (for example the Cactus labels fullTreeAnc238,
PrimatesAnc007).

Internal nodes are relabeled N0, N1, ... in breadth-first order from the root, children
in Newick order, so that a parent always carries a lower number than its descendants,
which the LINGUA scripts assume when they take the lowest-numbered node as the oldest.

Excluded nodes follow the rule used for the Brassicaceae analysis, stated relative to the
crown node of the ingroup (the densely sampled clade in which de novo genes are sought).
    ORF nodes excluded:     the ingroup crown, its ancestors, and every node outside it
    enabler nodes excluded: ancestors of the crown except its parent
On the Brassicaceae tree with crown N3 this gives ORF exclusions N0-N3 and enabler
exclusions N0-N1, the values hardcoded in the published scripts.
--allow-orf-node removes named nodes from the ORF exclusions (for the T2T primates the
crown N4 is allowed, because rodents and laurasiatherians polarize it; decided 2026-09-24).
--conditional-enabler-node NODE=REQ[,REQ] removes NODE from the enabler exclusions but lets it
count only for a disabler that every REQ node also carries (T2T: the root N0 only together with
N1, since the root alone is unpolarized; decided 2026-10-08).

Outputs (in --outdir):
    tree_config.json          species, aliases, node order, clades, exclusions
    species_tree_N.nwk        the tree with internal nodes relabeled N0..Nk
    node_name_map.tsv         N label, original internal-node name, descendant species

Usage:
    python3 make_tree_config.py TREE.nwk --ingroup-crown PrimatesAnc007 --outdir data \
        [--allow-orf-node PrimatesAnc007] [--conditional-enabler-node fullTreeAnc238=fullTreeAnc235]
"""

import argparse
import json
import os
import re
import sys
from collections import deque


class Node:
    def __init__(self, name=None, length=None):
        self.name, self.length, self.children, self.parent = name, length, [], None

    def is_tip(self):
        return not self.children


def parse_newick(text):
    text = text.strip().rstrip(";")
    pos = 0

    def label():
        nonlocal pos
        m = re.match(r"([^,():;]*)(?::([0-9.eE+-]+))?", text[pos:])
        pos += m.end()
        return (m.group(1).strip() or None), m.group(2)

    def node():
        nonlocal pos
        n = Node()
        if text[pos] == "(":
            pos += 1
            while True:
                c = node()
                c.parent = n
                n.children.append(c)
                if text[pos] == ",":
                    pos += 1
                elif text[pos] == ")":
                    pos += 1
                    break
                else:
                    sys.exit(f"ERROR: unexpected {text[pos]!r} at position {pos}")
        n.name, n.length = label()
        return n

    root = node()
    if pos != len(text):
        sys.exit(f"ERROR: trailing text after position {pos}")
    return root


def tips(n):
    return [n.name] if n.is_tip() else [t for c in n.children for t in tips(c)]


def ancestors(n):
    out = []
    while n.parent is not None:
        n = n.parent
        out.append(n)
    return out


def to_newick(n):
    s = n.name or ""
    if n.children:
        s = "(" + ",".join(to_newick(c) for c in n.children) + ")" + s
    if n.length is not None:
        s += f":{n.length}"
    return s


def abbreviation(species):
    parts = species.split("_")
    return parts[0][0] + "".join(parts[1:]) if len(parts) > 1 else species


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("tree")
    ap.add_argument("--ingroup-crown", required=True,
                    help="Original name of the ingroup crown node (for example PrimatesAnc007).")
    ap.add_argument("--outdir", default="data")
    ap.add_argument("--allow-orf-node", action="append", default=[],
                    help="Node (original or N name) to keep out of the ORF exclusions; repeatable.")
    ap.add_argument("--conditional-enabler-node", action="append", default=[],
                    help="NODE=REQ[,REQ]: NODE counts as an enabler node only when every REQ node "
                         "carries the same disabler; repeatable.")
    ap.add_argument("--stops-focal-frame-only", action="store_true",
                    help="Step 6: count only stop events in the focal reading frame, excluding the focal "
                         "terminal stop (decided 2026-10-09 for the T2T tree).")
    ap.add_argument("--reference-species", default=None,
                    help="Step 8: species whose relatives are listed first in the summary tables "
                         "(default: the first ingroup species).")
    args = ap.parse_args()
    if not os.path.isfile(args.tree):
        sys.exit(f"ERROR: tree not found: {args.tree}")

    root = parse_newick(open(args.tree).read())
    internal, q = [], deque([root])
    while q:
        n = q.popleft()
        if not n.is_tip():
            internal.append(n)
            q.extend(n.children)
    orig = {}
    for i, n in enumerate(internal):
        orig[f"N{i}"] = n.name
        n.name = f"N{i}"

    names = [n.name for n in internal]
    crown = [n for n in internal if orig[n.name] == args.ingroup_crown or n.name == args.ingroup_crown]
    if len(crown) != 1:
        sys.exit(f"ERROR: ingroup crown {args.ingroup_crown!r} not found among internal nodes")
    crown = crown[0]
    inside = {n.name for n in internal if crown in ancestors(n) or n is crown}
    excluded_orf = sorted({n for n in names if n not in inside} | {crown.name},
                          key=lambda x: int(x[1:]))
    for a in args.allow_orf_node:
        hit = [n.name for n in internal if n.name == a or orig[n.name] == a]
        if len(hit) != 1:
            sys.exit(f"ERROR: --allow-orf-node {a!r} not found among internal nodes")
        if hit[0] not in excluded_orf:
            sys.exit(f"ERROR: --allow-orf-node {a!r} ({hit[0]}) is not an excluded ORF node")
        excluded_orf.remove(hit[0])
    excluded_enabler = sorted({a.name for a in ancestors(crown)} - {crown.parent.name}
                              if crown.parent else set(), key=lambda x: int(x[1:]))

    def resolve(a):
        hit = [n.name for n in internal if n.name == a or orig[n.name] == a]
        if len(hit) != 1:
            sys.exit(f"ERROR: node {a!r} not found among internal nodes")
        return hit[0]
    conditional = {}
    for spec in args.conditional_enabler_node:
        if "=" not in spec:
            sys.exit(f"ERROR: --conditional-enabler-node {spec!r} is not NODE=REQ[,REQ]")
        node, reqs = spec.split("=", 1)
        node = resolve(node)
        if node not in excluded_enabler:
            sys.exit(f"ERROR: {node} is not an excluded enabler node")
        excluded_enabler.remove(node)
        conditional[node] = [resolve(r) for r in reqs.split(",") if r]

    species = tips(root)
    if len(set(species)) != len(species):
        sys.exit("ERROR: duplicate tip names")
    aliases = {}
    for sp in species:
        ab = abbreviation(sp)
        for a in (sp, sp.lower(), ab, ab.lower()):
            if a in aliases and aliases[a] != sp:
                sys.exit(f"ERROR: alias {a!r} maps to both {aliases[a]} and {sp}")
            aliases[a] = sp
    clades = {n.name: tips(n) for n in internal}

    cfg = {
        "source_tree": os.path.abspath(args.tree),
        "ingroup_crown": {"label": crown.name, "original": orig[crown.name]},
        "species": species,
        "ingroup_species": tips(crown),
        "aliases": aliases,
        "node_original_names": orig,
        "node_order_young_to_old": list(reversed(names)),
        "node_clades": clades,
        "excluded_orf_nodes": excluded_orf,
        "excluded_enabler_nodes": excluded_enabler,
        "conditional_enabler_nodes": conditional,
        "stops_focal_frame_only": args.stops_focal_frame_only,
    }
    if args.reference_species:
        cfg["reference_species"] = args.reference_species
    os.makedirs(args.outdir, exist_ok=True)
    json.dump(cfg, open(os.path.join(args.outdir, "tree_config.json"), "w"), indent=2)
    open(os.path.join(args.outdir, "species_tree_N.nwk"), "w").write(to_newick(root) + ";\n")
    with open(os.path.join(args.outdir, "node_name_map.tsv"), "w") as fh:
        fh.write("node\toriginal_name\tn_species\tdescendant_species\n")
        for n in internal:
            fh.write(f"{n.name}\t{orig[n.name]}\t{len(clades[n.name])}\t{','.join(clades[n.name])}\n")
    print(f"{len(species)} species, {len(internal)} internal nodes N0-N{len(internal) - 1}; "
          f"ingroup crown {crown.name} ({orig[crown.name]}), {len(tips(crown))} ingroup species",
          file=sys.stderr)
    print(f"excluded ORF nodes {excluded_orf}; excluded enabler nodes {excluded_enabler}; "
          f"conditional enabler nodes {conditional}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
