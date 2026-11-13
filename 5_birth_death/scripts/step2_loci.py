"""
step2_loci.py

Stage 5 step 2. Groups the classified genes of different species that are the same ancestral locus, and records in
which genomes each gene's ORF is found.

The same ancestral ORF is often annotated, and classified, separately in several descendant species. For every gene G
of species S whose ORF node N is an internal node at which ORFs may originate, the coordinates of G's locus in each
other species T of clade(N) (from step 1) are compared with the annotated loci of T; G is joined to the first locus of T
that overlaps them. Joined genes form one locus (union-find, so joins are transitive). A locus takes the OLDEST ORF node
among its members, and is a DNG when any member is a DNG, otherwise a PDNG. Species-specific genes never enter the
lookup, so they cannot pull an ancestral locus to their own node or type; each stays a locus of its own.

A gene's presence set is every genome with a step 1 row for it, namely genomes in which its locus aligns and tBLASTn
finds its ORF. Genes are keyed by (species, geneID), so the same ID in two species stays two genes.

This is ancestral_nodes_estimates_and_trees/anc_nodes_estimates.py of the Brassicaceae analysis with the tree read from
tree_config.json.

Usage:
  python step2_loci.py --tree-config tree_config.json --coords coords.tsv --outdir DIR
Output (in DIR): loci.tsv (ClusterID, Node, Gene_type, one column per species with that species' member geneIDs joined
        by ';'), genes.tsv (species, geneID, ORF_Node, Gene_type, ClusterID, Cluster_node, present_in = genomes with
        tBLASTn evidence joined by ';')
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s5_common import Tree, die, log, need, read_tsv, write_tsv  # noqa: E402


class UnionFind:
    def __init__(self):
        self.parent, self.rank = {}, {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x], self.rank[x] = x, 0
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree-config", required=True)
    ap.add_argument("--coords", required=True, help="step 1 output")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    need(a.tree_config, a.coords)
    os.makedirs(a.outdir, exist_ok=True)
    for f in ("loci.tsv", "genes.tsv"):
        if os.path.exists(os.path.join(a.outdir, f)):
            die(f"output exists: {os.path.join(a.outdir, f)}")
    tree = Tree(a.tree_config)
    anc = set(tree.origin_nodes)

    rows = read_tsv(a.coords)
    log(f"{len(rows):,} coordinate rows")
    by_focal = defaultdict(list)
    for r in rows:
        r["focal_species"] = tree.canon(r["focal_species"])
        r["species_name"] = tree.canon(r["species_name"])
        for k in ("focal_species", "species_name"):
            if r[k] not in tree.species:
                die(f"{a.coords}: {k} {r[k]} is not a species of tree_config")
        r["start"], r["end"] = int(r["start"]), int(r["end"])
        by_focal[r["focal_species"]].append(r)

    meta = {}                       # (species, gene) -> (ORF_Node, Gene_type), from the first own-genome row
    index = {}                      # species -> chrom -> sorted [(start, end, gene)] of its ancestral-node genes
    presence = defaultdict(set)     # (species, gene) -> genomes with a row
    for sp, rs in by_focal.items():
        idx = defaultdict(list)
        for r in rs:
            presence[(sp, r["geneID"])].add(r["species_name"])
            if r["species_name"] != sp:
                continue
            key = (sp, r["geneID"])
            if key not in meta:
                meta[key] = (r["ORF_Node"], r["Gene_type"])
            if r["ORF_Node"] in anc:
                idx[r["chromosome"].lower()].append((r["start"], r["end"], r["geneID"]))
        for c in idx:
            idx[c].sort()
        index[sp] = idx
    no_self = {k for k in presence if k not in meta}
    if no_self:
        log(f"  {len(no_self):,} genes have no row in their own genome and are left out (e.g. {sorted(no_self)[:3]})")

    def first_overlap(sp, chrom, s, e):
        for gs, ge, g in index.get(sp, {}).get(chrom.lower(), []):
            if gs > e:
                break
            if ge >= s:
                return g
        return None

    uf = UnionFind()
    for sp, rs in by_focal.items():
        rows_of = defaultdict(list)
        for r in rs:
            rows_of[r["geneID"]].append(r)
        for r in rs:
            if r["species_name"] != sp or r["ORF_Node"] not in anc:
                continue
            g = r["geneID"]
            for other in tree.clades[r["ORF_Node"]]:
                if other == sp or other not in index:
                    continue
                for o in rows_of[g]:
                    if o["species_name"] != other:
                        continue
                    hit = first_overlap(other, o["chromosome"], o["start"], o["end"])
                    if hit is not None:
                        uf.union((sp, g), (other, hit))

    clusters = defaultdict(list)
    for key, (node, gt) in meta.items():
        clusters[uf.find(key)].append((key, node, gt))

    def old_node(members):
        nodes = {m[1] for m in members if m[1] in anc}
        return tree.oldest(nodes) if nodes else members[0][1]

    def sort_key(root):
        m = clusters[root]
        return (tree.rank.get(old_node(m), 10 ** 6), min(k[1] for k, _, _ in m))

    loci, genes, cid_of = [], [], {}
    for i, root in enumerate(sorted(clusters, key=sort_key), 1):
        m = clusters[root]
        node = old_node(m)
        gt = "DNG" if any(x[2] == "DNG" for x in m) else "PDNG"
        cid = f"C{i}"
        cells = defaultdict(list)
        for (sp, g), _, _ in m:
            cells[sp].append(g)
            cid_of[(sp, g)] = (cid, node)
        loci.append([cid, node, gt] + [";".join(cells.get(sp, [])) for sp in tree.species])
    for (sp, g), (node, gt) in sorted(meta.items()):
        cid, cnode = cid_of[(sp, g)]
        genes.append([sp, g, node, gt, cid, cnode, ";".join(s for s in tree.species if s in presence[(sp, g)])])

    reassigned = sum(1 for (sp, g), (node, _) in meta.items() if cid_of[(sp, g)][1] != node)
    assert sum(len(v) for v in clusters.values()) == len(meta) == len(genes)
    write_tsv(os.path.join(a.outdir, "loci.tsv"), ["ClusterID", "Node", "Gene_type"] + tree.species, loci)
    write_tsv(os.path.join(a.outdir, "genes.tsv"),
              ["species", "geneID", "ORF_Node", "Gene_type", "ClusterID", "Cluster_node", "present_in"], genes)
    log(f"{len(meta):,} genes -> {len(loci):,} loci; {reassigned:,} genes moved to an older node by their locus")


if __name__ == "__main__":
    main()
