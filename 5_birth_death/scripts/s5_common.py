"""
s5_common.py

Shared helpers for LINGUA stage 5: the species tree read from the stage 4 tree configuration
(tree_config.json, written by 4_enablers/scripts/make_tree_config.py), branch lengths, species-name
aliases, gene-ID normalization and small TSV helpers. Not run directly.
"""
import csv
import json
import os
import re
import sys

csv.field_size_limit(10 ** 9)


def die(msg):
    sys.exit(f"ERROR: {msg}")


def need(*paths):
    for p in paths:
        if not os.path.exists(p):
            die(f"input not found: {p}")


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def strip_isoform(gid):
    """Base locus ID: a trailing '.p' (or '.p2'...) and then one trailing '.<number>' are removed. Anchored at the
    end so IDs with an internal period (Bostr.23794s0730.1.p) keep their scaffold number; an unanchored version
    merged different Boechera and Capsella genes into one key. In ID schemes where the trailing number is part of
    the gene name (evm.model.Chr1.12) different genes would share a key; check_id_collisions stops the run then."""
    return re.sub(r"\.\d+$", "", strip_p(gid))


def check_id_collisions(ids_by_species, what):
    """Stop when two different IDs of one species share a strip_isoform key (they would be joined as one gene).
    ids_by_species: species -> iterable of raw IDs."""
    bad = []
    for sp, ids in ids_by_species.items():
        seen = {}
        for g in set(ids):
            seen.setdefault(strip_isoform(g), set()).add(g)
        col = sorted(sorted(v) for v in seen.values() if len(v) > 1)
        if col:
            bad.append(f"{sp}: {len(col)} keys shared by different IDs (e.g. {col[:2]})")
    if bad:
        die(f"{what}: removing isoform suffixes ('.1', '.p') would join different IDs. Either several isoforms of "
            "one gene are present (keep one per gene, e.g. the longest), or the trailing number is part of the gene "
            "name in this ID scheme (rename the IDs so that an isoform suffix is the only trailing number):\n  " + "\n  ".join(bad))


def strip_p(gid):
    return re.sub(r"\.p\d*$", "", gid)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(["" if v is None else v for v in (r if isinstance(r, (list, tuple)) else [r[h] for h in header])])


def loose_id(gid, species=None):
    """An ID with its formatting removed (case, 'db|' prefixes, a leading species tag, '_' vs '.', version
    suffixes), used only to detect ID-form mismatches between files, never to join."""
    g = gid.strip().split("|")[-2 if gid.count("|") >= 2 else -1].lower()
    if species:
        for tag in (species.lower() + "_", species.lower().split("_")[0] + "_"):
            if g.startswith(tag):
                g = g[len(tag):]
    g = re.sub(r"[._]p\d*$", "", g)
    g = re.sub(r"[._](t)?\d+$", "", g)
    return g.replace("_", ".")


def check_hog_ids(tree, genes_by_species, hog_ids_by_species, key):
    """Stop when gene IDs and HOG-table IDs differ in form. genes_by_species: species -> set of geneIDs of the
    classified loci; hog_ids_by_species: species -> set of raw HOG-table IDs; key(species, id) -> the join key.
    A species stops the run when it has HOG-table IDs but none of its genes match, or when some of its genes match
    a HOG-table ID only after formatting is removed. Writes one line per species to stderr."""
    check_id_collisions(hog_ids_by_species, "HOG table")
    check_id_collisions(genes_by_species, "classified loci")
    bad = []
    log("species\tgenes\tin_HOG_table\tmatch_only_if_reformatted")
    for sp in tree.species:
        genes = genes_by_species.get(sp, set())
        hog = hog_ids_by_species.get(sp, set())
        if not genes:
            continue
        hk = {key(sp, h) for h in hog}
        hit = {g for g in genes if key(sp, g) in hk}
        loose = {loose_id(h, sp) for h in hog}
        near = sorted(g for g in genes - hit if loose_id(g, sp) in loose)
        log(f"{sp}\t{len(genes)}\t{len(hit)} ({len(hit) / len(genes):.0%})\t{len(near)}")
        if hog and not hit:
            bad.append(f"{sp}: none of {len(genes)} genes is in the HOG table (e.g. {sorted(genes)[:2]} vs "
                       f"{sorted(hog)[:2]})")
        elif near:
            bad.append(f"{sp}: {len(near)} genes match a HOG-table ID only after reformatting (e.g. {near[:3]})")
    if bad:
        die("gene IDs and HOG-table IDs differ in form; make them consistent before stage 5:\n  " + "\n  ".join(bad))


class Tree:
    """The species tree of tree_config.json: species (leaves), internal nodes with their clades, parent and
    children of every node, node age order, the nodes at which an ORF may originate, and the ingroup crown."""

    def __init__(self, config_path):
        need(config_path)
        c = json.load(open(config_path))
        self.species = list(c["species"])
        self.aliases = dict(c["aliases"])
        self.clades = {k: frozenset(v) for k, v in c["node_clades"].items()}
        self.young_to_old = list(c["node_order_young_to_old"])
        self.excluded_orf = set(c["excluded_orf_nodes"])
        self.crown = c["ingroup_crown"]["label"]
        sp = set(self.species)
        if set(self.young_to_old) != set(self.clades):
            die("node_order_young_to_old and node_clades list different nodes")
        roots = [n for n, cl in self.clades.items() if cl == sp]
        if len(roots) != 1:
            die(f"expected one node whose clade is every species, found {roots}")
        self.root = roots[0]
        self.parent = {}
        for x in list(self.clades) + self.species:
            members = self.clades.get(x, frozenset([x]))
            cands = [n for n, cl in self.clades.items() if n != x and members < cl]
            if not cands:
                continue
            self.parent[x] = min(cands, key=lambda n: len(self.clades[n]))
        self.children = {}
        for ch, p in self.parent.items():
            self.children.setdefault(p, []).append(ch)
        for n in self.clades:
            leaves = {s for s in self.subtree(n) if s in sp}
            if leaves != set(self.clades[n]):
                die(f"node {n}: clade in tree_config does not match the nesting of clades")
        self.depth = {}
        for x in list(self.clades) + self.species:
            d, y = 0, x
            while y in self.parent:
                y, d = self.parent[y], d + 1
            self.depth[x] = d
        old_to_young = list(reversed(self.young_to_old))
        self.rank = {n: i for i, n in enumerate(old_to_young)}
        self.origin_nodes = [n for n in old_to_young if n not in self.excluded_orf]
        self.ingroup = set(self.clades[self.crown])

    def canon(self, name):
        s = str(name).strip()
        return self.aliases.get(s, self.aliases.get(s.lower(), s))

    def subtree(self, n):
        out, stack = [], [n]
        while stack:
            x = stack.pop()
            out.append(x)
            stack.extend(self.children.get(x, []))
        return out

    def postorder(self, n):
        out = []
        for c in self.children.get(n, []):
            out.extend(self.postorder(c))
        out.append(n)
        return out

    def oldest(self, nodes):
        """Oldest of a set of internal nodes: fewest ancestors, ties by the node order of tree_config."""
        return min(nodes, key=lambda n: (self.depth[n], self.rank[n]))

    def branch_of(self, orf_node):
        """The branch on which a locus with this ORF node originated: an internal node at which an ORF may
        originate, or a species for a species-specific locus ('<species>-specific' or the species name).
        None for loci outside the window (ORF node excluded in tree_config)."""
        n = str(orf_node)
        if n in self.clades:
            return None if n in self.excluded_orf else n
        if n.endswith("-specific"):
            n = n[: -len("-specific")]
        s = self.canon(n)
        return s if s in self.species else None

    def report_branches(self):
        """Branches written to the rate table: the ingroup crown, every internal node below it, then every
        species (outgroups included, as in the published table)."""
        below = [n for n in self.subtree(self.crown) if n in self.clades]
        return sorted(below, key=lambda n: self.rank[n]) + list(self.species)


def read_branch_lengths(path, tree):
    """Branch lengths (time units of the user's choice, e.g. million years) from a TSV with columns
    branch, length; or node ages from a TSV with columns node, age (species age 0 unless listed)."""
    need(path)
    rows = read_tsv(path)
    cols = set(rows[0]) if rows else set()
    out = {}
    if {"branch", "length"} <= cols:
        for r in rows:
            b = r["branch"] if r["branch"] in tree.clades else tree.canon(r["branch"])
            if b not in tree.parent and b != tree.root:
                die(f"branch lengths: unknown branch {r['branch']}")
            v = r["length"].strip()
            out[b] = None if v in ("", "NA") else float(v)
    elif {"node", "age"} <= cols:
        age = {s: 0.0 for s in tree.species}
        for r in rows:
            n = r["node"] if r["node"] in tree.clades else tree.canon(r["node"])
            age[n] = float(r["age"])
        for x, p in tree.parent.items():
            if x in age and p in age:
                out[x] = age[p] - age[x]
                if out[x] < 0:
                    die(f"node ages: {x} is older than its parent {p}")
    else:
        die(f"{path}: need columns 'branch, length' or 'node, age'")
    return out
