#!/usr/bin/env python3
"""
db_report.py - OPTIONAL: check how well a homology database covers the tree of life outside your analysis clade.

Reads the species of a database (a species-tagged FASTA, e.g. from make-db, or a plain list of names), places each species
on the NCBI taxonomy, and sorts it into a TIER by how far it lies outside your clade (the closest outside group first: for
a family, its order; then its class; and so on out to the root). It reports how many species and sequences each tier holds,
which tiers are empty or thin, species that hold a very large share of the database, species that look too small, and any
species that sit INSIDE your analysis clade (these do not belong in a database that tests for outside homologs).

It reads local files only: the NCBI taxonomy dump downloaded once with setup/get_taxonomy.sh, or (with --tiers) a table you
write yourself, in which case no taxonomy is needed. It changes and downloads nothing. This is guidance, not a verdict:
public databases are uneven, and tier sizes are heuristics; see docs/databases.md.

  python3 db_report.py --db db/mosses.faa db/fungi.faa --clade Brassicaceae --taxonomy taxonomy/ -o report/
  python3 db_report.py --db db/all.faa --clade-inputs candidates/ --taxonomy taxonomy/ -o report/ --tree report/tree.nwk
  python3 db_report.py --species-list names.txt --tiers my_tiers.tsv -o report/     (no taxonomy needed)
"""

import argparse
import csv
import json
import os
import re
import sys
from array import array

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402

MODULE = "db-report"
NAME_CLASSES = ("scientific name", "synonym", "equivalent name")   # 'authority', 'type material' etc. are not used for matching
DEFAULT_RANKS = "family,order,class,phylum,kingdom,domain"
INSIDE = "INSIDE_CLADE"
UNMATCHED = "UNMATCHED"
SPECIES_TAG_RE = h.SPECIES_TAG_RE


# --------------------------------------------------------------------------
# The database's species
# --------------------------------------------------------------------------
looks_like_species = h.looks_like_species


def species_from_fasta(paths):
    """Count sequences per '__Species_Name' tag over the FASTA files. Returns (counts, untagged, odd_tags)."""
    counts, untagged, odd = {}, 0, 0
    for path in paths:
        if not os.path.isfile(path):
            raise h.HepError("database not found: %s" % path)
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                rid = line[1:].split(None, 1)[0] if line[1:].strip() else ""
                tag = h.split_id(rid)[1]
                if tag is None:
                    untagged += 1
                elif looks_like_species(tag):
                    counts[tag] = counts.get(tag, 0) + 1
                else:
                    odd += 1       # an '__' inside a random-looking ID: not read as a species
                    untagged += 1
    return counts, untagged, odd


def species_from_list(path):
    """One species per line, optionally followed by a tab and a sequence count. '#' lines are skipped."""
    if not os.path.isfile(path):
        raise h.HepError("species list not found: %s" % path)
    counts = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            name = re.sub(r"\s+", "_", parts[0].strip())
            n = None
            if len(parts) > 1 and parts[1].strip():
                try:
                    n = int(parts[1].replace(",", ""))
                except ValueError:
                    raise h.HepError("%s: the count for %s is not a number: %r" % (path, name, parts[1]))
            counts[name] = counts.get(name, 0) + (n or 0) if n is not None else counts.get(name)
    return counts


# --------------------------------------------------------------------------
# The taxonomy (NCBI names.dmp and nodes.dmp)
# --------------------------------------------------------------------------
def norm_name(name):
    return re.sub(r"\s+", " ", str(name).replace("_", " ").strip()).lower()


class Taxonomy(object):
    """Parent, rank and the names asked for, read from the dump. Fields are separated by TAB-pipe-TAB."""

    def __init__(self, folder):
        self.folder = folder
        names_p, nodes_p = os.path.join(folder, "names.dmp"), os.path.join(folder, "nodes.dmp")
        for p in (names_p, nodes_p):
            if not os.path.isfile(p):
                raise h.HepError("%s not found. Download the taxonomy once with: bash setup/get_taxonomy.sh %s" % (p, folder))
        self.names_p, self.nodes_p = names_p, nodes_p
        self.parent = array("i")
        self.rank_code = array("B")
        self.ranks = []           # code -> rank name
        self._rank_index = {}
        with open(nodes_p, "r", encoding="utf-8") as fh:
            for line in fh:
                f = line.split("\t|\t")
                if len(f) < 3:
                    raise h.HepError("%s: a line does not look like an NCBI nodes.dmp line: %r" % (nodes_p, line[:60]))
                tid, par, rank = int(f[0]), int(f[1]), f[2]
                if tid >= len(self.parent):
                    grow = tid + 1 - len(self.parent) + 500000
                    self.parent.extend([0] * grow)
                    self.rank_code.extend([0] * grow)
                if rank not in self._rank_index:
                    self._rank_index[rank] = len(self.ranks)
                    self.ranks.append(rank)
                self.parent[tid] = par
                self.rank_code[tid] = self._rank_index[rank]
        if not self.ranks:
            raise h.HepError("%s is empty" % nodes_p)
        self.source = self._source_note(folder)
        self._sci = {}            # taxid -> scientific name (only those asked for)

    @staticmethod
    def _source_note(folder):
        p = os.path.join(folder, "TAXONOMY_SOURCE.txt")
        info = {}
        if os.path.isfile(p):
            with open(p, "r", encoding="utf-8") as fh:
                for line in fh:
                    k, _, v = line.strip().partition("=")
                    info[k] = v
        return info or {"source": "unknown (no TAXONOMY_SOURCE.txt next to the dump)"}

    def has(self, tid):
        return 0 < tid < len(self.parent) and self.parent[tid] != 0

    def rank(self, tid):
        return self.ranks[self.rank_code[tid]]

    def lineage(self, tid):
        """taxid up to the root, inclusive."""
        out, seen = [tid], {tid}
        while True:
            par = self.parent[tid]
            if par == tid or par == 0 or par in seen:
                break
            out.append(par)
            seen.add(par)
            tid = par
            if len(out) > 400:
                raise h.HepError("the taxonomy has a loop or is malformed near taxid %d" % tid)
        return out

    def find_names(self, wanted):
        """{normalised name: [taxids]} for the wanted normalised names (scientific names, synonyms, equivalent names)."""
        wanted = set(wanted)
        found = {}
        with open(self.names_p, "r", encoding="utf-8") as fh:
            for line in fh:
                f = line.split("\t|\t")
                if len(f) < 4:
                    continue
                cls = f[3].rstrip("\t|\r\n")
                if cls not in NAME_CLASSES:
                    continue
                key = norm_name(f[1])
                if key in wanted:
                    tid = int(f[0])
                    if tid not in found.setdefault(key, []):
                        found[key].append(tid)
        return found

    def sci_names(self, tids):
        need = {t for t in tids if t not in self._sci}
        if need:
            with open(self.names_p, "r", encoding="utf-8") as fh:
                for line in fh:
                    f = line.split("\t|\t")
                    if len(f) >= 4 and f[3].rstrip("\t|\r\n") == "scientific name":
                        t = int(f[0])
                        if t in need:
                            self._sci[t] = f[1]
        return {t: self._sci.get(t, str(t)) for t in tids}


def resolve(tax, species, taxid_map):
    """
    Map each database species tag to a taxid. Tries the whole name, then the first two words (a strain or subspecies
    counts as its species), then the first word as a genus. A name that matches more than one taxon is ambiguous and is
    NOT guessed. Returns {species: (taxid or None, matched_as, level, note)}.
    """
    want = set()
    for sp in species:
        parts = norm_name(sp).split(" ")
        want.add(" ".join(parts))
        want.add(" ".join(parts[:2]))
        want.add(parts[0])
    found = tax.find_names(want)
    out = {}
    for sp in species:
        if sp in taxid_map:
            t = taxid_map[sp]
            out[sp] = (t, "given in --taxid-map", "given", None) if tax.has(t) else (None, None, "none", "taxid %d not in the taxonomy" % t)
            continue
        parts = norm_name(sp).split(" ")
        tries = [("exact", " ".join(parts)), ("species", " ".join(parts[:2])), ("genus", parts[0])]
        done = None
        seen = set()
        for level, key in tries:
            if not key or (level, key) in seen:
                continue
            seen.add((level, key))
            cands = found.get(key, [])
            if level == "species":
                sp_c = [t for t in cands if tax.rank(t) == "species"]
                cands = sp_c or cands
            if level == "genus":
                cands = [t for t in cands if tax.rank(t) == "genus"]
            if len(cands) == 1:
                done = (cands[0], key, level, None)
                break
            if len(cands) > 1:
                done = (None, key, "none", "ambiguous: '%s' matches %d taxa (%s); give the right one in --taxid-map"
                        % (key, len(cands), ", ".join(str(c) for c in cands[:4])))
                break
        out[sp] = done or (None, None, "none", "no taxon of this name in the taxonomy")
    return out


# --------------------------------------------------------------------------
# Tiers
# --------------------------------------------------------------------------
def lca(tax, tids):
    paths = [list(reversed(tax.lineage(t))) for t in tids]
    common = None
    for i in range(min(len(p) for p in paths)):
        if len({p[i] for p in paths}) == 1:
            common = paths[0][i]
        else:
            break
    return common


def build_ladder(tax, clade, ranks):
    """The clade's ancestors at the chosen ranks, nearest first, then the root as the last catch-all."""
    lin = tax.lineage(clade)[1:]
    ladder = [t for t in lin if ranks == "all" or tax.rank(t) in ranks]
    return ladder


def assign_tier(tax, tid, clade, ladder):
    """INSIDE_CLADE, or the nearest ladder ancestor whose subtree holds the species, or 'outside' (index len(ladder))."""
    lin = set(tax.lineage(tid))
    if clade in lin:
        return INSIDE, None
    for i, anc in enumerate(ladder):
        if anc in lin:
            return "ladder", i
    return "outside", len(ladder)


# --------------------------------------------------------------------------
# Tree (Newick) of the database species and the clade
# --------------------------------------------------------------------------
def _nw(name):
    return re.sub(r"[\s]+", "_", re.sub(r"[()\[\]:;,'\"]", "", str(name)))


def newick(tax, species_tids, clade, clade_name):
    """Taxonomy-induced tree over the matched species plus a marker leaf for the analysis clade. No branch lengths."""
    children, tip_names = {}, {}
    def add(path, leaf):
        for a, b in zip(path, path[1:]):
            children.setdefault(a, [])
            if b not in children[a]:
                children[a].append(b)
        children.setdefault(path[-1], [])
        tip_names[path[-1]] = leaf
    for sp, t in species_tids.items():
        add(list(reversed(tax.lineage(t))), sp)
    marker = -clade   # a synthetic child of the clade node
    add(list(reversed(tax.lineage(clade))) + [marker], "ANALYSIS_CLADE_" + _nw(clade_name))
    root = 1 if 1 in children else next(iter(children))
    names = tax.sci_names([t for t in children if t > 0])

    def collapse(node):
        kids = children.get(node, [])
        while len(kids) == 1 and node != clade and node != root and kids[0] != marker:
            node = kids[0]
            kids = children.get(node, [])
        return node

    def write(node):
        node = collapse(node) if node != root else node
        kids = children.get(node, [])
        if not kids:
            return tip_names.get(node, str(node))
        inner = ",".join(write(k) for k in kids)
        return "(%s)%s" % (inner, _nw(names.get(node, "")) if node > 0 else "")
    return write(root) + ";\n"


# --------------------------------------------------------------------------
def write_tsv(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for row in rows:
            w.writerow(["NA" if (c is None or c == "") else c for c in row])


def read_pairs(path, what):
    if not os.path.isfile(path):
        raise h.HepError("%s not found: %s" % (what, path))
    out = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 2 or not p[0].strip() or not p[1].strip():
                raise h.HepError("%s: expected 'name<TAB>value', got %r" % (path, line))
            out[re.sub(r"\s+", "_", p[0].strip())] = p[1].strip()
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="db-report", description=__doc__.strip().split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version="%s %s" % (h.PIPELINE_NAME, h.__version__))
    src = ap.add_argument_group("the database")
    src.add_argument("--db", nargs="+", metavar="FASTA", help="species-tagged FASTA file(s) (make-db writes these)")
    src.add_argument("--species-list", metavar="FILE", help="instead of FASTA: one species per line, optionally TAB and a sequence count")
    cl = ap.add_argument_group("your analysis clade (give one)")
    cl.add_argument("--clade", metavar="NAME", help="a taxon name, e.g. Brassicaceae")
    cl.add_argument("--clade-taxid", type=int, metavar="ID", help="an NCBI taxid (use if the name is ambiguous)")
    cl.add_argument("--clade-species", metavar="FILE", help="a list of your analysis species (one per line); the clade is their common ancestor")
    cl.add_argument("--clade-inputs", metavar="DIR", help="the folder of your analysis FASTA files; species come from the file names")
    tx = ap.add_argument_group("where tiers come from")
    tx.add_argument("--taxonomy", metavar="DIR", help="folder with names.dmp and nodes.dmp (setup/get_taxonomy.sh writes it)")
    tx.add_argument("--tiers", metavar="FILE", help="your own table: species<TAB>tier (no taxonomy needed; overrides it)")
    tx.add_argument("--expected-tiers", metavar="A,B,C", help="with --tiers: the tiers you expect, so empty ones are reported")
    tx.add_argument("--taxid-map", metavar="FILE", help="species<TAB>taxid, to fix names the taxonomy cannot match")
    tx.add_argument("--species-map", metavar="FILE", help="for --clade-inputs: two columns (species name, file basename), a header row, tab- or comma-separated: the same map Stage 1 writes")
    tx.add_argument("--ranks", default=DEFAULT_RANKS,
                    help="ranks that form the tiers, comma-separated, or 'all' (default %s)" % DEFAULT_RANKS)
    th = ap.add_argument_group("what counts as a warning (heuristics, adjustable)")
    th.add_argument("--thin", type=int, default=3, help="a tier with fewer species than this is THIN (default 3)")
    th.add_argument("--max-share", type=float, default=0.25, help="a species with a larger share of the sequences is DOMINANT (default 0.25)")
    th.add_argument("--min-seqs", type=int, default=500, help="a species with fewer sequences looks fragmentary: FEW_SEQUENCES (default 500)")
    ap.add_argument("-o", "--output", required=True, help="output folder")
    ap.add_argument("--tree", metavar="FILE", help="also write a Newick tree of the database species and your clade (needs the taxonomy)")
    ap.add_argument("--force", action="store_true", help="overwrite an earlier report in the output folder")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    warnings = []
    try:
        if bool(a.db) == bool(a.species_list):
            raise h.HepError("give the database with --db (FASTA files) or --species-list, not both and not neither")
        if a.tiers is None and a.taxonomy is None:
            raise h.HepError("give --taxonomy (folder from setup/get_taxonomy.sh) or --tiers (your own table)")
        if a.tiers is None and sum(x is not None for x in (a.clade, a.clade_taxid, a.clade_species, a.clade_inputs)) != 1:
            raise h.HepError("give exactly one of --clade, --clade-taxid, --clade-species, --clade-inputs")
        if a.tree and a.tiers is not None:
            raise h.HepError("--tree needs the taxonomy; it cannot be made from --tiers")
        h.check_output_folder(a.output, MODULE, MODULE, resume=False, force=a.force)
        if a.db:
            counts, untagged, odd = species_from_fasta(a.db)
            if not counts:
                raise h.HepError("no species tags (ID__Species_Name) found in the database; tag it with make-db, or use --species-list")
            if untagged:
                warnings.append("%d sequence IDs carry no usable species tag and are left out of this report%s"
                                % (untagged, " (%d of them contain '__' inside a random-looking ID, which is not read as a species)" % odd if odd else ""))
        else:
            counts, untagged, odd = species_from_list(a.species_list), 0, 0
        total = sum(v for v in counts.values() if v)     # over the species whose count is known

        rows, ladder_names, tiers_info, tree_tids, clade_desc = [], [], [], {}, None
        if a.tiers:
            tiermap = read_pairs(a.tiers, "tiers table")
            order = []
            for sp in counts:
                if sp not in tiermap:
                    warnings.append("%s is not in the tiers table" % sp)
            for t in list(tiermap.values()) + (a.expected_tiers.split(",") if a.expected_tiers else []):
                t = t.strip()
                if t and t not in order:
                    order.append(t)
            species_rows = {sp: (None, None, "given", tiermap.get(sp, UNMATCHED), None) for sp in counts}
            ladder = order
            tax = None
        else:
            tax = Taxonomy(a.taxonomy)
            taxid_map = {k: int(v) for k, v in read_pairs(a.taxid_map, "taxid map").items()} if a.taxid_map else {}
            if a.clade_taxid:
                clade = a.clade_taxid
                if not tax.has(clade):
                    raise h.HepError("taxid %d is not in the taxonomy" % clade)
            elif a.clade:
                f = tax.find_names({norm_name(a.clade)}).get(norm_name(a.clade), [])
                if not f:
                    raise h.HepError("the clade '%s' was not found in the taxonomy" % a.clade)
                if len(f) > 1:
                    raise h.HepError("the clade name '%s' matches %d taxa (%s); give one with --clade-taxid"
                                     % (a.clade, len(f), ", ".join(str(x) for x in f[:5])))
                clade = f[0]
            else:
                if a.clade_species:
                    with open(a.clade_species, "r", encoding="utf-8") as fh:
                        names = [re.sub(r"\s+", "_", l.strip()) for l in fh if l.strip() and not l.startswith("#")]
                else:
                    smap = h.load_species_map(a.species_map) if a.species_map else None
                    names = [re.sub(r"\s+", "_", sp) for sp, _f in h.discover_inputs(a.clade_inputs, species_map=smap)]
                if not names:
                    raise h.HepError("the list of analysis species is empty")
                res = resolve(tax, names, taxid_map)
                bad = [n for n in names if res[n][0] is None]
                if bad:
                    raise h.HepError("cannot place these analysis species on the taxonomy: %s. Use full names (--species-map for "
                                     "--clade-inputs) or --taxid-map" % ", ".join("%s (%s)" % (b, res[b][3]) for b in bad[:6]))
                clade = lca(tax, [res[n][0] for n in names])
            clade_desc = "%s (taxid %d, %s)" % (tax.sci_names([clade])[clade], clade, tax.rank(clade))
            ranks = "all" if a.ranks == "all" else {r.strip() for r in a.ranks.split(",") if r.strip()}
            ladder_ids = build_ladder(tax, clade, ranks)
            res = resolve(tax, list(counts), taxid_map)
            species_rows = {}
            for sp, (tid, matched, level, note) in res.items():
                if tid is None:
                    species_rows[sp] = (None, matched, level, UNMATCHED, note)
                    warnings.append("%s: %s" % (sp, note))
                    continue
                kind, idx = assign_tier(tax, tid, clade, ladder_ids)
                species_rows[sp] = (tid, matched, level, kind if kind == INSIDE else idx, None)
                if level == "genus":
                    warnings.append("%s was matched only at genus level (%s); check it" % (sp, matched))
                tree_tids[sp] = tid
            names_for = tax.sci_names(ladder_ids)
            ladder = ["%s (%s)" % (names_for[t], tax.rank(t)) for t in ladder_ids] + ["outside all of these (rest of the tree)"]
    except h.HepError as e:
        h.die(str(e))

    # ---- summaries
    tier_label = {}
    for sp, (tid, matched, level, tier, note) in species_rows.items():
        if a.tiers:
            tier_label[sp] = tier
        elif tier in (INSIDE, UNMATCHED):
            tier_label[sp] = tier
        else:
            tier_label[sp] = ladder[tier]
    order = ([INSIDE] if any(v == INSIDE for v in tier_label.values()) else []) + list(ladder) + \
            ([UNMATCHED] if any(v == UNMATCHED for v in tier_label.values()) else [])
    sp_out = [["species", "n_sequences", "pct_of_database", "taxid", "matched_as", "match_level", "tier", "flags"]]
    for sp in sorted(counts, key=lambda s: (order.index(tier_label[s]) if tier_label[s] in order else 99, s)):
        n = counts[sp]
        share = (n / total) if (n is not None and total) else None
        flags = []
        if tier_label[sp] == INSIDE:
            flags.append(INSIDE)
        if tier_label[sp] == UNMATCHED:
            flags.append(UNMATCHED)
        if share is not None and share > a.max_share:
            flags.append("DOMINANT")
        if n is not None and 0 < n < a.min_seqs and a.db:
            flags.append("FEW_SEQUENCES")
        tid, matched, level, _t, _n = species_rows[sp]
        sp_out.append([sp, n if n is not None else "NA", ("%.1f" % (100 * share)) if share is not None else "NA",
                       tid if tid else "NA", matched if matched else "NA", level, tier_label[sp], ";".join(flags) if flags else "NA"])
        for f in flags:
            if f in ("DOMINANT", "FEW_SEQUENCES", INSIDE):
                warnings.append("%s: %s%s" % (sp, f, {"DOMINANT": " (%.0f%% of all sequences)" % (100 * share) if share else "",
                                                       "FEW_SEQUENCES": " (%d sequences)" % n,
                                                       INSIDE: " (it belongs to your analysis clade; leave it out of the database)"}[f]))
    tier_rows = [["tier", "n_species", "n_sequences", "pct_of_database", "status"]]
    tier_status = {}
    for t in order:
        members = [sp for sp in counts if tier_label[sp] == t]
        n_sp = len(members)
        n_seq = sum(counts[s] for s in members) if all(counts[s] is not None for s in members) else None
        if t == INSIDE:
            status = "INSIDE_CLADE" if n_sp else "ok"
        elif t == UNMATCHED:
            status = UNMATCHED
        else:
            status = "EMPTY" if n_sp == 0 else ("THIN" if n_sp < a.thin else "ok")
        tier_status[t] = status
        tier_rows.append([t, n_sp, n_seq if n_seq is not None else "NA",
                          ("%.1f" % (100 * n_seq / total)) if (n_seq is not None and total) else "NA", status])
        if status in ("EMPTY", "THIN"):
            warnings.append("tier '%s' is %s (%d species)" % (t, status, n_sp))
    # ---- write
    os.makedirs(a.output, exist_ok=True)
    with open(os.path.join(a.output, h.RUN_MARKER), "w", encoding="utf-8") as fh:
        fh.write("module=%s\nlabel=%s\n" % (MODULE, MODULE))
    write_tsv(os.path.join(a.output, "species_tiers.tsv"), sp_out)
    write_tsv(os.path.join(a.output, "tier_summary.tsv"), tier_rows)
    if a.tree:
        try:
            tree_txt = newick(tax, tree_tids, clade, tax.sci_names([clade])[clade])
        except h.HepError as e:
            h.die(str(e))
        with open(a.tree, "w", encoding="utf-8") as fh:
            fh.write(tree_txt)
    info = {"module": MODULE, "pipeline": h.PIPELINE_NAME, "version": h.__version__, "clade": clade_desc,
            "taxonomy": (tax.source if tax else "not used (--tiers table)"),
            "databases": [os.path.abspath(p) for p in (a.db or [])], "species_list": a.species_list,
            "ranks": a.ranks, "thin": a.thin, "max_share": a.max_share, "min_seqs": a.min_seqs,
            "n_species": len(counts), "n_sequences": total,
            "tiers": {t: tier_status[t] for t in order}, "warnings": warnings}
    with open(os.path.join(a.output, "report.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
    if not a.quiet:
        h.log("%d species, %s sequences; clade: %s" % (len(counts), total if total else "NA", clade_desc or "(tiers table)"))
        for row in tier_rows[1:]:
            h.log("  %-45s %5s species  %s" % (row[0][:45], row[1], row[4]))
        for w in warnings[:15]:
            h.warn(w)
        if len(warnings) > 15:
            h.warn("... and %d more (all in report.json)" % (len(warnings) - 15))
        h.log("report written to %s. This is guidance, not a verdict: see docs/databases.md" % a.output)


if __name__ == "__main__":
    main()
