"""
step5_rates.py

Stage 5 step 5. Births and deaths of DNGs and PDNGs on every branch of the species tree, and their rates per unit of
branch length.

1. Filters. A gene is left out when a coding exon of its locus overlaps a coding exon of another annotated gene on
   the same strand (step 3: same_strand_cds_overlap; --overlap is required, --no-overlap only for tests and comparisons). Antisense loci and loci in an intron
   or UTR of a same-strand gene are counted; --exclude-overlap any leaves out every overlap class, as the paper did.
   A gene is also left out when it is a non-representative member of a family of
   unclustered duplicates (step 4; skipped without --families).
2. Origins. Loci (step 2) that share a hierarchical orthogroup (HOG) are copies of one origin followed by duplication,
   so they count once. A locus takes the HOG most frequent among its members; a locus with no member in the HOG table
   is an origin of its own. Within a HOG the DNG loci form one unit and the PDNG loci another. A PDNG locus is a
   duplicate of the DNG (its parent) and is dropped when every species in which it occurs also has a DNG of the same
   HOG; otherwise it stays in the PDNG unit. Each unit is born once, on the branch of its oldest locus (fewest
   ancestors; ties go to the locus listed first in loci.tsv, i.e. the older node, then the smaller geneID). A
   species-specific locus is born on that species' branch.
   Lineage split (--hog-span split, the default). An orthogroup groups proteins by similarity, not by position, so its
   loci can sit on branches that no single origin explains (e.g. one locus at node N5 and another in a species outside
   clade N5). The loci of a unit are taken oldest first; a locus joins the first origin whose clade contains its
   branch, otherwise it starts a new origin, with its own birth and deaths. Origins of one unit are labeled HOG,
   HOG#2, HOG#3... --hog-span merge counts one origin per unit, as in the paper.
3. Deaths (Dollo parsimony). For a unit born at internal node N, a species of clade(N) carries it when tBLASTn finds
   the ORF of any of the unit's genes in that genome (step 1); an internal node carries it when any descendant does;
   N carries it. A death is recorded on every branch below N whose parent carries the unit and which does not. This
   is the fewest deaths for one origin at N. Units born on a species branch have no deaths.
   With --unknown (T2T 2026-11-13), genomes listed there for a gene (alignment row missing or too sparse to call) take
   whichever state gives the fewest deaths (dollo_unknown): missing data never adds a death. The log gives the deaths
   this removed.
4. Rates = count / branch length, for every branch from the ingroup crown down, and for every species.

Gene IDs of loci.tsv are compared with the HOG table after removing isoform and '.p' suffixes, within a species
(--hog-id-match normalized, the default). The run stops when the IDs of a species differ in form from those of the
HOG table (no gene of the species matches, or some match only after case, prefixes, '_' vs '.' or version suffixes are
ignored). The Brassicaceae paper compared raw IDs across all species (exact), which missed the HOG of genes whose IDs
carry '.p' in the HOG table only (Boechera, Brassica rapa, Capsella, Thlaspi); use exact only to reproduce it.

This is ancestral_nodes_estimates_and_trees/generate_rates_noDNGparalog.py of the Brassicaceae analysis (the rate table
of the paper) with the tree, branch lengths and species read from inputs.

Usage:
  python step5_rates.py --tree-config tree_config.json --loci loci.tsv --genes genes.tsv --hogs N0.tsv \
      --branch-lengths lengths.tsv (--overlap overlap.tsv | --no-overlap) [--families unclustered_families.tsv] \
      [--hog-id-match normalized|exact] [--hog-span split|merge] [--exclude-overlap same_strand_cds|any] --outdir DIR
  lengths.tsv: columns branch, length (branch = internal node or species; length in the time unit of the rates,
  e.g. million years), or columns node, age (node ages; species are age 0).
Output (in DIR): rates_by_branch.tsv (Branch, Length, DNG/PDNG/TDNG births and deaths and their rates; rates NA where
        the length is missing or zero), origins.tsv (one row per counted unit: HOG, type, birth branch, death branches,
        number of loci), gene_status.tsv (every gene of loci.tsv and why it is counted or not)
"""
import argparse
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s5_common import (Tree, check_hog_ids, die, log, need, read_branch_lengths, read_tsv, strip_isoform,  # noqa: E402
                        write_tsv)

OVERLAP_REASON = {"same_strand_cds_overlap": "CDS overlaps a same-strand gene",
                  "same_strand_overlap": "overlaps a same-strand gene (intron or UTR)",
                  "antisense_overlap": "overlaps a gene on the opposite strand"}
HOG_META = {"HOG", "OG", "Gene Tree Parent Clade"}
LEAF_DEPTH = 10 ** 6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree-config", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--hogs", required=True)
    ap.add_argument("--branch-lengths", required=True)
    ov = ap.add_mutually_exclusive_group(required=True)
    ov.add_argument("--overlap", help="step 3 output")
    ov.add_argument("--no-overlap", action="store_true",
                    help="count without the gene-overlap filter (tests and comparisons only)")
    ap.add_argument("--unknown", default=None,
                    help="TSV species, geneID, unknown_in: genomes without a usable call, taken as carriers for deaths")
    ap.add_argument("--families")
    ap.add_argument("--hog-id-match", choices=("normalized", "exact"), default="normalized")
    ap.add_argument("--hog-span", choices=("split", "merge"), default="split")
    ap.add_argument("--exclude-overlap", choices=("same_strand_cds", "any"), default="same_strand_cds",
                    help="same_strand_cds (default): only CDS-CDS overlaps on the same strand are left out; "
                         "any: every overlap class (the paper)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    need(*[p for p in (a.tree_config, a.loci, a.genes, a.hogs, a.branch_lengths, a.overlap, a.families) if p])
    os.makedirs(a.outdir, exist_ok=True)
    for f in ("rates_by_branch.tsv", "origins.tsv", "gene_status.tsv"):
        if os.path.exists(os.path.join(a.outdir, f)):
            die(f"output exists: {os.path.join(a.outdir, f)}")
    tree = Tree(a.tree_config)
    blen = read_branch_lengths(a.branch_lengths, tree)

    presence = {}
    for r in read_tsv(a.genes):
        presence[(tree.canon(r["species"]), r["geneID"])] = set(filter(None, r["present_in"].split(";")))
    unknown = {}
    if a.unknown:
        need(a.unknown)
        for r in read_tsv(a.unknown):
            unknown[(tree.canon(r["species"]), r["geneID"])] = {tree.canon(x) for x in r["unknown_in"].split(";") if x}
        log(f"--unknown: {len(unknown):,} genes with genomes taken as carriers for deaths only")
    removed_by_unknown = Counter()

    excluded = {}
    overlap_classes = {"same_strand_cds_overlap"} | (
        {"antisense_overlap", "same_strand_overlap"} if a.exclude_overlap == "any" else set())
    if a.overlap:
        for r in read_tsv(a.overlap):
            if r["orientation"] in overlap_classes:
                excluded.setdefault((tree.canon(r["species"]), r.get("geneID", r.get("gene_id"))), r["orientation"])
    redundant = defaultdict(set)
    if a.families:
        for r in read_tsv(a.families):
            if str(r["is_representative"]).strip() == "False":
                redundant[tree.canon(r["species"])].add(r["gene"])

    hog_rows = read_tsv(a.hogs)
    cols = {}
    for c in hog_rows[0]:
        if c in HOG_META:
            continue
        s = tree.canon(c)
        if s not in tree.species:
            die(f"HOG table column {c!r} is not a species of tree_config (add it to the aliases)")
        cols[c] = s
    g2h, hog_ids, raw_seen = {}, defaultdict(set), {}
    for r in hog_rows:
        for c, sp in cols.items():
            for g in (r[c] or "").split(","):
                g = g.strip()
                if g:
                    if raw_seen.setdefault((sp, g), r["HOG"]) != r["HOG"]:
                        die(f"{a.hogs}: gene {g} ({c}) is in two HOGs, {raw_seen[(sp, g)]} and {r['HOG']}")
                    g2h[g if a.hog_id_match == "exact" else (sp, strip_isoform(g))] = r["HOG"]
                    hog_ids[sp].add(g)
    hkey = (lambda sp, g: g) if a.hog_id_match == "exact" else (lambda sp, g: (sp, strip_isoform(g)))

    loci = read_tsv(a.loci)
    sp_cols = [c for c in loci[0] if c not in ("ClusterID", "Node", "Gene_type")]
    for c in sp_cols:
        if tree.canon(c) not in tree.species:
            die(f"loci.tsv column {c} is not a species of tree_config")
    if a.hog_id_match == "normalized":
        genes_by_sp = defaultdict(set)
        for r in loci:
            for c in sp_cols:
                genes_by_sp[tree.canon(c)].update(g.strip() for g in (r[c] or "").replace(",", ";").split(";")
                                                  if g.strip())
        check_hog_ids(tree, genes_by_sp, hog_ids, lambda sp, g: (sp, strip_isoform(g)))

    status = {}
    units = defaultdict(lambda: {"dng": [], "pdng": [], "dng_species": set()})
    for r in loci:
        gt, node, cid = r["Gene_type"], r["Node"], r["ClusterID"]
        members, species = [], set()
        for c in sp_cols:
            sp = tree.canon(c)
            for g in (r[c] or "").replace(",", ";").split(";"):
                g = g.strip()
                if not g:
                    continue
                if (sp, g) in excluded:
                    status[(sp, g)] = [cid, node, gt, "", "excluded: " + OVERLAP_REASON[excluded[(sp, g)]]]
                elif strip_isoform(g) in redundant.get(sp, set()):
                    status[(sp, g)] = [cid, node, gt, "", "excluded: duplicate of an unclustered locus"]
                else:
                    members.append((sp, g))
                    species.add(sp)
        if gt not in ("DNG", "PDNG"):
            for m in members:
                status[m] = [cid, node, gt, "", f"not counted: gene type {gt or 'empty'}"]
            continue
        if not members:
            continue
        br = tree.branch_of(node)
        if br is None:
            for m in members:
                status[m] = [cid, node, gt, "", "excluded: ORF node outside the window of tree_config"]
            continue
        hogs = [g2h[hkey(*m)] for m in members if hkey(*m) in g2h]
        hk = Counter(hogs).most_common(1)[0][0] if hogs else f"SINGLETON_{cid}"
        depth = tree.depth[br] if br in tree.clades else LEAF_DEPTH
        u = units[hk]
        if gt == "DNG":
            u["dng"].append((depth, br, members, cid, species))
            u["dng_species"].update(species)
        else:
            u["pdng"].append((depth, br, species, members, cid))
        for m in members:
            status[m] = [cid, node, gt, hk, ""]

    def dollo(origin, pres, counter):
        if origin not in tree.clades:
            return []
        clade = tree.subtree(origin)
        leaves = [x for x in clade if x not in tree.clades]
        pa = {s: int(s in pres) for s in leaves}
        for x in tree.postorder(origin):
            if x in tree.clades:
                pa[x] = max(pa[c] for c in tree.children[x])
        pa[origin] = 1
        lost = [x for x in clade if x != origin and pa[tree.parent[x]] == 1 and pa[x] == 0]
        for x in lost:
            counter[x] += 1
        return lost

    def dollo_unknown(origin, pres, unk, counter):
        """Dollo with genomes of unknown state: each unknown genome takes the state (carrier or not) that gives the
        fewest deaths for the origin, so missing data never adds a death. c1/c0 = fewest deaths below a node that
        carries / does not carry (c0 possible only when no carrier genome is below it). Ties put the death on the
        higher branch, as plain Dollo does. With no unknown genome this equals dollo() (checked on every call)."""
        if origin not in tree.clades:
            return []
        INF = float("inf")
        c1, c0 = {}, {}
        for x in tree.postorder(origin):
            if x in tree.clades:
                c1[x] = sum(min(c1[ch], c0[ch] + 1) for ch in tree.children[x])
                c0[x] = 0 if all(c0[ch] == 0 for ch in tree.children[x]) else INF
            else:
                c1[x] = 0 if (x in pres or x in unk) else INF
                c0[x] = INF if x in pres else 0
        lost, stack = [], [origin]
        while stack:
            x = stack.pop()
            for ch in tree.children.get(x, ()) if x in tree.clades else ():
                if c0[ch] + 1 <= c1[ch]:
                    lost.append(ch)
                else:
                    stack.append(ch)
        order = {n: i for i, n in enumerate(tree.subtree(origin))}
        lost.sort(key=lambda n: order[n])
        for x in lost:
            counter[x] += 1
        return lost

    def pres_of(member_lists):
        out = set()
        for ms in member_lists:
            for m in ms:
                out |= presence.get(m, set())
        return out

    below = {n: set(tree.subtree(n)) for n in tree.clades}

    def inside(br, origin):
        return br == origin or br in below.get(origin, ())

    def lineages(items):
        """items: (depth, branch, members) in loci.tsv order. Returns [(birth branch, [members lists])]."""
        order = sorted(items, key=lambda t: t[0])
        if a.hog_span == "merge":
            return [(order[0][1], [t[2] for t in order])]
        groups = []
        for d, br, ms in order:
            for g in groups:
                if inside(br, g[0]):
                    g[1].append(ms)
                    break
            else:
                groups.append((br, [ms]))
        return groups

    births = {"DNG": Counter(), "PDNG": Counter()}
    deaths = {"DNG": Counter(), "PDNG": Counter()}
    origins = []
    for hk, u in units.items():
        surv = [t for t in u["pdng"] if t[2] - u["dng_species"]]
        for t in u["pdng"]:
            if t not in surv:
                for m in t[3]:
                    status[m][4] = "excluded: PDNG paralog of a DNG in the same species (DNG kept as parent)"
        for gt, items in (("DNG", [(t[0], t[1], t[2]) for t in u["dng"]]),
                          ("PDNG", [(t[0], t[1], t[3]) for t in surv])):
            if not items:
                continue
            for k, (br, mls) in enumerate(lineages(items), 1):
                label = hk if k == 1 else f"{hk}#{k}"
                births[gt][br] += 1
                pres = pres_of(mls)
                if unknown:
                    unk = set().union(*[unknown.get(m, set()) for ms in mls for m in ms]) - pres
                    strict = dollo(br, pres, Counter())
                    if sorted(dollo_unknown(br, pres, set(), Counter())) != sorted(strict):
                        die(f"dollo_unknown without unknown genomes differs from dollo for {hk} at {br}")
                    lost = dollo_unknown(br, pres, unk, deaths[gt])
                    if len(lost) > len(strict):
                        die(f"{hk}: unknown genomes added a death")
                    removed_by_unknown[gt] += len(strict) - len(lost)
                else:
                    lost = dollo(br, pres, deaths[gt])
                origins.append([label, gt, br, ";".join(lost), len(mls)])
                for ms in mls:
                    for m in ms:
                        status[m][3] = label
                        status[m][4] = f"counted: {gt} unit of its HOG"

    def fmt(v):
        return f"{v:.3f}"

    rows = []
    for b in tree.report_branches():
        bl = blen.get(b)
        db, pb, dd, pd_ = births["DNG"][b], births["PDNG"][b], deaths["DNG"][b], deaths["PDNG"][b]
        tb, td = db + pb, dd + pd_
        rates = ["NA"] * 6 if not bl else [fmt(x / bl) for x in (db, pb, dd, pd_, tb, td)]
        rows.append([b, "NA" if bl is None else bl, db, pb, dd, pd_, tb, td] + rates)
    off = {b for c in (births, deaths) for t in c.values() for b in t} - set(tree.report_branches())
    if off:
        die(f"events on branches outside the table: {sorted(off)}")
    write_tsv(os.path.join(a.outdir, "rates_by_branch.tsv"),
              ["Branch", "Length", "DNG_births", "PDNG_births", "DNG_deaths", "PDNG_deaths", "TDNG_births",
               "TDNG_deaths", "DNG_birth_rate", "PDNG_birth_rate", "DNG_death_rate", "PDNG_death_rate",
               "TDNG_birth_rate", "TDNG_death_rate"], rows)
    write_tsv(os.path.join(a.outdir, "origins.tsv"), ["HOG", "Gene_type", "birth_branch", "death_branches", "loci"],
              origins)
    write_tsv(os.path.join(a.outdir, "gene_status.tsv"),
              ["species", "geneID", "ClusterID", "Cluster_node", "Gene_type", "HOG_unit", "status"],
              [[sp, g] + v for (sp, g), v in sorted(status.items())])

    n_genes = sum(len([g for g in (r[c] or "").replace(",", ";").split(";") if g.strip()]) for r in loci for c in sp_cols)
    st = Counter(v[4] for v in status.values())
    assert sum(st.values()) == len(status) == n_genes, (len(status), n_genes)
    assert all(sum(births[t].values()) == sum(1 for o in origins if o[1] == t) for t in births)
    log(f"{n_genes:,} genes in {len(loci):,} loci; " + "; ".join(f"{k or 'not counted'}: {v:,}" for k, v in st.items()))
    for t in ("DNG", "PDNG"):
        log(f"{t}: {sum(births[t].values()):,} births, {sum(deaths[t].values()):,} deaths")
        if unknown:
            log(f"{t}: {removed_by_unknown[t]:,} deaths not counted because the genome had no usable call (--unknown)")
    log(f"wrote {a.outdir}/rates_by_branch.tsv, origins.tsv, gene_status.tsv")


if __name__ == "__main__":
    main()
