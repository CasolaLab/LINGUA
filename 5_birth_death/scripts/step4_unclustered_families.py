"""
step4_unclustered_families.py

Stage 5 step 4. Duplicates among the loci that OrthoFinder left out of every orthogroup.

Step 5 collapses loci that share a hierarchical orthogroup (HOG) into one origin. Loci absent from the HOG table
("unclustered") have no orthogroup, so duplicates among them would each count as an origin. This step finds them:
  select    writes, per species, the proteins of its unclustered DNG and PDNG loci (geneIDs compared without isoform
            or '.p' suffixes, anchored at the end of the ID)
  (BLASTp)  all against all within each species, E <= 1e-5, one HSP per pair (run_stage5.sh does this)
  families  joins loci connected by any BLASTp hit, directly or through other loci (connected components), and keeps
            the longest protein of each family as its representative (ties: alphabetically first ID); the other
            members are duplicates and step 5 leaves them out
Species with fewer than two unclustered loci found in their protein file are not searched. The run stops when the
gene IDs of a species differ in form from those of the HOG table (see step5_rates.py), and when fewer than 95% of a
species' unclustered loci are found in its protein file (--min-protein-match).

This is unclustered_paralog_blast.py and build_unclustered_families.py of the Brassicaceae analysis without their
species lists.

Usage:
  python step4_unclustered_families.py select --tree-config tree_config.json --loci loci.tsv --hogs N0.tsv \
      --proteins DIR --outdir DIR
  python step4_unclustered_families.py families --outdir DIR
  --hogs: an OrthoFinder Phylogenetic_Hierarchical_Orthogroups table (N0.tsv, or the table of the node that holds all
  species); its species columns must be species names or aliases of tree_config.
  --proteins: the stage 4 candidate protein FASTAs, one per species, headers '>geneID' or '>geneID__Species'.
Output (in DIR): <species>_unclustered.faa, unclustered_counts.tsv (select); unclustered_families.tsv (species, gene,
        family_id, family_size, is_representative, protein_length) and unclustered_family_summary.tsv (families)
"""
import argparse
import glob
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s5_common import Tree, check_hog_ids, die, log, need, read_tsv, strip_isoform, write_tsv  # noqa: E402

HOG_META = {"HOG", "OG", "Gene Tree Parent Clade"}


def hog_columns(tree, header):
    cols = {}
    for c in header:
        if c in HOG_META:
            continue
        s = tree.canon(c)
        if s not in tree.species:
            die(f"HOG table column {c!r} is not a species of tree_config (add it to the aliases)")
        cols[c] = s
    return cols


def species_of_file(tree, fname):
    parts = re.split(r"\.(faa|fa|fasta)$", os.path.basename(fname))[0].split("_")
    for k in range(len(parts), 0, -1):
        s = tree.canon("_".join(parts[:k]))
        if s in tree.species:
            return s
    return None


def read_fasta(path):
    """Sequences keyed by strip_isoform(ID); also returns the keys that more than one header maps to."""
    seqs, cur, buf, n = {}, None, [], defaultdict(int)
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if cur is not None:
                seqs[cur] = "".join(buf)
            cur, buf = strip_isoform(line[1:].split("__")[0].strip()), []
            n[cur] += 1
        elif line:
            buf.append(line)
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs, {k for k, v in n.items() if v > 1}


def select(a):
    need(a.tree_config, a.loci, a.hogs, a.proteins)
    tree = Tree(a.tree_config)
    os.makedirs(a.outdir, exist_ok=True)
    if glob.glob(os.path.join(a.outdir, "*_unclustered.faa")):
        die(f"{a.outdir} already holds *_unclustered.faa")
    hog_rows = read_tsv(a.hogs)
    cols = hog_columns(tree, hog_rows[0].keys())
    in_hog, raw_hog = defaultdict(set), defaultdict(set)
    for r in hog_rows:
        for c, sp in cols.items():
            for g in (r[c] or "").split(","):
                if g.strip():
                    in_hog[sp].add(strip_isoform(g.strip()))
                    raw_hog[sp].add(g.strip())
    unc, all_genes = defaultdict(set), defaultdict(set)
    for r in read_tsv(a.loci):
        if r["Gene_type"] not in ("DNG", "PDNG"):
            continue
        for sp in tree.species:
            for g in (r.get(sp) or "").replace(",", ";").split(";"):
                if not g.strip():
                    continue
                all_genes[sp].add(g.strip())
                if strip_isoform(g.strip()) not in in_hog[sp]:
                    unc[sp].add(strip_isoform(g.strip()))
    check_hog_ids(tree, all_genes, raw_hog, lambda sp, g: (sp, strip_isoform(g)))
    files = {}
    for f in sorted(glob.glob(os.path.join(a.proteins, "*"))):
        if os.path.basename(f).startswith("._") or not re.search(r"\.(faa|fa|fasta)$", f):
            continue
        s = species_of_file(tree, f)
        if s is None:
            die(f"cannot tell the species of {f} from its name")
        if s in files:
            die(f"two protein files for {s}: {files[s]}, {f}")
        files[s] = f
    counts = []
    for sp in tree.species:
        genes = unc.get(sp, set())
        if not genes:
            counts.append([sp, 0, 0, False])
            continue
        if sp not in files:
            die(f"no protein file for {sp}, which has {len(genes)} unclustered loci")
        seqs, multi = read_fasta(files[sp])
        clash = sorted(genes & multi)
        if clash:
            die(f"{files[sp]}: several sequences for the unclustered loci {clash[:5]} (isoforms, or IDs that differ "
                "only in a trailing number); keep one protein per locus")
        found = {g: seqs[g] for g in genes if g in seqs}
        if len(found) < len(genes):
            log(f"  {sp}: {len(genes) - len(found)} of {len(genes)} unclustered loci not in {files[sp]}")
        if len(found) < a.min_protein_match * len(genes):
            die(f"{sp}: only {len(found)} of {len(genes)} unclustered loci found in {files[sp]} "
                f"(e.g. missing {sorted(set(genes) - set(found))[:3]}; FASTA IDs {sorted(seqs)[:2]})")
        counts.append([sp, len(genes), len(found), len(found) >= 2])
        if len(found) >= 2:
            with open(os.path.join(a.outdir, f"{sp}_unclustered.faa"), "w") as fh:
                for g in sorted(found):
                    fh.write(f">{g}\n{found[g]}\n")
    write_tsv(os.path.join(a.outdir, "unclustered_counts.tsv"), ["species", "unclustered", "in_fasta", "blast_run"],
              counts)
    for c in counts:
        log("\t".join(map(str, c)))


def families(a):
    need(a.outdir)
    out_f = os.path.join(a.outdir, "unclustered_families.tsv")
    if os.path.exists(out_f):
        die(f"output exists: {out_f}")
    gene_rows, summary = [], []
    for faa in sorted(glob.glob(os.path.join(a.outdir, "*_unclustered.faa"))):
        sp = os.path.basename(faa)[: -len("_unclustered.faa")]
        bl = os.path.join(a.outdir, f"{sp}_unclustered.blastp.tsv")
        need(bl)
        lengths, cur = {}, None
        for line in open(faa):
            line = line.rstrip()
            if line.startswith(">"):
                cur = line[1:].strip()
                lengths[cur] = 0
            else:
                lengths[cur] += len(line)
        parent = {g: g for g in lengths}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        n_pairs = 0
        for line in open(bl):
            f = line.split("\t")
            q, s = f[0], f[1]
            if q == s or q not in lengths or s not in lengths:
                continue
            n_pairs += 1
            parent[find(q)] = find(s)
        comps = defaultdict(list)
        for g in lengths:
            comps[find(g)].append(g)
        comps = sorted((sorted(c) for c in comps.values()), key=lambda c: (-len(c), c[0]))
        for i, members in enumerate(comps, 1):
            rep = min(members, key=lambda g: (-lengths[g], g))
            for g in members:
                gene_rows.append([sp, g, f"{sp}_FAM{i:04d}", len(members), g == rep, lengths[g]])
        multi = [c for c in comps if len(c) > 1]
        red = sum(len(c) for c in multi) - len(multi)
        summary.append([sp, len(lengths), len(comps), len(comps) - len(multi), len(multi), sum(len(c) for c in multi),
                        red, round(100 * red / len(lengths), 1)])
        log(f"  {sp}: {len(lengths)} loci, {n_pairs} hit pairs -> {len(comps)} families, {red} duplicates")
    write_tsv(out_f, ["species", "gene", "family_id", "family_size", "is_representative", "protein_length"], gene_rows)
    write_tsv(os.path.join(a.outdir, "unclustered_family_summary.tsv"),
              ["species", "total_unclustered", "total_families", "singleton_families", "multi_gene_families",
               "genes_in_multi_fam", "redundant_genes", "pct_redundant"], summary)
    log(f"wrote {out_f}; duplicates left out in step 5: {sum(s[6] for s in summary):,}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    s = sub.add_parser("select")
    s.add_argument("--tree-config", required=True)
    s.add_argument("--loci", required=True, help="step 2 loci.tsv")
    s.add_argument("--hogs", required=True)
    s.add_argument("--proteins", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--min-protein-match", type=float, default=0.95)
    f = sub.add_parser("families")
    f.add_argument("--outdir", required=True)
    a = ap.parse_args()
    select(a) if a.mode == "select" else families(a)


if __name__ == "__main__":
    main()
