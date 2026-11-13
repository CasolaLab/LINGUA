"""
step3_overlap.py

Stage 5 step 3 (optional; needs the genome annotations). Classifies each classified locus by its position relative
to the protein-coding genes of the species' GFF3 that are not themselves classified loci.
Overlap classes (by priority):
  same_strand_cds_overlap  a coding exon of the locus overlaps a coding exon of a gene on the same strand. The locus
                           may be a fragment or an alternative frame of that gene; step 5 leaves it out.
  same_strand_overlap      the spans overlap on the same strand but the coding exons do not (the locus lies in an
                           intron or a UTR of the other gene); counted.
  antisense_overlap        the spans overlap on opposite strands, any part of the other gene; counted.
Step 5 with --exclude-overlap any leaves out all three (the paper). Loci that overlap no other gene are classified by
the nearest neighbor within 1 kb (up_same, down_same, divergent, convergent) or as intergenic; counted. Two classes
mean the locus could not be checked, and it is counted: unknown (no GFF gene at the locus, or its chromosome is not in
the GFF) and own_gene_not_found (GFF genes overlap the locus but none carries its geneID).

Own gene. The locus's own gene is the overlapping GFF gene that carries the locus geneID: as a transcript ID or CDS
protein_id (NP_..., XP_..., ENSMUST...); else as a transcript ID followed by a separator (Thhalv10029120m ->
Thhalv10029120m.v1.0), also after removing isoform suffixes (AT1G45163.1 -> AT1G45163); else as a whole word of the
gene ID (Thecc1EG008161t1 -> Thecc1EG008161). The own gene is never guessed from coordinates, because a long host
gene can cover the locus as well as its own gene. The run stops when more than --max-unmatched (default 5%) of a
species' loci overlapping GFF genes have no own gene, which means the IDs of the GFF and of the loci differ in form.
The coding exons of the locus are the CDS of the matched transcript (all transcripts of the gene when it was found by
gene ID, or when the locus isoform is not in the GFF; the log counts these), within the locus span; the span itself
when there are none. The strand of the locus is that of its own gene.

Coordinates from step 1 are 0-based starts, GFF3 is 1-based. Every gene overlapping a locus is found, however far
upstream it starts. A locus chromosome missing from the GFF is matched ignoring case (chr14 / Chr14). When a GFF names
chromosomes differently from the genome used in stage 3 (GenBank vs RefSeq accessions), give that species' .fai: a GFF
chromosome absent from the .fai is renamed to the .fai sequence of the same length (from the GFF ##sequence-region
lines) only when that length is unique in both files; genes on chromosomes whose length is shared are dropped
(logged), so their loci are unknown rather than compared with genes of another scaffold. A GFF chromosome without a
##sequence-region line keeps its name; a GFF with none of these lines stops the run when a .fai is given.
The run also stops when more than --max-unknown (default 5%) of a species' loci are unknown, which means the GFF
belongs to another assembly or names chromosomes differently.

Only protein-coding genes count: a GFF gene whose biotype is stated (gene_biotype, biotype or gene_type) and is not
protein_coding is ignored; genes without a stated biotype are kept.

This is orientation_dngs/orientation_all_species.py of the Brassicaceae analysis without its species lists, figure
and statistics, with these changes (2026-11-12): protein-coding genes only; a complete overlap search (the original
stopped at the first upstream gene that ends before the locus, so it missed loci inside long genes); the own gene by
ID instead of by overlap; CDS-level same-strand overlaps; the nearest upstream neighbor is the upstream gene that ends
closest to the locus; the .fai renaming refuses lengths shared by several sequences.

Usage:
  python step3_overlap.py --tree-config tree_config.json --coords coords.tsv --gff-table gffs.tsv --out overlap.tsv \
      [--max-unmatched 0.05] [--max-unknown 0.05]
  gffs.tsv: columns species, gff and optionally fai (paths relative to the table's folder or absolute)
Output: overlap.tsv (species, geneID, orientation, other_gene, chrom, start, end), one row per own-genome coordinate
        row of a locus; a locus is excluded in step 5 when any of its rows is an excluded overlap class
"""
import argparse
import bisect
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s5_common import Tree, die, log, need, read_tsv, strip_p, write_tsv  # noqa: E402

GAP_THRESH = 1000


def chr_map_from_fai(gff_path, fai_path):
    """GFF chromosome -> .fai name. Names present in the .fai are kept; others are renamed by length only when that
    length occurs once in the GFF and once in the .fai; the rest map to None (their genes are dropped)."""
    gff_sizes = {}
    with open(gff_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("##sequence-region"):
                p = line.split()
                if len(p) >= 4:
                    gff_sizes[p[1]] = int(p[3])
    if not gff_sizes:
        die(f"{gff_path}: a .fai was given but the GFF has no ##sequence-region lines, so its chromosomes cannot be "
            "renamed by length")
    fai_sizes = {l.split()[0]: int(l.split()[1]) for l in open(fai_path) if l.strip()}
    n_fai, n_gff = Counter(fai_sizes.values()), Counter(gff_sizes.values())
    by_size = {s: c for c, s in fai_sizes.items()}
    out = {}
    for c, s in gff_sizes.items():
        if c in fai_sizes:
            out[c] = c
        elif n_fai[s] == 1 and n_gff[s] == 1:
            out[c] = by_size[s]
        else:
            out[c] = None
    return out


BIOTYPE = re.compile(r"(?:gene_biotype|biotype|gene_type)=([^;]+)")
ATTR = {k: re.compile(r"(?:^|;)\s*%s=([^;]+)" % k) for k in ("ID", "Parent", "protein_id")}


class Genome(dict):
    """chromosome -> (genes sorted by start, starts, running maximum of ends); .lower maps lower-case names."""

    def __init__(self, *a):
        super().__init__(*a)
        self.lower = {c.lower(): c for c in self}


def load_gff(path, chr_map):
    """Protein-coding genes of a GFF3 as a Genome. Each gene is [start, end, strand, ID, CDS of all its transcripts,
    aliases, CDS by transcript]; aliases maps the lower-case IDs of its transcripts (with and without an rna-/cds-/
    transcript- prefix) and the protein_id of its CDS to the transcript they belong to. CDS are linked to genes
    through the Parent chain (CDS -> mRNA/transcript -> gene). Returns (genome, non-coding genes skipped, genes
    dropped because their chromosome could not be renamed)."""
    genes, skipped, dropped, parent, cds = {}, 0, 0, {}, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            ftype = p[2]
            if ftype == "exon":
                continue
            pa = ATTR["Parent"].search(p[8])
            if ftype == "CDS":
                if pa:
                    pid = ATTR["protein_id"].search(p[8])
                    cds.append((pa.group(1).split(",")[0], int(p[3]), int(p[4]), pid.group(1) if pid else None))
                continue
            m = ATTR["ID"].search(p[8])
            if m and pa:
                parent[m.group(1)] = pa.group(1).split(",")[0]
            if ftype != "gene":
                continue
            b = BIOTYPE.search(p[8])
            if b and b.group(1) != "protein_coding":
                skipped += 1
                continue
            chrom = chr_map.get(p[0], p[0]) if chr_map else p[0]
            if chrom is None:
                dropped += 1
                continue
            gid = m.group(1) if m else f"unnamed_gene_{p[0]}_{p[3]}_{p[4]}_{p[6]}"
            genes.setdefault(chrom, []).append([int(p[3]), int(p[4]), p[6], gid, [], {}, {}])
    gene_ids = {g[3]: g for gs in genes.values() for g in gs}

    def gene_of(x):
        seen = 0
        while x not in gene_ids and x in parent and seen < 5:
            x, seen = parent[x], seen + 1
        return gene_ids.get(x)

    for tid in parent:
        g = gene_of(tid)
        if g is not None:
            g[5][tid.lower()] = tid
            g[5][re.sub(r"^(rna|cds|transcript)[-:]", "", tid.lower())] = tid
    for par, s, e, pid in cds:
        g = gene_of(par)
        if g is not None:
            g[4].append((s, e))
            g[6].setdefault(par, []).append((s, e))
            if pid:
                g[5][pid.lower()] = par
    out = Genome()
    for c, gs in genes.items():
        gs.sort(key=lambda x: x[0])
        mx, run = [], 0
        for g in gs:
            run = max(run, g[1])
            mx.append(run)
        out[c] = (gs, [g[0] for g in gs], mx)
    out.lower = {c.lower(): c for c in out}
    return out, skipped, dropped


def resolve_chrom(chrom, genome):
    """The GFF chromosome of a locus chromosome: the same name, else the same name ignoring case (chr14 / Chr14)."""
    return chrom if chrom in genome else genome.lower.get(chrom.lower())


def overlapping(entry, q1, q_end):
    """Indices of every gene overlapping [q1, q_end] (1-based, inclusive), nearest start first."""
    gs, starts, mx = entry
    i = bisect.bisect_right(starts, q_end) - 1
    out = []
    while i >= 0 and mx[i] >= q1:
        if gs[i][1] >= q1:
            out.append(i)
        i -= 1
    return out


def id_forms(gid):
    """Forms of a locus geneID, most specific first (AT1G45163.1 -> at1g45163.1, at1g45163; Aa_G876160.h2.t1 ->
    aa_g876160.h2.t1, aa_g876160.h2; Thecc1EG008161t1 -> thecc1eg008161t1, thecc1eg008161)."""
    g = strip_p(gid).lower()
    out = [gid.lower(), g, re.sub(r"\.\d+$", "", g)]
    out.append(re.sub(r"[._]?t\d+$", "", out[-1]))
    return [x for x in dict.fromkeys(out) if len(x) >= 4]


def _prefix(alias, f):
    return alias.startswith(f) and (len(alias) == len(f) or not alias[len(f)].isalnum())


def _word(f, text):
    return re.search(r"(?:^|[^a-z0-9])" + re.escape(f) + r"(?:$|[^a-z0-9])", text) is not None


def own_gene(chrom, q_start0, q_end, genome, gid):
    """(status, gene, CDS of the locus, CDS source) with status no_chrom, no_gene, unmatched or named, and CDS source
    'one transcript' or 'several transcripts' (the locus isoform is not in the GFF, or the gene was matched by its
    gene ID, so the CDS of every matching transcript are used). See the module notes."""
    chrom = resolve_chrom(chrom, genome)
    if chrom is None:
        return "no_chrom", None, None, None
    entry = genome[chrom]
    idx = sorted(overlapping(entry, q_start0 + 1, q_end))
    if not idx:
        return "no_gene", None, None, None
    gs = [entry[0][i] for i in idx]
    g0 = gid.lower()
    for g in gs:
        if g0 in g[5]:
            t = g[5][g0]
            return ("named", g, g[6][t], "one transcript") if g[6].get(t) else ("named", g, g[4], "several transcripts")
    forms = id_forms(gid)
    for f in forms:
        for g in gs:
            tids = sorted({t for a, t in g[5].items() if _prefix(a, f)})
            with_cds = [t for t in tids if g[6].get(t)]
            if tids and with_cds:
                return "named", g, [x for t in with_cds for x in g[6][t]], \
                    "one transcript" if len(with_cds) == 1 else "several transcripts"
            if tids:
                return "named", g, g[4], "several transcripts"
    for f in forms:
        for g in gs:
            if _word(f, g[3].lower()):
                return "named", g, g[4], "one transcript" if len(g[6]) <= 1 else "several transcripts"
    return "unmatched", None, None, None


def classify_orientation(t_s, t_e, t_strand, g_s, g_e, g_strand):
    if not (t_e < g_s or t_s > g_e):
        return "antisense_overlap" if t_strand != g_strand else "same_strand_overlap"
    t_left = t_e < g_s
    if t_strand == g_strand:
        if g_strand == "+":
            return "down_same" if not t_left else "up_same"
        return "down_same" if t_left else "up_same"
    if g_strand == "+":
        return "divergent" if t_left else "convergent"
    return "convergent" if t_left else "divergent"


def classify(rows, genome, locus_gff_ids):
    """One row per locus row: (geneID, class, other gene, chrom, start, end). Overlap classes, by priority:
    same_strand_cds_overlap (a CDS of the locus overlaps a CDS of a same-strand gene), same_strand_overlap (the spans
    overlap on the same strand, CDS not: intronic or UTR), antisense_overlap (the spans overlap on opposite strands)."""
    out = []
    for gid, chrom0, t_s0, t_e in rows:
        t_s1 = t_s0 + 1
        status, own, own_cds, _ = own_gene(chrom0, t_s0, t_e, genome, gid)
        if status in ("no_chrom", "no_gene"):
            out.append((gid, "unknown", "", chrom0, t_s1, t_e))
            continue
        if status == "unmatched":
            out.append((gid, "own_gene_not_found", "", chrom0, t_s1, t_e))
            continue
        t_strand = own[2]
        entry = genome[resolve_chrom(chrom0, genome)]
        genes, starts = entry[0], entry[1]
        ov = [genes[i] for i in overlapping(entry, t_s1, t_e) if genes[i][3] not in locus_gff_ids]
        q = [(max(cs, t_s1), min(ce, t_e)) for cs, ce in own_cds if cs <= t_e and ce >= t_s1] or [(t_s1, t_e)]
        same = [g for g in ov if g[2] == t_strand]
        same_cds = [g for g in same if any(cs <= qe and ce >= qs for cs, ce in g[4] for qs, qe in q)]
        if same_cds:
            out.append((gid, "same_strand_cds_overlap", same_cds[0][3], chrom0, t_s1, t_e))
            continue
        if same:
            out.append((gid, "same_strand_overlap", same[0][3], chrom0, t_s1, t_e))
            continue
        if ov:
            out.append((gid, "antisense_overlap", ov[0][3], chrom0, t_s1, t_e))
            continue
        idx_r = bisect.bisect_right(starts, t_e)
        left = [g for g in genes[:idx_r] if g[1] < t_s1 and g[3] not in locus_gff_ids]
        left_cand = max(left, key=lambda g: g[1]) if left else None
        right_cand = next((genes[i] for i in range(idx_r, len(genes)) if genes[i][3] not in locus_gff_ids), None)
        nb = []
        if left_cand:
            nb.append((t_s1 - left_cand[1], left_cand))
        if right_cand:
            nb.append((right_cand[0] - t_e, right_cand))
        if not nb:
            o, other = "intergenic", ""
        else:
            gap, g = min(nb, key=lambda x: x[0])
            o, other = (classify_orientation(t_s1, t_e, t_strand, *g[:3]), g[3]) if gap <= GAP_THRESH \
                else ("intergenic", "")
        out.append((gid, o, other, chrom0, t_s1, t_e))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree-config", required=True)
    ap.add_argument("--coords", required=True, help="step 1 output")
    ap.add_argument("--gff-table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-gff", default="",
                    help="comma-separated species that have no annotation; their loci are not checked (any other "
                         "species of the coordinates without a GFF stops the run)")
    ap.add_argument("--max-unmatched", type=float, default=0.05,
                    help="stop when a larger share of a species' loci overlapping GFF genes has no own gene")
    ap.add_argument("--max-unknown", type=float, default=0.05,
                    help="stop when a larger share of a species' loci lies on no GFF gene or on a chromosome missing "
                         "from the GFF")
    a = ap.parse_args()
    need(a.tree_config, a.coords, a.gff_table)
    if os.path.exists(a.out):
        die(f"output exists: {a.out}")
    tree = Tree(a.tree_config)
    base = os.path.dirname(os.path.abspath(a.gff_table))
    table = {}
    for r in read_tsv(a.gff_table):
        sp = tree.canon(r["species"])
        if sp not in tree.species:
            die(f"{a.gff_table}: {r['species']} is not a species of tree_config")
        res = lambda p: p if not p or os.path.isabs(p) else os.path.join(base, p)
        table[sp] = (res(r["gff"]), res(r.get("fai", "") or ""))
        need(*[p for p in table[sp] if p])
    coords = read_tsv(a.coords)
    focal = [s for s in tree.species if any(tree.canon(r["focal_species"]) == s for r in coords)]
    opt_out = {tree.canon(x.strip()) for x in a.no_gff.split(",") if x.strip()}
    bad = sorted(x for x in opt_out if x not in tree.species)
    if bad:
        die(f"--no-gff: {bad} not species of tree_config")
    both = sorted(opt_out & set(table))
    if both:
        die(f"{both} have a GFF in {a.gff_table} and are also listed in --no-gff; choose one")
    no_gff = [s for s in focal if s not in table and s not in opt_out]
    if no_gff:
        die(f"no GFF for {no_gff} in {a.gff_table}; add their annotation, or name them in --no-gff (NO_GFF_SPECIES in "
            f"the driver) to count their loci without the overlap check")
    skipped_sp = [s for s in focal if s in opt_out]
    if skipped_sp:
        log(f"WARNING: no GFF by choice (--no-gff) for {skipped_sp}; their loci are counted without the overlap check")
    out = []
    for sp in focal:
        if sp not in table:
            continue
        gff, fai = table[sp]
        own = [r for r in coords if tree.canon(r["focal_species"]) == sp and tree.canon(r["species_name"]) == sp]
        genes, skipped, dropped = load_gff(gff, chr_map_from_fai(gff, fai) if fai else None)
        ids, cds_src = set(), Counter()
        for r in own:
            st, g, _, src = own_gene(r["chromosome"], int(r["start"]), int(r["end"]), genes, r["geneID"])
            if st == "named":
                ids.add(g[3])
                cds_src[src] += 1
        res = classify([(r["geneID"], r["chromosome"], int(r["start"]), int(r["end"])) for r in own], genes, ids)
        c = Counter(x[1] for x in res)
        log(f"{sp}: {len(res):,} rows ({skipped:,} non-protein-coding GFF genes ignored, {dropped:,} genes dropped on "
            "chromosomes the .fai could not rename unambiguously); " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
        log(f"  locus CDS from its own transcript: {cds_src['one transcript']:,}; from several transcripts of its gene "
            f"(isoform not in the GFF or gene matched by gene ID): {cds_src['several transcripts']:,}")
        if own and c["unknown"] / len(res) > a.max_unknown:
            ex = [x[0] for x in res if x[1] == "unknown"][:5]
            die(f"{sp}: {c['unknown']:,} of {len(res):,} loci lie on no GFF gene or on a chromosome missing from the GFF "
                f"(e.g. {ex}); the GFF probably belongs to another assembly or names chromosomes differently (give a .fai)")
        checkable = len(res) - c["unknown"]
        if checkable and c["own_gene_not_found"] / checkable > a.max_unmatched:
            ex = [x[0] for x in res if x[1] == "own_gene_not_found"][:5]
            die(f"{sp}: {c['own_gene_not_found']:,} of {checkable:,} loci overlapping GFF genes have no GFF gene carrying "
                f"their geneID (e.g. {ex}); the IDs of the GFF and of the loci probably differ in form")
        out.extend([sp] + list(x) for x in res)
    write_tsv(a.out, ["species", "geneID", "orientation", "other_gene", "chrom", "start", "end"], out)
    log(f"wrote {a.out}: {len(out):,} rows")


if __name__ == "__main__":
    main()
