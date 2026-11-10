"""
step8c_one_collinear_row.py  (LINGUA stage 3, new step 8c: after step 8 stitching, before step 9 trimming)

One row per species and per ancestor, from one place in its genome. Step 8 keys rows by "Species.contig", so exons
extracted from different contigs end up in several rows of one species, and exons from the same contig are merged
into one row however far apart or on whichever strand they lie. This step rebuilds every species' row from the exons
of its dominant collinear chain and turns every other exon of that species into gaps (missing data).

Rule, per (locus, species or ancestor), decided 2026-10-28:
  1. the source of each exon is read from the step 3 MAF of that exon: contig, strand, start, end (minus-strand
     coordinates converted to the plus strand); strand is the species' own strand in the MAF, not the reference's
  2. exons are grouped by (contig, strand); if a species has exons on more than one strand, the strand groups compete
     like any other groups (most aligned bases wins), so all exons are kept only when they share one strand
  3. within a group, exons are chained in reference exon order; a chain breaks when the species' coordinates stop
     moving in one direction (collinearity with the reference lost) or when the chain's genomic span would exceed
     MAXSPAN (50 kb, the limit step 3 applies within an exon)
  4. the dominant chain is the one with the most aligned (non-gap) bases in the stitched alignment; ties go to more
     exons, then the shorter span, then the contig name
  5. the output row, header "Species.<dominant contig>", keeps the dominant chain's exons at their columns and gaps
     all other exons of that species; the reference species is never touched; ancestors follow the same rule
  6. an exon whose MAF blocks for a species lie on both strands is split into one piece per strand (each piece =
     the columns of its blocks' reference positions); the pieces compete like exons, so only the pieces on the
     dominant strand are kept and the others are gapped (decided 2026-10-28)

Missing-data mask (decided 2026-10-28): <out_dir>/step8c_missing_data.tsv (or --mask PATH, also in --dry-run) lists,
per alignment and non-reference species, every reference range in which that species has no aligned sequence after
this step: exons it was absent from (step 8 padding) and every piece gapped here, each with its reason. Stage 3 step
14 reads it to stop calling those gap runs deletions (they are missing data, not disablers).

Checks, each aborting the run: one row per species after; species set unchanged; row lengths equal; every kept base
identical to step 8 at the same column; kept + gapped bases = step 8 bases per species; reference row byte-identical;
every kept row on one contig and strand, collinear, within MAXSPAN; per-exon widths sum to the alignment length;
alignments with nothing to change are copied byte for byte.

Inputs:  <species_dir>: a run_stage3.sh workdir holding <Species>.bed (the only .bed there), step3/, step4/ and
         step8/ (step 8 run WITHOUT --one-row-per-species, so every species.contig row is still present)
Outputs: <out_dir>/*.aln (all alignments), <out_dir>/master_summary.tsv (aligned_species updated), <log_tsv>
         (one line per (alignment, species) changed: kept chain and every gapped exon with its reason)
Usage:   python3 step8c_one_collinear_row.py <species_dir> <out_dir> <log_tsv> [--maxspan 50000] [--dry-run]
                                                [--tree <newick>]
         --dry-run writes the log and the summary only (for sensitivity to --maxspan)
         --tree (default $TREE_FILE) types each logged row living (tree tip), ancestor (labeled internal node) or
         not_in_tree, so any ancestral naming scheme works; without a tree the type is 'unknown (no tree)'
         python3 step8c_one_collinear_row.py --self-test   runs a synthetic case in a temporary folder
"""
import argparse
import collections
import csv
import glob
import os
import re
import shutil
import sys
import tempfile

WRAP = 60


def die(msg):
    sys.exit(f"ERROR: {msg}")


def read_fasta(path):
    recs = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            recs.append([line[1:].split()[0], []])
        elif recs:
            recs[-1][1].append(line.strip())
    return [(h, "".join(s)) for h, s in recs]


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for h, s in recs:
            fh.write(f">{h}\n")
            for i in range(0, len(s), WRAP):
                fh.write(s[i:i + WRAP] + "\n")


def species_of(header):
    return header.split(".", 1)[0]


def contig_of(header):
    return header.split(".", 1)[1] if "." in header else header


def read_bed(path):
    """locus -> (strand, [(chrom, start, end) in BED order])"""
    loci = collections.OrderedDict()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            # skip '#' comments and a header line (non-numeric start) as step 8 does; skipping the first line
            # unconditionally would lose the first exon of a BED without a header
            if line.startswith("#") or len(f) < 7 or not f[1].strip().isdigit():
                continue
            chrom, s, e, strand, tid = f[0].strip(), f[1].strip(), f[2].strip(), f[5].strip(), f[6].strip()
            loci.setdefault(tid, [strand, []])[1].append((chrom, s, e))
    return loci


def tree_node_types(tree_path):
    """name -> 'living' (tip) or 'ancestor' (labeled internal node), read from the Newick tree, so any naming scheme
    for ancestral genomes works (N12, fullTreeAnc235, PrimatesAnc045, ...). None when no tree is given."""
    if not tree_path:
        return None
    from Bio import Phylo
    tree = Phylo.read(tree_path, "newick")
    types = {c.name: "ancestor" for c in tree.get_nonterminals() if c.name}
    types.update({c.name: "living" for c in tree.get_terminals() if c.name})
    if not any(v == "ancestor" for v in types.values()):
        print(f"  warning: {tree_path} has no internal node labels; ancestors will be typed not_in_tree",
              file=sys.stderr)
    return types


def node_type(types, name):
    if types is None:
        return "unknown (no tree)"
    return types.get(name, "not_in_tree")


def exon_key(chrom, s, e):
    return (chrom.lower(), s, e)


def index_exon_files(species_dir):
    """(chrom lower, start, end) -> step 3 MAF path and step 4 per-exon alignment path"""
    maf, aln = {}, {}
    pat = re.compile(r"^[^.]+\.(.+)_(\d+)_(\d+)$")
    for p in glob.glob(os.path.join(species_dir, "step3", "*.txt")):
        m = pat.match(os.path.basename(p)[:-4])
        if m:
            maf[(m.group(1).lower(), m.group(2), m.group(3))] = p
    for d in glob.glob(os.path.join(species_dir, "step4", "*")):
        m = pat.match(os.path.basename(d))
        if m and os.path.isdir(d):
            fs = glob.glob(os.path.join(d, "*.aln"))
            if len(fs) == 1:
                aln[(m.group(1).lower(), m.group(2), m.group(3))] = fs[0]
    return maf, aln


def read_sources(maf_path, ref_species):
    """species -> list of pieces for one exon, one per (contig, strand):
    [contig, strand, sp_start, sp_end, bases, [(ref_start, ref_end), ...]] with the species' coordinates on its plus
    strand and the reference range (plus strand, 0-based, half-open) of every MAF block the piece comes from."""
    pieces = collections.defaultdict(dict)
    block, ref_range = [], None

    def flush():
        nonlocal block, ref_range
        if ref_range is None and block:
            die(f"{maf_path}: a block has no {ref_species} line")
        for sp, contig, start, end, size, strand in block:
            rec = pieces[sp].setdefault((contig, strand), [contig, strand, start, end, 0, []])
            rec[2], rec[3], rec[4] = min(rec[2], start), max(rec[3], end), rec[4] + size
            rec[5].append(ref_range)
        block, ref_range = [], None

    for line in open(maf_path, errors="replace"):
        if line.startswith("a"):
            flush()
            continue
        if not line.startswith("s "):
            continue
        f = line.split()
        sp, contig = species_of(f[1]), contig_of(f[1])
        start, size, strand, srcsize = int(f[2]), int(f[3]), f[4], int(f[5])
        if sp == ref_species:
            if strand != "+":
                die(f"{maf_path}: reference block on the minus strand")
            ref_range = (start, start + size)
            continue
        if strand == "-":
            start = srcsize - start - size
        block.append((sp, contig, start, start + size, size, strand))
    flush()
    return {sp: list(d.values()) for sp, d in pieces.items()}


def exon_width(aln_path):
    return max((len(s) for _, s in read_fasta(aln_path)), default=0)


def chains_for(exons, maxspan):
    """exons: list of (exon_index_in_reference_order, contig, strand, start, end, aligned_bases, exon_id).
    Returns the chains as lists of such exons, grouped by (contig, strand), broken at non-collinear steps or at
    maxspan. Direction: for a '+' strand group coordinates must increase with the reference exon order, for '-'
    they must decrease (reference exon order = ascending reference coordinate)."""
    groups = collections.defaultdict(list)
    for ex in exons:
        groups[(ex[1], ex[2])].append(ex)
    chains = []
    for (contig, strand), exs in groups.items():
        exs.sort(key=lambda x: x[0])
        cur = []
        for ex in exs:
            ok = True
            if cur:
                prev = cur[-1]
                if strand == "+" and ex[3] < prev[4]:
                    ok = False
                if strand == "-" and ex[4] > prev[3]:
                    ok = False
                span = max(ex[4], max(c[4] for c in cur)) - min(ex[3], min(c[3] for c in cur))
                if span > maxspan:
                    ok = False
            if ok:
                cur.append(ex)
            else:
                chains.append((contig, strand, cur))
                cur = [ex]
        if cur:
            chains.append((contig, strand, cur))
    return chains


def dominant(chains):
    def key(ch):
        contig, strand, exs = ch
        span = max(e[4] for e in exs) - min(e[3] for e in exs)
        return (-sum(e[5] for e in exs), -len(exs), span, contig, strand)
    return sorted(chains, key=key)[0]


def gap_reason(ex, dom, maxspan):
    contig, strand, exs = dom
    if ex[1] != contig:
        return "other contig"
    if ex[2] != strand:
        return "opposite strand"
    lo, hi = min(e[3] for e in exs), max(e[4] for e in exs)
    if max(ex[4], hi) - min(ex[3], lo) > maxspan:
        return f"beyond {maxspan} bp"
    return "not collinear"


def process_locus(tid, strand, bed_exons, maf_idx, aln_idx, stitched_recs, ref_species, maxspan, log_rows, fname,
                  mask_rows=None, types=None):
    """Returns (new_recs, changed, stats) for one stitched alignment."""
    order = sorted(range(len(bed_exons)), key=lambda i: int(bed_exons[i][1]), reverse=(strand == "-"))
    ref_order = sorted(range(len(bed_exons)), key=lambda i: int(bed_exons[i][1]))   # ascending reference coordinate
    rank = {i: r for r, i in enumerate(ref_order)}
    # step 8 silently skips an exon without a step 4 alignment, so the same exons are skipped here; the widths-sum
    # check below catches any disagreement
    order = [i for i in order if exon_key(*bed_exons[i]) in aln_idx]
    widths, cols, c = [], {}, 0
    for i in order:
        k = exon_key(*bed_exons[i])
        if k not in maf_idx:
            die(f"{fname}: exon {bed_exons[i]} has a step 4 alignment but no step 3 MAF")
        w = exon_width(aln_idx[k])
        widths.append(w)
        cols[i] = (c, c + w)
        c += w
    L = len(stitched_recs[0][1])
    if c != L:
        die(f"{fname}: per-exon widths sum to {c}, alignment length is {L}")
    if any(len(s) != L for _, s in stitched_recs):
        die(f"{fname}: rows of unequal length")
    rows_by_sp = collections.defaultdict(list)
    for idx, (h, s) in enumerate(stitched_recs):
        rows_by_sp[species_of(h)].append(idx)
    if len(rows_by_sp.get(ref_species, [])) != 1:
        die(f"{fname}: reference species has {len(rows_by_sp.get(ref_species, []))} rows")
    sources = {i: read_sources(maf_idx[exon_key(*bed_exons[i])], ref_species) for i in order}
    ref_row = stitched_recs[rows_by_sp[ref_species][0]][1]
    # reference genomic position of every column of every exon; an insertion column (reference gap) takes the
    # genomically preceding reference position, which is the previous column for a plus-strand gene and the next
    # column for a minus-strand gene (step 4 reverse-complements those exons); the padding base added by step 4
    # sits at the exon's genomic end
    col_pos = {}
    for i in order:
        a, b = cols[i]
        ref_cols = [j for j in range(a, b) if ref_row[j] != "-"]
        start = int(bed_exons[i][1])
        positions = [start + k for k in range(len(ref_cols))]
        if strand == "-":
            positions = positions[::-1]
        pos_at = dict(zip(ref_cols, positions))
        scan = range(a, b) if strand == "+" else range(b - 1, a - 1, -1)
        last = min(positions) if positions else start
        for j in scan:
            if j in pos_at:
                last = pos_at[j]
            col_pos[j] = last

    new_recs, changed = [], False
    stats = collections.Counter()
    for sp, idxs in rows_by_sp.items():
        if sp == ref_species:
            new_recs.append(stitched_recs[idxs[0]])
            continue
        entries = []          # one per (exon, contig, strand) piece: chain fields + the piece's columns and row
        for i in order:
            a, b = cols[i]
            segs = [(idx, stitched_recs[idx][1][a:b]) for idx in idxs]
            with_bases = [(idx, seg) for idx, seg in segs if any(ch != "-" for ch in seg)]
            if not with_bases:
                if mask_rows is not None:
                    mask_rows.append([fname[:-4], sp, f"{bed_exons[i][0]}:{bed_exons[i][1]}-{bed_exons[i][2]}",
                                      bed_exons[i][1], bed_exons[i][2], "absent"])
                continue
            # a species whose only bases in the exon sit at the padding position (beyond the exon end, trimmed by
            # step 9) has no sequence in the exon after trimming: masked as absent, chaining left as it is so the
            # alignments do not change
            if mask_rows is not None:
                exon_end = int(bed_exons[i][2])
                if not any(seg[j - a] != "-" and col_pos[j] < exon_end for idx, seg in with_bases for j in range(a, b)):
                    mask_rows.append([fname[:-4], sp, f"{bed_exons[i][0]}:{bed_exons[i][1]}-{bed_exons[i][2]}",
                                      bed_exons[i][1], bed_exons[i][2], "absent (padding only)"])
            if len(with_bases) > 1:
                die(f"{fname}: {sp} has bases in two rows within exon {bed_exons[i]}")
            idx, seg = with_bases[0]
            pcs = sources[i].get(sp)
            if not pcs:
                die(f"{fname}: {sp} has bases in exon {bed_exons[i]} but no step 3 MAF block")
            if len({pc[0] for pc in pcs}) > 1:
                die(f"{fname}: {sp} on 2 contigs within exon {bed_exons[i]} (step 3 should prevent this)")
            if contig_of(stitched_recs[idx][0]) != pcs[0][0]:
                die(f"{fname}: {sp} exon {bed_exons[i]} is in row {stitched_recs[idx][0]} but its MAF contig is {pcs[0][0]}")
            if len(pcs) > 1:
                stats["exons with pieces on both strands"] += 1
            # assign every column of the exon with a base to the piece whose blocks cover its reference position,
            # else to the piece with the nearest block
            piece_cols = {k: [] for k in range(len(pcs))}
            for j in range(a, b):
                if seg[j - a] == "-":
                    continue
                pos = col_pos[j]
                hit = [k for k, pc in enumerate(pcs) if any(r0 <= pos < r1 for r0, r1 in pc[5])]
                if not hit:
                    hit = [min(range(len(pcs)), key=lambda k: min(min(abs(pos - r0), abs(pos - (r1 - 1)))
                                                                    for r0, r1 in pcs[k][5]))]
                piece_cols[hit[0]].append(j)
            for k, pc in enumerate(pcs):
                if not piece_cols[k]:
                    continue
                entries.append((rank[i], pc[0], pc[1], pc[2], pc[3], len(piece_cols[k]),
                                f"{bed_exons[i][0]}:{bed_exons[i][1]}-{bed_exons[i][2]}", (i, k), idx, piece_cols[k]))
        if not entries:
            die(f"{fname}: {sp} has rows but no bases")
        chains = chains_for(entries, maxspan)
        dom = dominant(chains)
        kept_ids = {e[7] for e in dom[2]}
        row = ["-"] * L
        kept_b = gapped_b = 0
        gapped = []
        for e in entries:
            src_row = stitched_recs[e[8]][1]
            if e[7] in kept_ids:
                for j in e[9]:
                    row[j] = src_row[j]
                kept_b += e[5]
            else:
                gapped_b += e[5]
                reason = gap_reason(e, dom, maxspan)
                gapped.append(f"{e[6]}|{e[1]}|{e[2]}|{reason}")
                if mask_rows is not None:
                    i_ex = e[7][0]
                    pc = sources[i_ex][sp][e[7][1]]
                    whole = len(sources[i_ex][sp]) == 1
                    ranges = [(bed_exons[i_ex][1], bed_exons[i_ex][2])] if whole else sorted(set(pc[5]))
                    for r0, r1 in ranges:
                        mask_rows.append([fname[:-4], sp, e[6], str(r0), str(r1), reason])
        header = f"{sp}.{dom[0]}"
        new_row = "".join(row)
        total_before = sum(sum(ch != "-" for ch in stitched_recs[idx][1]) for idx in idxs)
        if kept_b + gapped_b != total_before:
            die(f"{fname}: {sp} base accounting failed ({kept_b} + {gapped_b} != {total_before})")
        if len(idxs) == 1 and not gapped:
            if new_row != stitched_recs[idxs[0]][1] or header != stitched_recs[idxs[0]][0]:
                die(f"{fname}: {sp} unchanged species would be rewritten differently")
            new_recs.append(stitched_recs[idxs[0]])
            continue
        changed = True
        stats["species rows rebuilt"] += 1
        stats["rows removed"] += len(idxs) - 1
        stats["exon pieces gapped"] += len(gapped)
        stats["bases gapped"] += gapped_b
        for g in gapped:
            stats["gapped: " + g.split("|")[-1]] += 1
        new_recs.append((header, new_row))
        log_rows.append([fname, sp, node_type(types, sp), str(len(idxs)),
                         dom[0], dom[1], str(len(dom[2])), str(len(entries)), str(kept_b), str(gapped_b),
                         ";".join(gapped) or "-"])
    # final checks on the rebuilt alignment
    sp_after = collections.Counter(species_of(h) for h, _ in new_recs)
    if set(sp_after) != set(rows_by_sp) or any(v != 1 for v in sp_after.values()):
        die(f"{fname}: species set or row count wrong after rebuild")
    if any(len(s) != L for _, s in new_recs):
        die(f"{fname}: row length changed")
    before = {species_of(h): 0 for h, _ in stitched_recs}
    for h, s in stitched_recs:
        before[species_of(h)] += sum(ch != "-" for ch in s)
    for sp, idxs in rows_by_sp.items():
        new = next(s for h, s in new_recs if species_of(h) == sp)
        merged = ["-"] * L                       # the species' old rows are disjoint, so they merge into one
        for idx in idxs:
            for j, o in enumerate(stitched_recs[idx][1]):
                if o != "-":
                    if merged[j] != "-":
                        die(f"{fname}: {sp} has bases in two rows at column {j}")
                    merged[j] = o
        if any(n != "-" and n != o for n, o in zip(new, merged)):
            die(f"{fname}: {sp} kept base differs from step 8")
    return new_recs, changed, stats


def run(species_dir, out_dir, log_tsv, maxspan, dry_run, mask_path=None, tree_path=None):
    types = tree_node_types(tree_path)
    beds = glob.glob(os.path.join(species_dir, "*.bed"))
    if len(beds) != 1:
        die(f"expected one BED in {species_dir}, found {len(beds)}")
    ref_species = os.path.basename(beds[0])[:-4]
    step8 = os.path.join(species_dir, "step8")
    for p in (step8, os.path.join(species_dir, "step3"), os.path.join(species_dir, "step4")):
        if not os.path.isdir(p):
            die(f"missing {p}")
    if not dry_run:
        if os.path.abspath(out_dir) == os.path.abspath(step8):
            die("out_dir must differ from step8")
        if os.path.exists(out_dir) and os.listdir(out_dir):
            die(f"{out_dir} is not empty; refusing to overwrite")
        os.makedirs(out_dir, exist_ok=True)
    loci = read_bed(beds[0])
    maf_idx, aln_idx = index_exon_files(species_dir)
    log_rows, totals = [], collections.Counter()
    rows_after = {}
    if mask_path is None and not dry_run:
        mask_path = os.path.join(out_dir, "step8c_missing_data.tsv")
    mask_rows = [] if mask_path else None
    files = sorted(f for f in os.listdir(step8) if f.endswith(".aln"))
    for f in files:
        tid = None
        for cand in loci:
            if f.startswith(f"{ref_species}_{cand}_"):
                tid = cand if tid is None or len(cand) > len(tid) else tid
        if tid is None:
            die(f"no BED locus matches {f}")
        strand, bed_exons = loci[tid]
        recs = read_fasta(os.path.join(step8, f))
        new_recs, changed, stats = process_locus(tid, strand, bed_exons, maf_idx, aln_idx, recs, ref_species,
                                                 maxspan, log_rows, f, mask_rows, types)
        totals.update(stats)
        totals["alignments"] += 1
        totals["alignments changed"] += changed
        rows_after[f[:-4]] = len(new_recs)
        if not dry_run:
            if changed:
                write_fasta(os.path.join(out_dir, f), new_recs)
            else:
                shutil.copy2(os.path.join(step8, f), os.path.join(out_dir, f))
    if mask_rows is not None:
        with open(mask_path, "w") as mk:
            mk.write("alignment\tspecies\texon\tref_start\tref_end\treason\n")
            for r in mask_rows:
                mk.write("\t".join(map(str, r)) + "\n")
        totals["missing-data mask rows"] = len(mask_rows)
    with open(log_tsv, "w") as lg:
        lg.write("alignment\tspecies\ttype\trows_before\tkept_contig\tkept_strand\tkept_pieces\tpieces_with_bases\t"
                 "kept_bases\tgapped_bases\tgapped_pieces(exon|contig|strand|reason)\n")
        for r in log_rows:
            lg.write("\t".join(r) + "\n")
    if not dry_run:
        ms = os.path.join(step8, "master_summary.tsv")
        if os.path.exists(ms):
            with open(ms) as fi, open(os.path.join(out_dir, "master_summary.tsv"), "w") as fo:
                rd = csv.DictReader(fi, delimiter="\t")
                wr = csv.DictWriter(fo, fieldnames=rd.fieldnames, delimiter="\t", lineterminator="\n")
                wr.writeheader()
                for r in rd:
                    if r["file_name"] in rows_after:
                        r["aligned_species"] = str(rows_after[r["file_name"]])
                    wr.writerow(r)
        for f in files:
            if not os.path.exists(os.path.join(out_dir, f)):
                die(f"{f} missing from {out_dir}")
    for k, v in sorted(totals.items()):
        print(f"  {k}: {v:,}", file=sys.stderr)
    return totals


def self_test():
    """Synthetic locus: 3 exons; species A split on two contigs (exons 1-2 on c1, exon 3 on c2, c1 has more bases);
    species B on one contig but exon 3 is 2 Mb away; species C on one contig, exons 1-2 on + and exon 3 on -;
    species D clean. Expected: A keeps c1 and gaps exon 3; B gaps exon 3 (beyond span); C keeps + and gaps exon 3;
    D untouched; species E has exon 2 half on + and half on -: the - half is gapped; reference untouched."""
    d = tempfile.mkdtemp(prefix="step8c_test_")
    ref = "Ref_sp"
    os.makedirs(os.path.join(d, "step3"))
    os.makedirs(os.path.join(d, "step8"))
    exons = [("Chr1", "100", "109"), ("Chr1", "200", "209"), ("Chr1", "300", "309")]
    with open(os.path.join(d, f"{ref}.bed"), "w") as fh:
        fh.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tProteinID\n")
        for c, s, e in exons:
            fh.write(f"{c}\t{s}\t{e}\tT1__CDS\t0\t+\tT1\n")
    seqs = {0: "ACGTACGTAC", 1: "GGGGCCCCAA", 2: "TTTTAAAACCG"}   # exon 3 carries a padding base (11th column)
    src = {  # species -> per exon (contig, start, strand)
        "A": [("c1", 1000, "+"), ("c1", 1100, "+"), ("c2", 5000, "+")],
        "B": [("c9", 1000, "+"), ("c9", 1100, "+"), ("c9", 2_001_100, "+")],
        "C": [("c3", 1000, "+"), ("c3", 1100, "+"), ("c3", 1200, "-")],
        "D": [("c4", 1000, "-"), ("c4", 900, "-"), ("c4", 800, "-")],
        "E": [("c5", 1000, "+"), ("c5", 1100, "+"), ("c5", 1200, "+")],   # exon 2 rewritten below: 5 bp +, 5 bp -
        "F": [("c6", 1000, "+"), ("c6", 1100, "+"), None],                  # absent from exon 3
        "G": [("c7", 1000, "+"), ("c7", 1100, "+"), ("c7", 1200, "+")],   # exon 3: only the padding base (row below)
    }
    for i, (c, s, e) in enumerate(exons):
        base = f"{ref}.{c}_{s}_{e}"
        os.makedirs(os.path.join(d, "step4", base))
        rows = [(f"{ref}.{c}", seqs[i])] + [(f"{sp}.{src[sp][i][0]}", "-" * 10 + "G" if (sp == "G" and i == 2) else seqs[i])
                                             for sp in src if src[sp][i]]
        write_fasta(os.path.join(d, "step4", base, base + ".aln"), rows)
        with open(os.path.join(d, "step3", base + ".txt"), "w") as fh:
            fh.write("##maf version=1\n\na score=0\n")
            fh.write(f"s {ref}.{c} {s} 10 + 1000000 {seqs[i][:10]}\n")
            for sp in src:
                if not src[sp][i]:
                    continue
                ct, st, strand = src[sp][i]
                if sp == "E" and i == 1:
                    continue
                st_field = st if strand == "+" else 10_000_000 - st - 10
                fh.write(f"s {sp}.{ct} {st_field} 10 {strand} 10000000 {seqs[i][:10]}\n")
            if i == 1:
                fh.write(f"\na score=0\ns {ref}.{c} {s} 5 + 1000000 {seqs[i][:5]}\ns E.c5 1100 5 + 10000000 {seqs[i][:5]}\n")
                fh.write(f"\na score=0\ns {ref}.{c} {int(s) + 5} 5 + 1000000 {seqs[i][5:]}\n"
                         f"s E.c5 {10_000_000 - 1105 - 5} 5 - 10000000 {seqs[i][5:]}\n")
    full = seqs[0] + seqs[1] + seqs[2]
    gap = "-" * 11
    stitched = [(f"{ref}.Chr1", full), ("A.c1", seqs[0] + seqs[1] + gap), ("A.c2", "-" * 20 + seqs[2]),
                ("B.c9", full), ("C.c3", full), ("D.c4", full), ("E.c5", full), ("F.c6", seqs[0] + seqs[1] + gap), ("G.c7", seqs[0] + seqs[1] + "-" * 10 + "G")]
    fname = f"{ref}_T1_chr1_100_309"
    write_fasta(os.path.join(d, "step8", fname + ".aln"), stitched)
    with open(os.path.join(d, "step8", "master_summary.tsv"), "w") as fh:
        fh.write("file_name\taligned_species\n" + f"{fname}\t9\n")
    out = os.path.join(d, "step8c")
    totals = run(d, out, os.path.join(d, "log.tsv"), 50000, False)
    got = dict(read_fasta(os.path.join(out, fname + ".aln")))
    assert got[f"{ref}.Chr1"] == full
    assert got["A.c1"] == seqs[0] + seqs[1] + gap and "A.c2" not in got
    assert got["B.c9"] == seqs[0] + seqs[1] + gap
    assert got["C.c3"] == seqs[0] + seqs[1] + gap
    assert got["D.c4"] == full
    assert got["E.c5"] == seqs[0] + seqs[1][:5] + "-" * 5 + seqs[2], got["E.c5"]
    assert totals["rows removed"] == 1 and totals["exon pieces gapped"] == 4
    assert totals["exons with pieces on both strands"] == 1
    reasons = {r.split("|")[-1] for line in open(os.path.join(d, "log.tsv")).read().split("\n")[1:] if line
               for r in line.split("\t")[-1].split(";")}
    assert reasons == {"other contig", "beyond 50000 bp", "opposite strand"}, reasons
    ms = open(os.path.join(out, "master_summary.tsv")).read()
    assert f"{fname}\t8" in ms
    assert got["G.c7"] == seqs[0] + seqs[1] + "-" * 10 + "G"            # alignment untouched
    mask = [l.split("\t") for l in open(os.path.join(out, "step8c_missing_data.tsv")).read().strip().split("\n")[1:]]
    got_mask = {(m[1], m[2], m[3], m[4], m[5]) for m in mask}
    assert ("F", "Chr1:300-309", "300", "309", "absent") in got_mask, got_mask
    assert ("A", "Chr1:300-309", "300", "309", "other contig") in got_mask
    assert ("B", "Chr1:300-309", "300", "309", "beyond 50000 bp") in got_mask
    assert ("C", "Chr1:300-309", "300", "309", "opposite strand") in got_mask
    assert ("E", "Chr1:200-209", "205", "210", "opposite strand") in got_mask, got_mask
    assert ("G", "Chr1:300-309", "300", "309", "absent (padding only)") in got_mask, got_mask
    assert len(mask) == 6, mask
    shutil.rmtree(d)
    print("  self-test passed", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("species_dir", nargs="?")
    ap.add_argument("out_dir", nargs="?")
    ap.add_argument("log_tsv", nargs="?")
    ap.add_argument("--maxspan", type=int, default=50000)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--mask", default=None, help="missing-data mask path (default <out_dir>/step8c_missing_data.tsv; "
                                                 "in --dry-run written only when given)")
    ap.add_argument("--tree", default=os.environ.get("TREE_FILE") or None,
                    help="Newick tree (default $TREE_FILE): tips are typed living and labeled internal nodes ancestor "
                         "in the log's type column; without a tree the type is 'unknown (no tree)'")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return
    if not (a.species_dir and a.out_dir and a.log_tsv):
        ap.error("species_dir, out_dir and log_tsv are required")
    if a.tree and not os.path.isfile(a.tree):
        ap.error(f"tree not found: {a.tree}")
    run(a.species_dir, a.out_dir, a.log_tsv, a.maxspan, a.dry_run, a.mask, a.tree)


if __name__ == "__main__":
    main()
