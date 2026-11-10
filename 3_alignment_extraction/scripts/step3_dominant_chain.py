#!/usr/bin/env python3
"""
step3_dominant_chain.py  (Option A replacement for LINGUA stage 3 step3_filtering_species.py)

Per exon (one mafExtractor file), for every non-reference species (ancestral nodes included):
  1. group its blocks by (chromosome, strand);
  2. within each group, find the best collinear chain: blocks taken in MAF order (the reference order) whose species
     coordinates also increase, without overlapping, with at most JUMP bp between consecutive blocks
     (same distance measure as the repo step 3), maximizing the number of aligned bases (dynamic programming);
  3. keep the group whose best chain has the most aligned bases, and only that chain;
  4. drop that species' sequence from every other block (other chromosomes or strands, distant, out-of-order or
     overlapping blocks).
The species is kept with one collinear orthologous chain, instead of being removed from the whole exon (strict) or
keeping every block (relaxed). Ties are broken deterministically (earliest blocks, first group seen). Reference-species blocks are never touched; blocks that lose the reference are not
written (as in the repo step 3). Output MAFs have the same format, so steps 4-17 run unchanged.

Usage: python3 step3_dominant_chain.py -i <mafextractor_dir> -o <output_dir> -r <reference_species>
                                       -s <summary_csv> [-j 50000] [-t 8]
Summary CSV: file_name,species,status,blocks_kept,blocks_dropped,bases_kept,bases_dropped
  status = clean (nothing dropped) | trimmed (outlier blocks dropped)
"""
import argparse
import collections
import concurrent.futures
import os
import sys


def parse_blocks(path):
    blocks, cur = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("a"):
                if cur is not None:
                    blocks.append(cur)
                cur = [line.rstrip("\n")]
            elif line.startswith("s") and cur is not None:
                cur.append(line.rstrip("\n"))
    if cur is not None:
        blocks.append(cur)
    return blocks


def s_fields(line):
    p = line.split()
    full = p[1]
    sp, chrom = (full.split(".", 1) + ["unknown"])[:2] if "." in full else (full, "unknown")
    start, size, strand, seq = int(p[2]), int(p[3]), p[4], p[6]
    return sp, chrom, start, size, strand, seq


def best_collinear_chain(g, jump):
    """g: blocks of one (chromosome, strand) in MAF order, as (block index, chrom, strand, start, end, bases).
    Weighted longest chain with start[i] >= end[j] and start[i] - end[j] <= jump; returns the blocks in order."""
    n = len(g)
    score = [g[i][5] for i in range(n)]
    prev = [-1] * n
    for i in range(n):
        for j in range(i):
            gap = g[i][3] - g[j][4]
            if 0 <= gap <= jump and score[j] + g[i][5] > score[i]:
                score[i] = score[j] + g[i][5]
                prev[i] = j
    i = max(range(n), key=lambda k: (score[k], -k))
    chain = []
    while i != -1:
        chain.append(g[i])
        i = prev[i]
    return chain[::-1]


def process(args):
    in_path, out_dir, ref, jump = args
    base = os.path.basename(in_path)
    blocks = parse_blocks(in_path)
    occ = collections.defaultdict(list)          # species -> [(block index, chrom, strand, start, end, bases)]
    for bi, b in enumerate(blocks):
        for line in b[1:]:
            sp, chrom, start, size, strand, seq = s_fields(line)
            if sp != ref:
                occ[sp].append((bi, chrom, strand, start, start + size, sum(c != "-" for c in seq)))
    keep = {}                                    # species -> set of block indices to keep
    rows = []
    for sp, lst in occ.items():
        by_group = collections.defaultdict(list)
        for o in lst:
            by_group[(o[1], o[2])].append(o)
        best = max((best_collinear_chain(g, jump) for g in by_group.values()),
                   key=lambda c: sum(o[5] for o in c))
        keep[sp] = {o[0] for o in best}
        kb = sum(o[5] for o in best)
        tb = sum(o[5] for o in lst)
        status = "clean" if len(best) == len(lst) else "trimmed"
        rows.append(f"{base},{sp},{status},{len(best)},{len(lst) - len(best)},{kb},{tb - kb}")
    out = ["##maf version=1", ""]
    for bi, b in enumerate(blocks):
        lines = [b[0]]
        has_ref = False
        for line in b[1:]:
            sp = s_fields(line)[0]
            if sp == ref:
                lines.append(line)
                has_ref = True
            elif bi in keep.get(sp, ()):
                lines.append(line)
        if has_ref:
            out.extend(lines + [""])
    with open(os.path.join(out_dir, base), "w") as fh:
        fh.write("\n".join(out) + "\n")
    return rows or [f"{base},None,clean,0,0,0,0"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("-i", "--mafextractor_dir", required=True)
    ap.add_argument("-o", "--output_dir", required=True)
    ap.add_argument("-r", "--ref", required=True)
    ap.add_argument("-s", "--summary", required=True)
    ap.add_argument("-j", "--jump", type=int, default=50000)
    ap.add_argument("-t", "--threads", type=int, default=8)
    a = ap.parse_args()
    if not os.path.isdir(a.mafextractor_dir):
        sys.exit(f"missing {a.mafextractor_dir}")
    os.makedirs(a.output_dir, exist_ok=True)
    files = sorted(f for f in os.listdir(a.mafextractor_dir) if f.endswith(".txt"))
    tasks = [(os.path.join(a.mafextractor_dir, f), a.output_dir, a.ref, a.jump) for f in files]
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.threads) as ex, open(a.summary, "w") as out:
        out.write("file_name,species,status,blocks_kept,blocks_dropped,bases_kept,bases_dropped\n")
        n = collections.Counter()
        for rows in ex.map(process, tasks, chunksize=20):
            for r in rows:
                out.write(r + "\n")
                n[r.split(",")[2]] += 1
    print(f"step3_dominant_chain: {len(files)} files; species-exon entries: {dict(n)}")


if __name__ == "__main__":
    main()
