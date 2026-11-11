#!/usr/bin/env python3
"""
step14b_deletion_evidence.py  (runs after step 14; reads its output, changes nothing in step14/)

For every deletion step 14 reported, records where the gap came from and what the alignment says about whether it is
a real deletion, in one table outside step14/ (stage 4 reads every .tsv in step14/ as indel calls).

Where the gap came from (gap_source), from the step 1 and step 3 MAFs of each exon the deletion covers:
  cactus_gap      every position is a '-' inside a block step 3 kept: the aligner's own gap
  no_alignment    some position has no block of the species at all, even in step 1
  dropped_piece   some position was covered only by a block step 3 (dominant chain) dropped from the species
For no_alignment and dropped_piece gaps within one exon, the in-place check compares the kept blocks just before and
after the gap: reference_gap_bp = reference bases between them, species_bases_in_place = species bases between them
(on the species' own MAF strand). sides_touch = yes when the species has no bases there (a real deletion looks like
this); in_place_difference_bp = reference_gap_bp - species_bases_in_place, the likely real size of the change;
in_place_frameshift = yes when that difference is not a multiple of 3.

evidence:
  aligner gap (kept as called)   cactus_gap
  likely real deletion           sides_touch = yes
  possible smaller frameshift    sides do not touch, in_place_difference_bp is not a multiple of 3
  probably not a frameshift      sides do not touch, in_place_difference_bp is 0 or a multiple of 3
  can't tell                     the in-place check could not run (cant_tell_reason says why)

Usage:
  python3 step14b_deletion_evidence.py <step14_master_indels.tsv> <step1_dir> <step3_dir> <step9_dir> <bed_file>
                                       <reference_species> <output.tsv> [--tree <newick>]
  --tree (default $TREE_FILE): species_type (living / ancestor / not_in_tree) and the sharing columns; without a
  tree those columns are NA.
"""
import argparse
import collections
import csv
import os
import sys


def die(msg):
    sys.exit(f"ERROR: {msg}")


def read_maf(path, ref):
    """[(ref_start, ref_len, [(sp, contig, start, size, strand, src_size, seq), ...])]"""
    blocks, cur = [], None
    for line in open(path):
        if line.startswith("a"):
            cur = [None, None, []]
            blocks.append(cur)
        elif line.startswith("s ") and cur is not None:
            f = line.split()
            sp, contig = f[1].split(".", 1) if "." in f[1] else (f[1], "NA")
            if sp == ref:
                if "-" in f[6] or f[4] != "+":
                    die(f"{path}: reference row with gaps or on the minus strand")
                cur[0], cur[1] = int(f[2]), len(f[6])
            else:
                cur[2].append((sp, contig, int(f[2]), int(f[3]), f[4], int(f[5]), f[6]))
    return [b for b in blocks if b[0] is not None]


class Exon:
    """Per species: status of each reference position (base / gap_in_block / dropped; absent when missing), kept and
    dropped blocks as (ref_start, ref_end, contig, strand, maf_start, maf_end, plus_start, plus_end)."""
    RANK = {"base": 3, "gap_in_block": 2, "dropped": 1}

    def __init__(self, step1_path, step3_path, ref):
        b3 = read_maf(step3_path, ref) if os.path.exists(step3_path) else []
        kept_lines = collections.Counter((b[0], b[1], r) for b in b3 for r in b[2])
        self.status = collections.defaultdict(dict)
        self.kept = collections.defaultdict(list)
        self.dropped = collections.defaultdict(list)
        for ref_start, ref_len, recs in read_maf(step1_path, ref):
            for r in recs:
                key = (ref_start, ref_len, r)
                is_kept = kept_lines[key] > 0
                if is_kept:
                    kept_lines[key] -= 1
                sp, contig, start, size, strand, src, seq = r
                for i, ch in enumerate(seq):
                    st = ("base" if ch != "-" else "gap_in_block") if is_kept else "dropped"
                    old = self.status[sp].get(ref_start + i)
                    if old is None or self.RANK[st] > self.RANK[old]:
                        self.status[sp][ref_start + i] = st
                plus = (start, start + size) if strand == "+" else (src - start - size, src - start)
                blk = (ref_start, ref_start + ref_len, contig, strand, start, start + size) + plus
                (self.kept if is_kept else self.dropped)[sp].append(blk)
        if sum(kept_lines.values()):
            die(f"{step3_path}: {sum(kept_lines.values())} step 3 block lines not found in step 1")
        for v in self.kept.values():
            v.sort()


def tree_relations(path):
    if not path:
        return None, None, None
    from Bio import Phylo
    tree = Phylo.read(path, "newick")
    types, parent, children = {}, {}, {}
    for cl in tree.find_clades(order="level"):
        if not cl.name:
            continue
        types[cl.name] = "living" if cl.is_terminal() else "ancestor"
        children[cl.name] = [c.name for c in cl.clades if c.name]
        for c in cl.clades:
            if c.name:
                parent[c.name] = cl.name
    return types, parent, children


def main():
    ap = argparse.ArgumentParser(description="Step 14b: evidence for each step 14 deletion")
    for name in ("master_indels", "step1_dir", "step3_dir", "step9_dir", "bed_file", "reference_species", "out_tsv"):
        ap.add_argument(name)
    ap.add_argument("--tree", default=os.environ.get("TREE_FILE") or None)
    a = ap.parse_args()
    for p in (a.master_indels, a.bed_file):
        if not os.path.isfile(p):
            die(f"missing file {p}")
    for p in (a.step1_dir, a.step3_dir, a.step9_dir):
        if not os.path.isdir(p):
            die(f"missing folder {p}")
    out_dir = os.path.dirname(os.path.abspath(a.out_tsv))
    if os.path.basename(out_dir) == "step14":
        die("write the output outside step14/ (stage 4 reads every .tsv there as indel calls)")
    ref = a.reference_species
    types, parent, children = tree_relations(a.tree)

    exons, chrom_of, strand_of = collections.defaultdict(list), {}, {}
    for line in open(a.bed_file):
        f = line.rstrip("\r\n").split("\t")
        if line.startswith("#") or len(f) < 7 or not f[1].strip().isdigit():
            continue
        tid = f[6].strip()
        exons[tid].append((int(f[1]), int(f[2])))
        chrom_of[tid], strand_of[tid] = f[0].strip(), f[5].strip()
    order = {t: sorted(v, reverse=(strand_of[t] == "-")) for t, v in exons.items()}

    def cds_to_genomic(tid, p):
        cum = 0
        for s0, e0 in order[tid]:
            if p <= cum + e0 - s0:
                return (s0, e0), (s0 + p - cum - 1 if strand_of[tid] == "+" else e0 - (p - cum))
            cum += e0 - s0
        die(f"{tid}: CDS position {p} is beyond the CDS")

    cache = {}

    def exon_data(tid, ex):
        k = (tid, ex)
        if k not in cache:
            name = f"{ref}.{chrom_of[tid]}_{ex[0]}_{ex[1]}.txt"
            p1 = os.path.join(a.step1_dir, name)
            if not os.path.exists(p1):
                die(f"no step 1 MAF {p1} for a reported deletion")
            cache[k] = Exon(p1, os.path.join(a.step3_dir, name), ref)
        return cache[k]

    dels = list(csv.DictReader(open(a.master_indels), delimiter="\t"))
    by_aln_sp = collections.defaultdict(list)
    for d in dels:
        by_aln_sp[(d["file_name"], d["Species"])].append((int(d["IndelStart_RefPos"]), int(d["IndelEnd_RefPos"])))
    present = {}
    for f in os.listdir(a.step9_dir):
        if f.endswith(".aln"):
            present[f[:-4]] = {l[1:].split(".")[0] for l in open(os.path.join(a.step9_dir, f)) if l.startswith(">")}

    cols = ["file_name", "gff_transcriptID", "Species", "species_type", "IndelStart_RefPos", "IndelEnd_RefPos",
            "IndelLength", "gap_source", "bases_cactus_gap", "bases_no_alignment", "dropped_covers_bp",
            "dropped_from", "dropped_distance_bp", "reference_gap_bp", "species_bases_in_place", "sides_touch",
            "in_place_difference_bp", "in_place_frameshift", "exon_edge", "cant_tell_reason",
            "shared_with_parent", "shared_with_sister", "evidence"]
    rows, stats = [], collections.Counter()
    for d in dels:
        aln, sp, tid = d["file_name"], d["Species"], d["gff_transcriptID"]
        a0, b0 = int(d["IndelStart_RefPos"]), int(d["IndelEnd_RefPos"])
        if tid not in order:
            die(f"{aln}: transcript {tid} is not in the BED")
        st = collections.Counter()
        by_exon = collections.defaultdict(list)
        for p in range(a0, b0 + 1):
            ex, g = cds_to_genomic(tid, p)
            by_exon[ex].append(g)
            s = exon_data(tid, ex).status[sp].get(g, "absent")
            st[s] += 1
        if st["base"]:
            die(f"{aln} {sp} {a0}-{b0}: {st['base']} deletion positions have a base in a kept block")
        src = "dropped_piece" if st["dropped"] else "no_alignment" if st["absent"] else "cactus_gap"
        r = dict.fromkeys(cols, "")
        r.update({"file_name": aln, "gff_transcriptID": tid, "Species": sp,
                  "species_type": "NA" if types is None else types.get(sp, "not_in_tree"),
                  "IndelStart_RefPos": str(a0), "IndelEnd_RefPos": str(b0), "IndelLength": d["IndelLength"],
                  "gap_source": src, "bases_cactus_gap": str(st["gap_in_block"]),
                  "bases_no_alignment": str(st["absent"]), "dropped_covers_bp": str(st["dropped"])})

        # exon edge, in gene direction
        if len(by_exon) > 1:
            r["exon_edge"] = f"spans {len(by_exon)} exons"
        else:
            (s0, e0), gs = next(iter(by_exon.items()))
            first, last = (s0, e0 - 1) if strand_of[tid] == "+" else (e0 - 1, s0)
            hit_first, hit_last = first in gs, last in gs
            r["exon_edge"] = ("both" if hit_first and hit_last else "start" if hit_first else
                              "end" if hit_last else "none")

        if src != "cactus_gap":
            if len(by_exon) == 1:
                ex, gs = next(iter(by_exon.items()))
                ed = exon_data(tid, ex)
                kept = ed.kept.get(sp, [])
                g0, g1 = min(gs), max(gs)
                if src == "dropped_piece":
                    drop = [k for k in ed.dropped.get(sp, []) if any(k[0] <= g < k[1] for g in gs)]
                    origin, dists = set(), []
                    for k in drop:
                        if not kept:
                            origin.add("no kept piece in exon")
                            continue
                        same = [x for x in kept if x[2] == k[2]]
                        if not same:
                            origin.add("other chromosome")
                            continue
                        origin.add("same chromosome and strand" if any(x[3] == k[3] for x in same)
                                   else "same chromosome, other strand")
                        lo, hi = min(x[6] for x in same), max(x[7] for x in same)
                        dists.append(0 if (k[6] < hi and k[7] > lo) else min(abs(k[6] - hi), abs(lo - k[7])))
                    r["dropped_from"] = ";".join(sorted(origin))
                    r["dropped_distance_bp"] = str(min(dists)) if dists else ""
                left = max((k for k in kept if k[0] < g0), key=lambda k: k[0], default=None)
                right = min((k for k in kept if k[1] > g1 + 1 and (left is None or k[0] > left[0])),
                            key=lambda k: k[0], default=None)
                minus = strand_of[tid] == "-"
                if not kept:
                    reason = "no kept piece anywhere in this exon"
                elif left is None:
                    reason = "no kept piece after the gap in this exon" if minus else \
                        "no kept piece before the gap in this exon"
                elif right is None:
                    reason = "no kept piece before the gap in this exon" if minus else \
                        "no kept piece after the gap in this exon"
                elif (left[2], left[3]) != (right[2], right[3]):
                    reason = "the two sides are on different chromosomes or strands"
                else:
                    ref_gap, sp_gap = right[0] - left[1], right[4] - left[5]
                    reason = "the species' two sides overlap" if sp_gap < 0 else ""
                    if not reason:
                        diff = ref_gap - sp_gap
                        r.update({"reference_gap_bp": str(ref_gap), "species_bases_in_place": str(sp_gap),
                                  "sides_touch": "yes" if sp_gap == 0 else "no",
                                  "in_place_difference_bp": str(diff),
                                  "in_place_frameshift": "yes" if diff % 3 else "no"})
                        r["evidence"] = ("likely real deletion" if sp_gap == 0 else
                                         "possible smaller frameshift" if diff % 3 else "probably not a frameshift")
                r["cant_tell_reason"] = reason
            else:
                r["cant_tell_reason"] = "gap crosses into the next exon"
                if src == "dropped_piece":
                    r["dropped_from"] = "not assessed (gap spans exons)"
            if not r["evidence"]:
                r["evidence"] = "can't tell"
        else:
            r["evidence"] = "aligner gap (kept as called)"

        if types is None:
            r["shared_with_parent"] = r["shared_with_sister"] = "NA"
        else:
            pres = present.get(aln, set())
            par = parent.get(sp)

            def shared(other):
                return any(x <= b0 and y >= a0 for x, y in by_aln_sp.get((aln, other), []))
            r["shared_with_parent"] = ("no parent" if par is None else "no data" if par not in pres else
                                       "yes" if shared(par) else "no")
            sis = [s for s in children.get(par, []) if s != sp and s in pres] if par else []
            r["shared_with_sister"] = "no data" if not sis else "yes" if any(shared(s) for s in sis) else "no"
        rows.append(r)
        stats[(src, r["evidence"])] += 1

    if sum(stats.values()) != len(dels):
        die("rows do not reconcile with step 14")
    with open(a.out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in cols) + "\n")
    print(f"step 14b: {len(dels):,} deletions from {a.master_indels}")
    for (src, ev), n in sorted(stats.items()):
        print(f"   {src:<14} {ev:<30} {n:,}")
    print(f"   written to {a.out_tsv}")


if __name__ == "__main__":
    main()
