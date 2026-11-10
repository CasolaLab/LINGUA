#!/usr/bin/env python3
"""
step8_stitching_exons.py  (Python port of LINGUA stage 3 step8_stitching_exons.sh, with an optional
one-row-per-species mode that replaces the former step8b_one_row_per_species.py)

Stitches the per-exon alignments of step 4 into one alignment per transcript, in transcript order, and writes
master_summary.tsv. Without --one-row-per-species it reproduces step8_stitching_exons.sh: rows are keyed by the
full header ("species.chromosome"), so a species whose exons come from different chromosomes gets one partial row
per chromosome. With --one-row-per-species, each species keeps only its dominant locus, the header with the most
aligned (non-gap) bases across the gene's exons (ties: the header seen first); exons from its other headers become
gaps. The species set, the most distant species and the MRCA node are unchanged by this; aligned_species then
equals the number of rows. Every species that had more than one row is logged in
<output_dir>/one_row_per_species_log.tsv (file, species, rows_before, kept_row, kept_bases, dropped_rows,
dropped_bases; the same columns as the former step 8b log; dropped rows separated by ';'), and the dropped rows
themselves, in the same columns as the alignment, are written to <output_dir>_one_row_dropped/<file>.aln for
inspection.

run_stage3.sh runs this step WITHOUT --one-row-per-species: step 8c (step8c_one_collinear_row.py) then keeps one
row per species from its dominant collinear chain (one contig and strand, collinear, within a span limit), which
--one-row-per-species does not check, and writes the missing-data mask step 14 uses. Step 8c needs every
species.contig row, so do not give it the output of --one-row-per-species.

Matching the bash script exactly:
  * exon alignments are matched to BED rows by the start_end at the end of the file name and by the chromosome
    appearing (case-insensitively) in the file name;
  * exons are ordered by chromStart, descending on the minus strand;
  * a species missing from an exon gets dashes of that exon's length; a species first seen in a later exon is
    prefixed with dashes for the exons before it;
  * rows are written in order of first appearance, 60 bases per line;
  * distant_species is the species farthest from the reference on the tree (first one on ties); distant_node is
    the MRCA of the species present, from Bio.Phylo, as in the bash script;
  * master_summary.tsv is sorted like `sort -t$'\\t' -k2,2V -k3,3 -k10,10n`.

Usage: TREE_FILE=<newick> python3 step8_stitching_exons.py <step4_output_dir> <bed_file> <output_dir>
                                                           <reference_species> [--one-row-per-species]
"""
import collections
import itertools
import os
import re
import sys

from Bio import Phylo


def sanitize_tid(t):
    return t.replace("\r", "").replace("\n", "").strip()


def sanitize_chr(c):
    c = c.replace("\r", "").replace("\n", "")
    c = re.sub(r"^[Cc]hr", "chr", c)
    c = re.sub(r"\s", "", c)
    while c[:1] in ("-", "_", "."):
        c = c[1:]
    return c


def version_key(s):
    """GNU sort -V approximation: digit runs compare numerically, others lexically."""
    return [(0, int(p), "") if p.isdigit() else (1, 0, p) for p in re.findall(r"\d+|\D+", s)]


def read_aln(path):
    """[(header, seq)] in file order; header is the first word after '>', '\\r' removed from sequences."""
    recs = []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            recs.append([line[1:].split()[0] if line[1:].split() else "", ""])
        elif recs:
            recs[-1][1] += line.replace("\r", "")
    return [(h, s) for h, s in recs if s]


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    one_row = "--one-row-per-species" in sys.argv[1:]
    if len(args) != 4:
        sys.exit(__doc__)
    aln_root, bed, out_root, ref = args
    tree_file = os.environ.get("TREE_FILE")
    if not tree_file or not os.path.exists(tree_file):
        sys.exit("Please set TREE_FILE to a Newick file")
    for p in (aln_root, bed):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
    species_name_clean = ref.replace(" ", "_")
    os.makedirs(out_root, exist_ok=True)
    dropped_dir = out_root.rstrip("/") + "_one_row_dropped"
    if one_row and os.path.isdir(dropped_dir) and os.listdir(dropped_dir):
        sys.exit(f"{dropped_dir} already has files; use a new output folder")

    coord = collections.defaultdict(list)
    for root, _, files in os.walk(aln_root):
        for f in files:
            if f.endswith(".aln"):
                base = f[:-4]
                parts = base.split("_")
                coord[f"{parts[-2]}:{parts[-1]}"].append(os.path.join(root, f))

    t_exons, t_chr, t_strand, validated = collections.defaultdict(list), {}, {}, {}
    with open(bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            # skip '#' comments and a header line (non-numeric start), as prepare_species_inputs.py does; the bash
            # step 8 always skipped the first line, which lost the first exon of a BED without a header
            if line.startswith("#") or len(f) < 7 or not f[1].strip().isdigit():
                continue
            chrom, cs, ce, strand, pid = (f[0].replace("\r", "").replace(" ", ""), f[1].strip(), f[2].strip(),
                                          f[5], f[6].replace("\r", "").replace("\n", "").replace(" ", ""))
            if not chrom or not cs or not pid:
                continue
            paths = coord.get(f"{cs}:{ce}")
            if not paths:
                continue
            cl = chrom.lower()
            aln_path = next((p for p in paths if cl in os.path.basename(p)[:-4].lower()), None)
            if aln_path is None:
                continue
            key = f"{cl}:{cs}:{ce}"
            validated[key] = aln_path
            t_exons[pid].append((int(cs), key))
            t_chr[pid] = chrom
            t_strand[pid] = strand

    tree = Phylo.read(tree_file, "newick")
    tips = {t.name for t in tree.get_terminals() if t.name}
    dist = {}
    for a, b in itertools.combinations([t.name for t in tree.get_terminals() if t.name], 2):
        try:
            d = float(tree.distance(a, b))
        except Exception:
            continue
        dist[(a, b)] = dist[(b, a)] = d

    def mrca(names):
        present = [s for s in names if s in tips]
        if not present:
            return "NA"
        m = tree.common_ancestor(present)
        return str(getattr(m, "name", f"MRCA_unlabeled;clade_size={len(present)}"))   # as in the bash script

    rows_out = []
    log = None
    if one_row:
        log = open(os.path.join(out_root, "one_row_per_species_log.tsv"), "w")
        log.write("file\tspecies\trows_before\tkept_row\tkept_bases\tdropped_rows\tdropped_bases\n")
    for tid, exons in t_exons.items():
        ordered = sorted(exons, key=lambda x: x[0], reverse=(t_strand[tid] == "-"))
        keys = [k for _, k in ordered]
        tid_clean, chr_clean = sanitize_tid(tid), sanitize_chr(t_chr[tid])
        starts = [int(k.split(":")[1]) for k in keys]
        ends = [int(k.split(":")[2]) for k in keys]
        total_len = sum(e - s + 1 for s, e in zip(starts, ends))
        min_start, max_end = min(starts), max(ends)
        out_base = f"{species_name_clean}_{tid_clean}_{chr_clean}_{min_start}_{max_end}"

        exon_recs = [read_aln(validated[k]) for k in keys]
        keep = None
        if one_row:
            bases = collections.Counter()
            first = {}
            for recs in exon_recs:
                for h, s in recs:
                    bases[h] += sum(c != "-" for c in s)
                    first.setdefault(h, len(first))
            best = {}
            for h in first:
                sp = h.split(".")[0]
                if sp not in best or (bases[h], -first[h]) > (bases[best[sp]], -first[best[sp]]):
                    best[sp] = h
            keep = set(best.values())
            for sp, kept in best.items():
                gone = [h for h in first if h.split(".")[0] == sp and h != kept]
                if gone:
                    gone.sort(key=lambda h: first[h])
                    log.write(f"{out_base}.aln\t{sp}\t{len(gone) + 1}\t{kept}\t{bases[kept]}\t"
                              f"{';'.join(gone)}\t{sum(bases[h] for h in gone)}\n")

        # Every row is stitched; with --one-row-per-species only the kept rows go to the alignment and the dropped
        # rows (same columns) go to <output_dir>_one_row_dropped/ (outside <output_dir>, because step 9 searches its
        # input folder recursively). The kept rows are identical to stitching them alone.
        stitched, all_order, stitched_len = {}, [], 0
        for i, recs in enumerate(exon_recs):
            exon_len = max((len(s) for _, s in recs), default=0)
            seqs = dict(recs)
            for h in all_order:
                if h not in seqs:
                    stitched[h] += "-" * exon_len
            for h, _ in recs:
                if h not in stitched:
                    stitched[h] = "-" * stitched_len + seqs[h]
                    all_order.append(h)
                else:
                    stitched[h] += seqs[h]
            stitched_len += exon_len
        order = [h for h in all_order if keep is None or h in keep]
        dropped = [h for h in all_order if keep is not None and h not in keep]

        def write_rows(path, headers):
            with open(path, "w") as o:
                for h in headers:
                    o.write(f">{h}\n")
                    s = stitched[h]
                    for j in range(0, len(s), 60):
                        o.write(s[j:j + 60] + "\n")

        write_rows(os.path.join(out_root, out_base + ".aln"), order)
        if dropped:
            os.makedirs(dropped_dir, exist_ok=True)
            write_rows(os.path.join(dropped_dir, out_base + ".aln"), dropped)

        sp_list = [h.split(".")[0] for h in order]
        if ref in sp_list:
            maxd, far = -1.0, "NA"
            for s in sp_list:
                if s == ref or (ref, s) not in dist:
                    continue
                if dist[(ref, s)] > maxd:
                    maxd, far = dist[(ref, s)], s
            distant = far
        else:
            maxd, pair = -1.0, "NA"
            for i in range(len(sp_list)):
                for j in range(i + 1, len(sp_list)):
                    d = dist.get((sp_list[i], sp_list[j]))
                    if d is not None and d > maxd:
                        maxd, pair = d, f"{sp_list[i]}|{sp_list[j]}"
            distant = pair
        rows_out.append([out_base, tid_clean, chr_clean, t_strand[tid], str(total_len), str(len(keys)),
                         str(len(order)), distant, mrca(sp_list), str(min_start), str(max_end)])

    if log:
        log.close()
    rows_out.sort(key=lambda r: (version_key(r[1]), r[2], int(r[9])))
    with open(os.path.join(out_root, "master_summary.tsv"), "w") as o:
        o.write("file_name\tgff_transcriptID\tchromosome\tstrand\tlength\tnumber_of_exons\taligned_species\t"
                "distant_species\tdistant_node\tstart\tend\n")
        for r in rows_out:
            o.write("\t".join(r) + "\n")
    print(f"Done. {len(rows_out)} transcripts; summary at {os.path.join(out_root, 'master_summary.tsv')}"
          + (" (one row per species)" if one_row else ""))


if __name__ == "__main__":
    main()
