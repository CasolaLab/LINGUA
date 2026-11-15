#!/usr/bin/env python3
"""
drop_unaligned_end_indels.py  (stage 4, step 5c; method options indel_unaligned_end_from_alignment,
                               indel_mask_remainder_from_alignment)

Removes from the stage 3 indel table of one focal species every row (one species or ancestor at one indel) whose gap,
in the alignment that stage 4 reads, belongs to a gap run reaching the first or last base of the focal CDS
(2026-11-15). Such a row is an unaligned end of that sequence, missing data, not a deletion. Stage 3 step 14 trims a
deletion that overlaps a missing-data exon to its part outside the exon (rule of 2026-10-29), so where the gap is an
unaligned end the reported remainder starts at an exon boundary and the step 6b terminal-gap test, which reads the
table coordinates, keeps it (Brassicaceae AT5G03710.1: table deletion 14-78 of 246, ancestral rows unaligned over
1-78). Removing the row, not the whole event, keeps the event at the carriers whose gap is internal, and removes the
same artefact where it would count as extant support.

indel_mask_remainder_from_alignment also removes a row whose gap run in the alignment is internal but longer than the
table interval: the only way step 14 reports less than the whole run is the missing-data trim, so the row is the
remainder of a gap that overlaps missing data (Brassicaceae Iin25571.t1: table deletion 391-541, the ancestral rows
lack about 210 positions around it). This overrides, for stage 4, the stage 3 rule of 2026-10-29 that counts such a
remainder as a deletion.

The gap run is read on the focal CDS positions of the row (row bases in focal-gap columns do not break it). Rows whose
interval is not a gap in the alignment (insertion rows) are kept unchanged.

The original table is kept once as <table>.pre_unaligned (no .tsv suffix, so step 6 does not read it) and is the input
of every rerun; the removed rows go to <table>.unaligned_end_rows (TSV, also without the .tsv suffix).

Usage: python3 drop_unaligned_end_indels.py INDEL_TABLE ALN_DIR --focal-species Arabidopsis_thaliana
"""
import argparse
import collections
import csv
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lingua_config as cfg  # noqa: E402

csv.field_size_limit(10 ** 9)


def die(m):
    sys.exit(f"ERROR: {m}")


def read_aln(p):
    rows, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split(".")[0]
            rows[k] = []
        elif k:
            rows[k].append(line.upper())
    return {a: "".join(b) for a, b in rows.items()}


def unaligned_end(ref, row, start, end):
    """"unaligned_end" when focal positions start..end (1-based) are a gap in `row` that is part of a gap run reaching
    CDS position 1 or L, "mask_remainder" when the run is internal but extends beyond start..end, "" when the run is
    exactly start..end, None when the interval is not a gap in the row."""
    hc = [c for c, ch in enumerate(ref) if ch != "-"]
    L = len(hc)
    if not (1 <= start <= end <= L) or any(row[hc[i]] != "-" for i in range(start - 1, end)):
        return None
    s, e = start - 1, end - 1
    while s > 0 and row[hc[s - 1]] == "-":
        s -= 1
    while e < L - 1 and row[hc[e + 1]] == "-":
        e += 1
    if s == 0 or e == L - 1:
        return "unaligned_end"
    return "mask_remainder" if (s, e) != (start - 1, end - 1) else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("table")
    ap.add_argument("aln_dir")
    ap.add_argument("--focal-species", required=True)
    a = ap.parse_args()
    F = a.focal_species
    if not os.path.isdir(a.aln_dir):
        die(f"not found: {a.aln_dir}")
    base = a.table + ".pre_unaligned"
    if not os.path.exists(base):
        if not os.path.exists(a.table):
            die(f"not found: {a.table}")
        shutil.copy2(a.table, base)
    rd = csv.DictReader(open(base), delimiter="\t")
    fields = rd.fieldnames
    rows = list(rd)
    cache = {}
    on = {"unaligned_end": bool(cfg.opt("indel_unaligned_end_from_alignment")),
          "mask_remainder": bool(cfg.opt("indel_mask_remainder_from_alignment"))}
    kept, removed, n_notgap = [], [], 0
    per = collections.Counter()
    for r in rows:
        fname = r["file_name"]
        if fname not in cache:
            p = os.path.join(a.aln_dir, fname + ".aln")
            if not os.path.exists(p):
                hits = glob.glob(os.path.join(a.aln_dir, f"{F}_{r['gff_transcriptID']}_*.aln"))
                p = hits[0] if len(hits) == 1 else None
            cache[fname] = read_aln(p) if p else None
        aln = cache[fname]
        if aln is None or F not in aln or r["Species"] not in aln:
            die(f"no alignment row for {r['Species']} of {fname} in {a.aln_dir}")
        u = unaligned_end(aln[F], aln[r["Species"]], int(r["IndelStart_RefPos"]), int(r["IndelEnd_RefPos"]))
        if u and on[u]:
            removed.append({**r, "reason": u})
            per[u] += 1
        else:
            n_notgap += u is None
            kept.append(r)
    for path, out, fl in ((a.table, kept, fields), (a.table + ".unaligned_end_rows", removed, fields + ["reason"])):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fl, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(out)
    assert len(kept) + len(removed) == len(rows)
    print(f"  {F}: {len(rows)} indel rows, {len(removed)} removed ({per['unaligned_end']} unaligned ends, "
          f"{per['mask_remainder']} mask remainders), {n_notgap} not a gap in the alignment (kept); original kept as "
          f"{os.path.basename(base)}", file=sys.stderr)


if __name__ == "__main__":
    main()
