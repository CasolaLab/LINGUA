#!/usr/bin/env python3
"""
focal_frame_stops.py  (stage 4, step 5a; method option stops_focal_frame_table)

Rebuilds the stop table of one focal species from the alignments, with every stop codon in the FOCAL reading frame
(2026-11-15). Stage 3 step 13 reads each row in frame 0 of its own gap-free sequence, so in a row that starts off a
focal codon boundary, or after an indel shifts its frame, the stops it reports are not those of the focal frame and
the focal-frame stops are missing (Brassicaceae, 23 focal species: 116,562 focal-frame stops in ancestral rows, 54,762
in the step 13 tables). Here a stop is the triplet of row bases aligned to the three bases of one focal codon, all
three present (no gap), read whatever happens upstream in the row; insertions in the row relative to the focal CDS
are skipped. Every focal codon is read, the focal stop codon included (drop_terminal_stop in step 6b handles it).

Columns are those of the step 13 table: StopCodonStart/End_RefPos are the focal CDS positions of the codon,
cds_length the gap-free length of the row, reference_length the focal CDS length. Chromosome, Start, End and
gff_strand_ref come from the step 13 rows of the same alignment when there are any (else from the file name, strand
NA); the gene ID is the step 13 one when present, else the AncSeq geneID without ".p".

The step 13 table is kept once as <table>.stage3 (no .tsv suffix, so no later step reads it) and is the source of the
metadata on every rerun. A step 5b backup (<table>.pre_start) made from the step 13 table is removed, so step 5b adds
its start-codon rows to the new table.

Usage: python3 focal_frame_stops.py STOP_TABLE ALN_DIR ANCSEQ_TSV --focal-species Arabidopsis_thaliana
"""
import argparse
import collections
import csv
import glob
import os
import shutil
import sys

csv.field_size_limit(10 ** 9)
STOP = {"TAA", "TAG", "TGA"}
FIELDS = ["file_name", "gff_transcriptID", "Chromosome", "Start", "End", "Species", "StopCodonStart_RefPos",
          "StopCodonEnd_RefPos", "StopCodon", "gff_strand_ref", "cds_length", "reference_length"]


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


def focal_frame_stops(ref, row):
    """[(focal position of codon base 1, codon)] for the stop triplets of `row` on the focal codons of `ref`."""
    hc = [i for i, c in enumerate(ref) if c != "-"]
    out = []
    for j in range(0, len(hc) - 2, 3):
        tri = "".join(row[hc[j + x]] for x in range(3))
        if tri in STOP:
            out.append((j + 1, tri))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("table")
    ap.add_argument("aln_dir")
    ap.add_argument("ancseq")
    ap.add_argument("--focal-species", required=True)
    a = ap.parse_args()
    F = a.focal_species
    for p in (a.aln_dir, a.ancseq):
        if not os.path.exists(p):
            die(f"not found: {p}")
    backup = a.table + ".stage3"
    if not os.path.exists(backup):
        if not os.path.exists(a.table):
            die(f"not found: {a.table}")
        pre = a.table + ".pre_start"
        shutil.copy2(pre if os.path.exists(pre) else a.table, backup)  # the step 13 table, without step 5b rows
    old = list(csv.DictReader(open(backup), delimiter="\t"))
    meta = {}
    for r in old:
        meta.setdefault(r["file_name"], r)
    genes = [r["geneID"] for r in csv.DictReader(open(a.ancseq), delimiter="\t")]
    out, n_aln, no_aln, per_kind = [], 0, 0, collections.Counter()
    for g in genes:
        gk = g[:-2] if g.endswith(".p") else g  # Phytozome IDs carry ".p" in the protein tables, not in alignments
        paths = glob.glob(os.path.join(a.aln_dir, f"{F}_{gk}_*.aln"))
        if len(paths) != 1:
            if len(paths) > 1:
                die(f"{len(paths)} alignments for {g} in {a.aln_dir}")
            no_aln += 1
            continue
        n_aln += 1
        fname = os.path.basename(paths[0])[:-4]
        aln = read_aln(paths[0])
        if F not in aln:
            die(f"no {F} row in {paths[0]}")
        ref = aln[F]
        L = len(ref.replace("-", ""))
        m = meta.get(fname)
        if m:
            tid, chrom, start, end, strand = m["gff_transcriptID"], m["Chromosome"], m["Start"], m["End"], m["gff_strand_ref"]
        else:
            chrom, start, end = fname[len(F) + 1 + len(gk) + 1:].rsplit("_", 2)
            tid, strand = gk, "NA"
        for sp, row in aln.items():
            if sp == F:
                continue
            if len(row) != len(ref):
                die(f"row {sp} of {paths[0]} is {len(row)} columns, the focal row {len(ref)}")
            clen = str(len(row.replace("-", "")))
            for pos, cod in focal_frame_stops(ref, row):
                out.append({"file_name": fname, "gff_transcriptID": tid, "Chromosome": chrom, "Start": start, "End": end,
                            "Species": sp, "StopCodonStart_RefPos": str(pos), "StopCodonEnd_RefPos": str(pos + 2),
                            "StopCodon": cod, "gff_strand_ref": strand, "cds_length": clen, "reference_length": str(L)})
                per_kind["ancestral" if sp[:1] == "N" and sp[1:].isdigit() else "extant"] += 1
    with open(a.table, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    pre = a.table + ".pre_start"
    if os.path.exists(pre):
        os.remove(pre)
    print(f"  {F}: {len(out)} focal-frame stops ({per_kind['ancestral']} ancestral, {per_kind['extant']} extant rows) "
          f"in {n_aln} alignments ({no_aln} loci without one); step 13 table ({len(old)} stops) kept as "
          f"{os.path.basename(backup)}", file=sys.stderr)


if __name__ == "__main__":
    main()
