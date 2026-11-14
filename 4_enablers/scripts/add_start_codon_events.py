#!/usr/bin/env python3
"""
add_start_codon_events.py  (stage 4, step 5b; method option start_codon_events, T2T 2026-11-13)

Start-codon gain as an enabler. For every locus whose focal start codon is ATG, every other alignment row (ancestor
or living species) with an aligned non-ATG codon at the focal start-codon columns (no gap, no N) gets one row in the
step 13 stop table, at StopCodonStart_RefPos 1 with StopCodon "noATG_<codon>". Position 1 is in the focal frame and is
never the focal terminal stop, so from here on the event follows every stop rule unchanged: shared by >= 2 consecutive
lineage ancestors, conditional N0, step 6c, synteny, ORF-clade consistency, option B (living carriers = living species
with a non-ATG codon there), step 7g. Any non-ATG codon counts as the same event. A row already present at position 1
for that taxon (a stop codon there) is not duplicated.

The original table is kept once as <table>.pre_start (no .tsv suffix, so step 6 does not read it) and is the input of
every later run, which makes the step idempotent.

Usage: python3 add_start_codon_events.py STOPS_DIR ALN_DIR ANCSEQ_TSV --focal-species Homo_sapiens
"""
import argparse
import csv
import glob
import os
import shutil
import sys

csv.field_size_limit(10 ** 9)


def read_aln(path):
    rows, k = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split(".")[0]
            rows[k] = []
        elif k:
            rows[k].append(line.upper())
    return {a: "".join(b) for a, b in rows.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stops_dir")
    ap.add_argument("aln_dir")
    ap.add_argument("ancseq")
    ap.add_argument("--focal-species", required=True)
    ap.add_argument("--table", default=None, help="step 13 stop table (default <stops_dir>/<F>_DNG_stop_codon_summary.tsv)")
    a = ap.parse_args()
    F = a.focal_species
    table = a.table or os.path.join(a.stops_dir, f"{F}_DNG_stop_codon_summary.tsv")
    base = table + ".pre_start"
    for p in (a.aln_dir, a.ancseq):
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")
    if not os.path.exists(base):
        if not os.path.exists(table):
            sys.exit(f"ERROR: not found: {table}")
        shutil.copy2(table, base)
    rows = list(csv.DictReader(open(base), delimiter="\t"))
    fields = list(csv.DictReader(open(base), delimiter="\t").fieldnames)
    have = {(r["gff_transcriptID"], r["Species"], r["StopCodonStart_RefPos"]) for r in rows}
    genes = [r["geneID"] for r in csv.DictReader(open(a.ancseq), delimiter="\t")]
    added, loci, skipped, noaln = 0, 0, 0, 0
    new = []
    for g in genes:
        g = g[:-2] if g.endswith(".p") else g  # Phytozome IDs carry ".p" in the protein tables, not in alignments
        paths = glob.glob(os.path.join(a.aln_dir, f"{F}_{g}_*.aln"))
        if not paths:
            noaln += 1
            continue
        if len(paths) != 1:
            sys.exit(f"ERROR: {len(paths)} alignments for {g} in {a.aln_dir}")
        aln = read_aln(paths[0])
        ref = aln[F]
        cols = [i for i, c in enumerate(ref) if c != "-"]
        if "".join(ref[c] for c in cols[:3]) != "ATG":
            skipped += 1
            continue
        fname = os.path.basename(paths[0])[:-4]
        tail = fname[len(f"{F}_{g}_"):].rsplit("_", 2)
        chrom, start, end = (tail + ["NA", "NA", "NA"])[:3] if len(tail) == 3 else ("NA", "NA", "NA")
        n_here = 0
        for taxon, row in aln.items():
            if taxon == F:
                continue
            cod = "".join(row[c] for c in cols[:3])
            if cod == "ATG" or set(cod) & set("-N") or len(cod) < 3:
                continue
            if (g, taxon, "1") in have:
                continue
            rec = {k: "NA" for k in fields}
            rec.update({"file_name": fname, "gff_transcriptID": g, "Chromosome": chrom, "Start": start, "End": end,
                        "Species": taxon, "StopCodonStart_RefPos": "1", "StopCodonEnd_RefPos": "3",
                        "StopCodon": f"noATG_{cod}", "reference_length": str(len(cols))})
            new.append(rec)
            n_here += 1
        added += n_here
        loci += n_here > 0
    with open(table, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows + new)
    print(f"  {F}: {added} start-codon rows added in {loci} loci ({skipped} loci without a focal ATG start skipped, {noaln} without an alignment); "
          f"original kept as {os.path.basename(base)}")


if __name__ == "__main__":
    main()
