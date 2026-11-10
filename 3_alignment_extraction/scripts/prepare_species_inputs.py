#!/usr/bin/env python3
"""
prepare_species_inputs.py

Writes the exon list for step 1 and the reference CDS FASTA for one species, as the 2026-10-06 stage 3 runner did:
unique exon coordinates in BED order, and one CDS record per BED transcript in BED order (Phytozome CDS IDs that
carry a '.p' suffix the BED IDs lack are matched without it). Stops if a transcript has no CDS.

Usage: python3 prepare_species_inputs.py <bed> <cds_fasta> <coords_out> <n_exons_out> <ref_raw_out>
"""
import sys


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    bed, cds_src, coords_out, n_out, ref_out = sys.argv[1:]
    exons, tx = [], []
    seen_e, seen_t = set(), set()
    for line in open(bed):
        line = line.rstrip("\n")
        if line.startswith("#") or not line:
            continue
        f = line.split("\t")
        if len(f) < 7 or not f[1].strip().isdigit():   # header line without '#'
            continue
        if (f[0], f[1], f[2]) not in seen_e:
            seen_e.add((f[0], f[1], f[2]))
            exons.append((f[0], f[1], f[2]))
        if f[6] not in seen_t:
            seen_t.add(f[6])
            tx.append(f[6])
    with open(coords_out, "w") as out:
        out.write("seq\tstart\tstop\n")
        for e in exons:
            out.write("\t".join(e) + "\n")
    open(n_out, "w").write(str(len(exons)))
    seqs, k = {}, None
    for line in open(cds_src):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]
            seqs[k] = []
        elif k:
            seqs[k].append(line)
    for t in tx:
        if t not in seqs and t + ".p" in seqs:
            seqs[t] = seqs[t + ".p"]
    absent = [t for t in tx if t not in seqs]
    if absent:
        sys.exit(f"{len(absent)} transcripts missing from the CDS FASTA, e.g. {absent[:3]}")
    with open(ref_out, "w") as out:
        for t in tx:
            out.write(f">{t}\n{''.join(seqs[t])}\n")
    print(f"prepared: {len(exons)} exons for {len(tx)} transcripts")


if __name__ == "__main__":
    main()
