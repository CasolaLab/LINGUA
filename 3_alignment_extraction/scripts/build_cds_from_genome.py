#!/usr/bin/env python3
"""
build_cds_from_genome.py

Build the reference CDS nucleotide FASTA for stage 3 from the focal genome and its CDS BED, and
check it against the annotated proteins. Exons are cut with `samtools faidx` (BED start taken as
0-based), joined in transcript order (descending start on the minus strand) and reverse-complemented
on the minus strand. Every CDS must translate (standard code, terminal stop removed) to its protein;
a non-ATG start whose first residue the annotation writes as Met (e.g. CUG-initiated POLGARF) is
accepted and reported. Writes nothing unless every transcript passes.

Usage: python3 build_cds_from_genome.py <genome.fa(.fai beside it)> <bed> <proteins.faa> <out.fna> <report.tsv>
  BED: header line, columns chrom, start, end, name, score, strand, transcript ID (= protein ID).
  Protein headers: >ID or >ID__Species.
"""
import collections
import os
import subprocess
import sys

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def die(msg):
    sys.exit(f"ERROR: {msg}")


def read_fasta(path):
    recs, k = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0].split("__")[0]
            if k in recs:
                die(f"duplicate record {k} in {path}")
            recs[k] = []
        elif k:
            recs[k].append(line)
    return {k: "".join(v) for k, v in recs.items()}


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    genome, bed, prot_path, out, report = sys.argv[1:]
    for p in (genome, genome + ".fai", bed, prot_path):
        if not os.path.exists(p):
            die(f"not found: {p}")
    if os.path.exists(out):
        die(f"refusing to overwrite {out}")

    tx = collections.OrderedDict()
    n_rows = 0
    for line in open(bed):
        line = line.rstrip("\r\n")
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 7 or not f[1].isdigit():
            continue
        n_rows += 1
        tx.setdefault(f[6], []).append((f[0], int(f[1]), int(f[2]), f[5]))
    print(f"BED: {n_rows} exon rows, {len(tx)} transcripts", file=sys.stderr)

    regions = []
    for t, ex in tx.items():
        if len({e[0] for e in ex}) != 1 or len({e[3] for e in ex}) != 1:
            die(f"{t} has exons on several chromosomes or strands")
        for c, s, e, _ in ex:
            regions.append(f"{c}:{s + 1}-{e}")
    res = subprocess.run(["samtools", "faidx", genome] + sorted(set(regions)),
                         capture_output=True, text=True)
    if res.returncode != 0:
        die(f"samtools faidx failed: {res.stderr.strip()}")
    seqs = {}
    k = None
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            k = line[1:].strip()
            seqs[k] = []
        else:
            seqs[k].append(line.strip().upper())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    missing = set(regions) - set(seqs)
    if missing:
        die(f"{len(missing)} regions not returned by samtools, e.g. {sorted(missing)[:3]}")

    prots = read_fasta(prot_path)
    rows, cds_out, counts = [], {}, collections.Counter()
    for t, ex in tx.items():
        strand = ex[0][3]
        ex = sorted(ex, key=lambda e: e[1], reverse=(strand == "-"))
        parts = [seqs[f"{c}:{s + 1}-{e}"] for c, s, e, _ in ex]
        cds = "".join(p.translate(RC)[::-1] for p in parts) if strand == "-" else "".join(parts)
        tr = "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))
        notes = []
        if len(cds) % 3:
            notes.append(f"length {len(cds)} not a multiple of 3")
        if tr.endswith("*"):
            tr = tr[:-1]
        else:
            notes.append("no terminal stop")
        p = prots.get(t)
        if p is None:
            status = "no_protein"
        elif tr == p:
            status = "identical"
        elif len(tr) == len(p) and tr[1:] == p[1:] and p[0] == "M" and not cds.startswith("ATG"):
            status = "identical_non_ATG_start"
            notes.append(f"start codon {cds[:3]}")
        else:
            status = "mismatch"
            if len(tr) != len(p):
                notes.append(f"translated {len(tr)} aa vs protein {len(p)} aa")
            else:
                d = [i + 1 for i, (a, b) in enumerate(zip(tr, p)) if a != b]
                notes.append(f"{len(d)} differing residues, first at {d[:5]}")
        counts[status] += 1
        rows.append((t, ex[0][0], strand, len(ex), len(cds), status, "; ".join(notes)))
        cds_out[t] = cds

    with open(report, "w") as fh:
        fh.write("transcript\tchrom\tstrand\tn_exons\tcds_nt\tstatus\tnotes\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    if sum(counts.values()) != len(tx):
        die("status counts do not reconcile")
    print(f"status: {dict(counts)}", file=sys.stderr)
    bad = len(tx) - counts["identical"] - counts["identical_non_ATG_start"]
    if bad:
        die(f"{bad} transcripts fail (see {report}); {out} not written")
    with open(out, "w") as fh:
        for t, s in cds_out.items():
            fh.write(f">{t}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    print(f"{len(tx)} transcripts in, {len(cds_out)} CDS out -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
