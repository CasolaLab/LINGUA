#!/usr/bin/env python3
"""
remove_literature_etes.py  (stage 4, step 7e; T2T mammals only; method option remove_literature_etes)

Remove candidates (DNG and PDNG) that are exapted TE genes (ETEs) described in the literature, e.g. suppressyn (ERVH48-1).
TE overlap at the DNA level is NOT a reason for removal (his rule, 2026-11-12): many de novo genes arise inside TEs; a candidate
goes only if its protein matches a TE protein (stage 2, 50/50 rule) or is a known ETE (this step).

ETE proteins: data/ete/human_ETE_proteins_20261112.faa (header ">protein_id SYMBOL"), every CHM13 RefSeq isoform of the 213 genes of
data/ete/human_ETE_consensus_catalog.tsv plus UniProt/Ensembl sequences for the 17 catalog entries with no CHM13 protein.
Each focal species' candidate proteins are searched with BLASTp (-seg no, E <= 1e-5, max_target_seqs 1000, as the stage 2
paralog rule); a candidate is an ETE when, for one ETE protein, identity >= 50% (identical positions / aligned length over all
HSPs) over >= 50% of the candidate OR of the ETE protein (union of the HSP intervals on that protein; the ETE side catches
orthologs with extra sequence, e.g. the 297-aa M. fascicularis suppressyn, 89.7% identical over 90.6% of human ERVH48-1 but only
48.8% of itself). Orthologs in other mammals are caught this way when they keep >= 50% identity to the human ETE.

Reads and rewrites <out_prefix>_nodes_summary_filtered.tsv (after step 7c) and writes <out_prefix>_ete_removed.tsv (removed genes,
the ETE matched, identity, coverage) and <out_prefix>_ete_blastp.tsv (all hits at E <= 1e-5, for review).

Usage: python3 remove_literature_etes.py NODES_DIR QUERY_FAA OUT_PREFIX ETE_FAA [--blastp BLASTP] [--makeblastdb MAKEBLASTDB]
"""
import argparse
import collections
import csv
import os
import subprocess
import sys
import tempfile


def die(m):
    sys.exit(f"ERROR: {m}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("nodes_dir")
    ap.add_argument("query_faa")
    ap.add_argument("out_prefix")
    ap.add_argument("ete_faa")
    ap.add_argument("--blastp", default="blastp")
    ap.add_argument("--makeblastdb", default="makeblastdb")
    ap.add_argument("--min-ident", type=float, default=50.0)
    ap.add_argument("--min-cov", type=float, default=50.0)
    a = ap.parse_args()
    summ = os.path.join(a.nodes_dir, f"{a.out_prefix}_nodes_summary_filtered.tsv")
    for p in (summ, a.query_faa, a.ete_faa):
        if not os.path.exists(p):
            die(f"not found: {p}")
    sym = {}
    for line in open(a.ete_faa):
        if line.startswith(">"):
            f = line[1:].split()
            if len(f) < 2:
                die(f"ETE header without symbol: {line.strip()}")
            sym[f[0]] = f[1]
    raw = os.path.join(a.nodes_dir, f"{a.out_prefix}_ete_blastp.tsv")
    with tempfile.TemporaryDirectory(dir=a.nodes_dir) as tmp:
        db = os.path.join(tmp, "ete")
        subprocess.run([a.makeblastdb, "-in", a.ete_faa, "-dbtype", "prot", "-out", db], check=True, stdout=subprocess.DEVNULL)
        subprocess.run([a.blastp, "-query", a.query_faa, "-db", db, "-seg", "no", "-evalue", "1e-5", "-max_target_seqs", "1000",
                        "-num_threads", "1", "-out", raw,
                        "-outfmt", "6 qseqid sseqid pident length nident qstart qend sstart send evalue bitscore qlen slen"], check=True)
    pairs = collections.defaultdict(list)
    for line in open(raw):
        f = line.rstrip("\n").split("\t")
        pairs[(f[0].split("__")[0], f[1])].append(f)
    best = {}
    for (q, s), hs in pairs.items():
        qlen, slen = int(hs[0][11]), int(hs[0][12])
        qcov = 100.0 * len({i for h in hs for i in range(int(h[5]), int(h[6]) + 1)}) / qlen
        scov = 100.0 * len({i for h in hs for i in range(int(h[7]), int(h[8]) + 1)}) / slen
        cov = max(qcov, scov)
        ident = 100.0 * sum(int(h[4]) for h in hs) / sum(int(h[3]) for h in hs)
        ok = ident >= a.min_ident and cov >= a.min_cov
        rec = (ok, ident * cov, s, ident, qcov, min(float(h[9]) for h in hs), scov)
        if q not in best or rec[:2] > best[q][:2]:
            best[q] = rec
    with open(summ, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        fields, rows = rd.fieldnames, list(rd)
    kept = [r for r in rows if not (r["geneID"] in best and best[r["geneID"]][0])]
    gone = [r for r in rows if r["geneID"] in best and best[r["geneID"]][0]]
    tmp_out = summ + ".tmp"
    with open(tmp_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(kept)
    os.replace(tmp_out, summ)
    with open(os.path.join(a.nodes_dir, f"{a.out_prefix}_ete_removed.tsv"), "w") as fh:
        fh.write("geneID\tGene_type\tORF_Node\tete_gene\tete_protein\tidentity_pct\tcandidate_cov_pct\tete_cov_pct\tevalue\n")
        for r in gone:
            b = best[r["geneID"]]
            fh.write(f"{r['geneID']}\t{r['Gene_type']}\t{r['ORF_Node']}\t{sym.get(b[2], '?')}\t{b[2]}\t{b[3]:.1f}\t{b[4]:.1f}\t{b[6]:.1f}\t{b[5]:.2g}\n")
    assert len(kept) + len(gone) == len(rows)
    n_hit = sum(1 for r in rows if r["geneID"] in best)
    print(f"{os.path.basename(summ)}: in={len(rows)}, kept={len(kept)}, removed={len(gone)} "
          f"{dict(collections.Counter(r['Gene_type'] for r in gone))} (literature ETEs, >= {a.min_ident:g}% identity over "
          f">= {a.min_cov:g}% of the candidate or the ETE); {n_hit} with any ETE hit at E <= 1e-5", file=sys.stderr)


if __name__ == "__main__":
    main()
