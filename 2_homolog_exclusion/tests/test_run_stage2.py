"""
test_run_stage2.py - scripts/run_stage2.sh end to end on a hand-built example, with real BLAST+ and HMMER.

Two species, Testa alpha and Testb beta, eight candidates. Expected fate, worked out by hand:
  cand1  (Testa)  = protein F, which the database holds in Homo sapiens and Mus musculus    -> removed by BLASTp/jackhmmer
  cand3  (Testb)  = first 40 aa of F + 80 unrelated aa; 33% of it matches F, below the 50%
                    query cover, so no hit of its own; it shares HOG H1 with cand1 and aligns
                    to cand1's hit region over 40 aa (>= 20)                                -> removed by the HOG rule
  cand4  (Testa)  = protein G; the Testa proteome holds G2 (G with 10 changes, not a candidate) -> removed by the paralog rule
  cand7  (Testb)  = a TE protein                                                            -> removed by the TE screen
  cand8  (Testa)  = encoded by the organelle genome                                         -> removed by the organelle screen
  cand5, cand6    unrelated                                                                 -> pass
CDD is left out with SKIP_SCREENS=cdd (and the status says so); InterProScan is a stand-in that finds nothing.
Also checked: SKIP_SCREENS cannot name the HOG or paralog rule; a missing HOGS stops the run before it starts.
Needs blastp, tblastn, makeblastdb and jackhmmer on PATH (prints SKIP otherwise).
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DRIVER = os.path.join(HERE, "..", "scripts", "run_stage2.sh")
AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA",
         "L": "CTT", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT",
         "W": "TGG", "Y": "TAT"}


def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def fasta(recs):
    return "".join(f">{i}\n{s}\n" for i, s in recs)


def main():
    missing = [t for t in ("blastp", "tblastn", "makeblastdb", "jackhmmer") if not shutil.which(t)]
    if missing:
        print(f"SKIP: {missing} not on PATH")
        return
    rnd = random.Random(7)
    seq = lambda n: "M" + "".join(rnd.choice(AA) for _ in range(n - 1))
    F, G, T, O, U5, U6, TAIL = seq(120), seq(120), seq(120), seq(100), seq(120), seq(120), seq(81)[1:]
    G2 = list(G)
    for i in range(5, 120, 12):
        G2[i] = "W" if G2[i] != "W" else "Y"
    G2 = "".join(G2)
    cand = {"Testa_alpha": [("cand1", F), ("cand4", G), ("cand5", U5), ("cand8", O)],
            "Testb_beta": [("cand3", F[:40] + TAIL), ("cand6", U6), ("cand7", T)]}
    with tempfile.TemporaryDirectory() as d:
        for sp, recs in cand.items():
            w(f"{d}/cand/{sp}.faa", fasta(recs))
        w(f"{d}/prot/Testa_alpha.faa", fasta(cand["Testa_alpha"] + [("G2", G2)]))
        w(f"{d}/prot/Testb_beta.faa", fasta(cand["Testb_beta"]))
        db = fasta([("hs1__Homo_sapiens", F), ("mm1__Mus_musculus", F[:60] + "W" + F[61:]), ("x1__Danio_rerio", seq(150))])
        w(f"{d}/db/animals.faa", db)
        w(f"{d}/db/te.faa", fasta([("te1__Repeat_element", T)]))
        w(f"{d}/db/org.fna", fasta([("chloro__Plastid_genome", "".join(rnd.choice("ACGT") for _ in range(300))
                                     + "".join(CODON[a] for a in O) + "TAA"
                                     + "".join(rnd.choice("ACGT") for _ in range(300)))]))
        for name, typ in (("animals.faa", "prot"), ("te.faa", "prot"), ("org.fna", "nucl")):
            subprocess.run(["makeblastdb", "-in", f"{d}/db/{name}", "-dbtype", typ, "-out", f"{d}/db/{name}"],
                           check=True, capture_output=True)
        w(f"{d}/N0.tsv", "HOG\tOG\tGene Tree Parent Clade\tTesta_alpha\tTestb_beta\n"
          "N0.HOG1\tOG1\tn0\tcand1\tcand3\nN0.HOG2\tOG2\tn0\tcand4, G2\t\n")
        stub = f"{d}/fake_interproscan.sh"
        w(stub, r'''#!/usr/bin/env bash
out=""
while [[ $# -gt 0 ]]; do
  case "$1" in -o) out="$2"; shift 2;; --version) echo "InterProScan version 5.99-test"; exit 0;; *) shift;; esac
done
: > "$out"
echo "Running the following analyses:"; echo "[Pfam-38.2]"; echo "100% done:  InterProScan analyses completed"
''')
        os.chmod(stub, 0o755)
        base = (f"CANDIDATES={d}/cand\nHOGS={d}/N0.tsv\nPROTEOMES={d}/prot\nBLAST_DBS={d}/db/animals.faa\n"
                f"JACKHMMER_DB={d}/db/animals.faa\nTE_DB={d}/db/te.faa\nTE_CONFIG={HERE}/../presets/brassicaceae_2026_te.conf\n"
                f"ORGANELLE_DB={d}/db/org.fna\nORGANELLE_CONFIG={HERE}/../presets/organelle_tblastn.conf\n"
                f"IPRSCAN_BIN={stub}\nTHREADS=1\n")

        def drive(name, extra, want_ok, want_text):
            w(f"{d}/{name}.conf", base + f"OUTDIR={d}/{name}\n" + extra)
            p = subprocess.run(["bash", DRIVER, f"{d}/{name}.conf"], capture_output=True, text=True)
            ok = p.returncode == 0
            assert ok == want_ok and want_text in p.stdout + p.stderr, f"{name}: exit {p.returncode}\n{p.stdout}\n{p.stderr}"
            return p

        drive("bad_skip", "SKIP_SCREENS=cdd,hog\n", False, "cannot be skipped")
        drive("no_hogs", f"SKIP_SCREENS=cdd\nHOGS={d}/missing.tsv\n", False, "HOGS not found")
        drive("full", "SKIP_SCREENS=cdd\n", True, "SUCCEEDED steps 1-9 (screens left out: cdd)")
        passed = sorted(open(f"{d}/full/paralog/passed_ids.tsv").read().split("\n")[1:])
        got = sorted(x.split("\t")[1] if "\t" in x else x for x in passed if x)
        assert got == ["cand5", "cand6"], (got, open(f"{d}/full/paralog/passed_ids.tsv").read())
        import csv
        gm = {r["gene_id"]: (r["verdict"], r["removed_by"])
              for r in csv.DictReader(open(f"{d}/full/combine/gene_matrix.tsv"), delimiter="\t")}
        assert gm == {"cand1": ("REMOVED", "blast;jackhmmer"), "cand3": ("PASS", "NA"), "cand4": ("PASS", "NA"),
                      "cand5": ("PASS", "NA"), "cand6": ("PASS", "NA"), "cand7": ("REMOVED", "te"),
                      "cand8": ("REMOVED", "organelle")}, gm
        hog = open(f"{d}/full/hog/hog_rule.tsv").read()
        assert "cand3" in hog and "REMOVED_HOG" in hog, hog
        par = open(f"{d}/full/paralog/paralog_rule.tsv").read()
        assert "REMOVED_PARALOG" in par and "G2" in par, par
    print("test_run_stage2: ALL PASSED")


if __name__ == "__main__":
    main()
