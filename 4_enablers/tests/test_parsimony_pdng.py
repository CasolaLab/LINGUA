"""
test_parsimony_pdng.py

Toy check of step 7i (parsimony_pdng.py, method option orf_clade_parsimony_pdng) on the Brassicaceae tree of
data/tree_config.json, answers worked out by hand before running. Focal Arabidopsis_thaliana; every DNG has ORF node N6
(A. thaliana, A. halleri, A. lyrata, Boechera, Camelina x2, Capsella x2) and one chosen stop event at CDS 10 (TAA in
carriers, GCT elsewhere), Stops_Node N2|N3. "Outside" = the N5 clade (11 species), Aethionema, Tarenaya and Theobroma.
N4 = N6 + Cardamine; N6 = N9 (Arabidopsis) + N10 (Boechera, Camelina, Capsella).
  gA  carriers: outside + Cardamine; N6 intact. N6 = 0 costs 1 (change on the N6 branch); N6 = 1 costs 2 (N9 and N10
      branches) -> {0} -> stays DNG.
  gB  carriers: outside + N10; Cardamine and N9 intact. N4 = N6 = 0 costs 2 (N4 branch, N10 branch); N4 = N6 = 1 costs 2
      (Cardamine, N9 branches); mixed assignments cost 3 -> N6 {0,1} -> PDNG (a sister lineage shares the ancestral state).
  gC  carriers: outside; Cardamine and N10 missing (no rows); N9 intact. N6 = 0 costs 1 (N4 or N6 branch); N6 = 1 costs 1
      (N9 branch) -> {0,1} -> PDNG (missing data only: with Cardamine and N10 intact, N6 = 1 costs 3, N6 = 0 costs 1).
  gD  carriers: outside + Cardamine + N10; N9 intact. N6 = 1 costs 1 (N9 branch), N6 = 0 costs 2 -> {1} -> PDNG.
  gE  ORF node Arabidopsis_thaliana-specific, carriers as gD: the focal tip is fixed at 0 -> stays DNG.
  gF  a PDNG row: untouched.
  gG.p  as gB, with the Phytozome ".p" suffix in the summary table only (alignment and detail table without it) -> PDNG.
Usage: python 4_enablers/tests/test_parsimony_pdng.py   (or pytest)
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(os.path.dirname(HERE), "scripts")
TREE = os.path.join(os.path.dirname(HERE), "data", "tree_config.json")
F = "Arabidopsis_thaliana"
INTACT = "ATG" + "GCT" * 8 + "TAA"   # 30 nt; codon 4 (positions 10-12) is GCT
STOP = INTACT[:9] + "TAA" + INTACT[12:]
CL = {k: set(v) for k, v in json.load(open(TREE))["node_clades"].items()}
TAXA = sorted(CL["N0"] - {F})
OUT = (CL["N5"] | {"Aethionema_arabicum", "Tarenaya_hassleriana", "Theobroma_cacao"})
CARD = {"Cardamine_amara_amara"}
N10 = CL["N10"]
CASES = {  # gene: (carriers, missing, ORF node)
    "gA": (OUT | CARD, set(), "N6"),
    "gB": (OUT | N10, set(), "N6"),
    "gC": (OUT, CARD | N10, "N6"),
    "gD": (OUT | CARD | N10, set(), "N6"),
    "gE": (OUT | CARD | N10, set(), f"{F}-specific"),
    "gG.p": (OUT | N10, set(), "N6"),  # as gB; Phytozome ".p" in the summary only, not in alignment or detail table
}


def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def test_parsimony_pdng():
    assert CL["N4"] == CL["N6"] | CARD and CL["N6"] == CL["N9"] | N10 and F in CL["N9"]
    with tempfile.TemporaryDirectory() as d:
        hdr_s = ["file_name", "gff_transcriptID", "Chromosome", "Start", "End", "Species", "StopCodonStart_RefPos",
                 "StopCodonEnd_RefPos", "StopCodon", "gff_strand_ref", "cds_length", "reference_length"]
        hdr_i = ["file_name", "gff_transcriptID", "Chromosome", "Start", "End", "Species", "IndelStart_RefPos",
                 "IndelEnd_RefPos", "IndelLength", "gff_strand_ref", "cds_length", "reference_length"]
        stops = []
        for g, (car, miss, _) in list(CASES.items()) + [("gF", (OUT, set(), "N6"))]:
            fa = f">{F}.c\n{INTACT}\n>N2.x\n{STOP}\n>N3.x\n{STOP}\n"
            for t in TAXA:
                if t not in miss:
                    fa += f">{t}.x\n{STOP if t in car else INTACT}\n"
            gk = g[:-2] if g.endswith(".p") else g
            w(f"{d}/aln/{F}_{gk}_chr1_1_30.aln", fa)
            for t in sorted(car) + ["N2", "N3"]:
                stops.append([f"{F}_{gk}_chr1_1_30", gk, "chr1", "1", "30", t, "10", "12", "TAA", "+", "30", "30"])
        w(f"{d}/stp.tsv", "\t".join(hdr_s) + "\n" + "".join("\t".join(x) + "\n" for x in stops))
        w(f"{d}/ind.tsv", "\t".join(hdr_i) + "\n")
        body = "geneID\tGene_type\tORF_Node\tEnablers_Nodes\tIndels_Node\tStops_Node\tIndels_event_RefPos\tStops_event_RefPos\n"
        for g, (_, _, orf) in CASES.items():
            body += f"{g}\tDNG\t{orf}\tN3\tNA\tN2|N3\tNA\t10\n"
        body += "gF\tPDNG\tN6\tN3\tNA\tN2|N3\tNA\t10\n"
        nd = f"{d}/final/nodes_summary_by_species"
        w(f"{nd}/athaliana_nodes_summary_filtered.tsv", body)
        w(f"{nd}/athaliana_nodes_summary.tsv", body)
        w(f"{d}/inputs.tsv", "\t".join([F, "Athaliana", "athaliana", "athaliana", f"{d}/aln", f"{d}/ind.tsv",
                                        f"{d}/stp.tsv", "anc"]) + "\n")
        w(f"{d}/opts.json", '{"orf_clade_parsimony_pdng": true}')
        p = subprocess.run([sys.executable, os.path.join(S, "parsimony_pdng.py"), f"{d}/final", f"{d}/inputs.tsv"],
                           capture_output=True, text=True,
                           env=dict(os.environ, LINGUA_CONFIG=TREE, LINGUA_OPTIONS=f"{d}/opts.json"))
        assert p.returncode == 0, p.stderr
        want = {"gA": "DNG", "gB": "PDNG", "gC": "PDNG", "gD": "PDNG", "gE": "DNG", "gF": "PDNG", "gG.p": "PDNG"}
        for suf in ("_nodes_summary_filtered.tsv", "_nodes_summary.tsv"):
            got = {r["geneID"]: r["Gene_type"] for r in csv.DictReader(open(f"{nd}/athaliana{suf}"), delimiter="\t")}
            assert got == want, got
        log = {r["geneID"]: r for r in csv.DictReader(open(f"{nd}/athaliana_parsimony_pdng.log"), delimiter="\t")}
        assert {g: r["orf_node_states"] for g, r in log.items()} == {"gA": "0", "gB": "01", "gC": "01", "gD": "1", "gE": "0", "gG.p": "01"}, log
        assert set(log["gC"]["living_missing"].split(",")) == CARD | N10


if __name__ == "__main__":
    test_parsimony_pdng()
    print("test_parsimony_pdng PASSED")
