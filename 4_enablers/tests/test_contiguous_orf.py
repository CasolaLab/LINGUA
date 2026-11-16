"""
test_contiguous_orf.py

Toy checks of orf_frame_measure "contiguous" (T2T, 2026-11-15), each answer worked out by hand before running. The focal CDS is
ATG + GCT x 30 + TAA (96 nt; GCT repeats have no stop in any frame). Thresholds: span >= 80 %, row bases >= 80 %.
  1. identical row: span 100, bases 100, ORF.
  2. TAA at focal codon 11 (positions 31-33): stretch 1-30 = 31.25 %, no later ATG: no ORF.
  3. TAA at codon 3 (7-9), ATG at codon 5 (13-15), intact to the end: stretch 13-96 = 87.5 %: ORF (a start after the stop counts).
  4. 1-nt deletion at 41 and 2-nt deletion at 61-62 (frame restored after 62): stretch 1-40 = 41.67 %, no ORF; the "orf" measure
     with a 50 % cutoff calls an ORF on the same row (pieces added up), which is the case the new measure removes.
  5. codon-aligned 3-nt deletion at 31-33: not a break; span 100, bases 93/96 = 96.875 %: ORF.
  6. codon-aligned 30-nt deletion at 31-60: span 100 but bases 66/96 = 68.75 %: no ORF.
  7. focal start CTG (CTG + GCT x 30 + TAA), identical row: ORF from the focal start codon.
  8. 1-nt deletion at 90 (inside the last 20 %): stretch 1-89 = 92.7 %: ORF (a late frameshift leaves the ORF).
  9. add_start_codon_events with contiguous and start_codon_extension_max 80: N4 (ATA, ATG at codon 3, intact; stretch 7-96 = 93.75 %)
     is not a carrier; N6 (ATA, TAA at codon 4, no later ATG) is.
 10. orf_length_node.py end to end with the contiguous options: cells N4 = 1 (row of test 3), N6 = 0 (row of test 4), N9 = 1
     (identical); the details table carries the stretch span as orf_in_frame_pct.
Usage: python stage4/tests/test_contiguous_orf.py   (or pytest); exits non-zero on failure
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
CDS = "ATG" + "GCT" * 30 + "TAA"
os.environ["LINGUA_CONFIG"] = TREE
sys.path.insert(0, S)
from orf_length_node import focal_frame_orf, orf_call  # noqa: E402

ST = {"ATG"}


def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def call(row, ref=CDS, starts=ST):
    has, span, bases = focal_frame_orf(ref, row, starts, 80, 80)
    return has, round(span, 3), round(bases, 3)


def test_rows():
    assert call(CDS) == (True, 100.0, 100.0)
    r2 = CDS[:30] + "TAA" + CDS[33:]
    assert call(r2) == (False, 31.25, 31.25), call(r2)
    r3 = CDS[:6] + "TAA" + CDS[9:12] + "ATG" + CDS[15:]
    assert call(r3) == (True, 87.5, 87.5), call(r3)
    r4 = CDS[:40] + "-" + CDS[41:60] + "--" + CDS[62:]
    assert call(r4)[:2] == (False, 41.667), call(r4)
    assert orf_call(CDS, r4, 80, 80, 50, True, "orf")[0] is True
    assert orf_call(CDS, r4, 80, 80, 80, True, "contiguous")[0] is False
    r5 = CDS[:30] + "---" + CDS[33:]
    assert call(r5) == (True, 100.0, 96.875), call(r5)
    r6 = CDS[:30] + "-" * 30 + CDS[60:]
    assert call(r6) == (False, 100.0, 68.75), call(r6)
    c7 = "CTG" + CDS[3:]
    assert call(c7, ref=c7, starts={"ATG", "CTG"})[0] is True
    assert call(c7, ref=c7, starts={"ATG"})[0] is False
    r8 = CDS[:89] + "-" + CDS[90:]
    assert call(r8)[:2] == (True, 92.708), call(r8)


def test_start_extension_contiguous():
    with tempfile.TemporaryDirectory() as d:
        opts = f"{d}/opts.json"
        w(opts, json.dumps({"start_codon_events": True, "start_codon_extension_max": 80, "orf_count_stop": True,
                            "orf_frame_measure": "contiguous", "orf_min_len": 80}))
        n4 = "ATA" + "GCT" + "ATG" + CDS[9:]
        n6 = "ATA" + "GCT" * 2 + "TAA" + "GCT" * 27 + "TAA"
        w(f"{d}/aln/{F}_g1_chr1_1_96.aln", f">{F}.chr1\n{CDS}\n>N4.x\n{n4}\n>N6.x\n{n6}\n")
        fields = ["file_name", "gff_transcriptID", "Chromosome", "Start", "End", "Species", "StopCodonStart_RefPos",
                  "StopCodonEnd_RefPos", "StopCodon", "gff_strand_ref", "cds_length", "reference_length"]
        w(f"{d}/stops/{F}_stops.tsv", "\t".join(fields) + "\n")
        w(f"{d}/anc.tsv", "geneID\tNode_AncSeq_ORF\tN4\tN6\ng1\tN9\t0\t0\n")
        p = subprocess.run([sys.executable, os.path.join(S, "add_start_codon_events.py"), f"{d}/stops", f"{d}/aln", f"{d}/anc.tsv",
                            "--focal-species", F, "--table", f"{d}/stops/{F}_stops.tsv"], capture_output=True, text=True,
                           env=dict(os.environ, LINGUA_CONFIG=TREE, LINGUA_OPTIONS=opts))
        assert p.returncode == 0, p.stderr
        added = {r["Species"] for r in csv.DictReader(open(f"{d}/stops/{F}_stops.tsv"), delimiter="\t")}
        ext = {r["taxon"]: r["orf_in_frame_pct"]
               for r in csv.DictReader(open(f"{d}/stops/start_codon_extension_rows_{F}.tsv"), delimiter="\t")}
        assert added == {"N6"}, added
        assert ext == {"N4": "93.8"}, ext


def test_orf_length_node_end_to_end():
    with tempfile.TemporaryDirectory() as d:
        opts = f"{d}/opts.json"
        w(opts, json.dumps({"orf_test": "length", "orf_count_stop": True, "orf_frame_measure": "contiguous",
                            "orf_min_in_frame": 80, "orf_node_contiguous": True}))
        r3 = CDS[:6] + "TAA" + CDS[9:12] + "ATG" + CDS[15:]
        r4 = CDS[:40] + "-" + CDS[41:60] + "--" + CDS[62:]
        w(f"{d}/aln/{F}_g1_chr1_1_96.aln", f">{F}.chr1\n{CDS}\n>N4.x\n{r3}\n>N6.x\n{r4}\n>N9.x\n{CDS}\n")
        w(f"{d}/anc.tsv", "geneID\tNode_AncSeq_ORF\tN4\tN6\tN9\ng1\tN9\t0\t0\t0\n")
        p = subprocess.run([sys.executable, os.path.join(S, "orf_length_node.py"), f"{d}/anc.tsv", f"{d}/aln",
                            "--focal-species", F, "--min-len", "80", "--min-cov", "80"], capture_output=True, text=True,
                           env=dict(os.environ, LINGUA_CONFIG=TREE, LINGUA_OPTIONS=opts))
        assert p.returncode == 0, p.stderr
        row = list(csv.DictReader(open(f"{d}/anc.tsv"), delimiter="\t"))[0]
        assert (row["N4"], row["N6"], row["N9"]) == ("1", "0", "1"), row
        det = {r["node"]: r for r in csv.DictReader(open(f"{d}/anc.orf_length_details.tsv"), delimiter="\t")}
        assert det["N4"]["orf_in_frame_pct"] == "87.5" and det["N6"]["orf_in_frame_pct"] == "41.7", det


if __name__ == "__main__":
    for t in (test_rows, test_start_extension_contiguous, test_orf_length_node_end_to_end):
        t()
        print(f"{t.__name__} PASSED")
