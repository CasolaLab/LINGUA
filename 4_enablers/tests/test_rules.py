"""
test_rules.py

Toy checks of the stage 4 rules added on 2026-11-14, each with the answer worked out by hand.
  1. orf_call: an ancestor with a 1-nt deletion at CDS position 41 of a 96-nt CDS whose three frames are all open. Its
     ORF starts at the focal ATG and shifts frame at the gap, so 14 of its codons (42 of 96 nt, 43.8%) start on a focal codon boundary, while
     the longest unbroken focal-frame stretch of the row is the 55 nt after the gap (57.3%). With orf_min_in_frame 50,
     measure "row" (the 2026-11-13 runs) calls an ORF and measure "orf" does not.
  2. filter_stop_events: in a 96-nt CDS, stops at 1 (start-codon event), 10 (focal frame), 5 (outside the focal frame)
     and 94 (the focal stop codon). Kept 1 and 10; removed 5 (outside_focal_frame) and 94 (focal_terminal_stop).
  3. add_start_codon_events with start_codon_extension_max 50: N4 has ATA at the start and an ATG at codon 3 opening an
     in-frame ORF to the focal stop (reads 90 of 96 nt in frame), so it is not a carrier; N6 has ATA and a TAA at codon 4
     with no later in-frame ATG, so it is a carrier.

Usage: python 4_enablers/tests/test_rules.py   (or pytest); exits non-zero on failure
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
CDS = "ATG" + "GCT" * 30 + "TAA"  # 96 nt; GCT repeats have no stop in any frame


def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def run(script, *args, opts):
    env = dict(os.environ, LINGUA_CONFIG=TREE, LINGUA_OPTIONS=opts)
    p = subprocess.run([sys.executable, os.path.join(S, script)] + list(args), capture_output=True, text=True, env=env)
    assert p.returncode == 0, f"{script} failed:\n{p.stderr}"
    return p


def test_orf_call_frame():
    os.environ["LINGUA_CONFIG"] = TREE
    sys.path.insert(0, S)
    from orf_length_node import orf_call
    row = CDS[:40] + "-" + CDS[41:]
    has_row, lp, cov, inf, fr = orf_call(CDS, row, 80, 80, 50, True, "row")
    has_orf, *_ = orf_call(CDS, row, 80, 80, 50, True, "orf")
    assert round(inf, 1) == 57.3 and round(fr, 1) == 43.8, (inf, fr)
    assert has_row and not has_orf, (has_row, has_orf)
    assert orf_call(CDS, CDS, 80, 80, 50, True, "orf")[0]


def test_filter_stop_events():
    with tempfile.TemporaryDirectory() as d:
        opts = f"{d}/opts.json"
        w(opts, json.dumps({"drop_terminal_stop": True, "stops_focal_frame_only": True}))
        w(f"{d}/aln/{F}_g1_chr1_1_96.aln", f">{F}.chr1\n{CDS}\n>N9.x\n{CDS}\n")
        w(f"{d}/stops/all_stops_nodes_{F}.tsv",
          "geneID\tStopCodonStart_RefPos\tNodes_stops\n" + "".join(f"g1\t{p}\tN6|N4\n" for p in (1, 5, 10, 94)))
        run("filter_stop_events.py", f"{d}/stops", f"{d}/aln", "--focal-species", F, opts=opts)
        kept = [r["StopCodonStart_RefPos"] for r in csv.DictReader(open(f"{d}/stops/all_stops_nodes_{F}.tsv"), delimiter="\t")]
        rem = {r["StopCodonStart_RefPos"]: r["reason"]
               for r in csv.DictReader(open(f"{d}/stops/removed_stop_events_{F}.tsv"), delimiter="\t")}
        assert kept == ["1", "10"], kept
        assert rem == {"5": "outside_focal_frame", "94": "focal_terminal_stop"}, rem


def test_start_codon_extension():
    with tempfile.TemporaryDirectory() as d:
        opts = f"{d}/opts.json"
        w(opts, json.dumps({"start_codon_events": True, "start_codon_extension_max": 50, "orf_count_stop": True}))
        n4 = "ATA" + "GCT" + "ATG" + CDS[9:]
        n6 = "ATA" + "GCT" * 2 + "TAA" + "GCT" * 27 + "TAA"
        w(f"{d}/aln/{F}_g1_chr1_1_96.aln", f">{F}.chr1\n{CDS}\n>N4.x\n{n4}\n>N6.x\n{n6}\n")
        fields = ["file_name", "gff_transcriptID", "Chromosome", "Start", "End", "Species", "StopCodonStart_RefPos",
                  "StopCodonEnd_RefPos", "StopCodon", "gff_strand_ref", "cds_length", "reference_length"]
        w(f"{d}/stops/{F}_stops.tsv", "\t".join(fields) + "\n")
        w(f"{d}/anc.tsv", "geneID\tNode_AncSeq_ORF\tN4\tN6\ng1\tN9\t0\t0\n")
        run("add_start_codon_events.py", f"{d}/stops", f"{d}/aln", f"{d}/anc.tsv", "--focal-species", F,
            "--table", f"{d}/stops/{F}_stops.tsv", opts=opts)
        added = {r["Species"] for r in csv.DictReader(open(f"{d}/stops/{F}_stops.tsv"), delimiter="\t")}
        ext = {r["taxon"]: r["orf_in_frame_pct"]
               for r in csv.DictReader(open(f"{d}/stops/start_codon_extension_rows_{F}.tsv"), delimiter="\t")}
        assert added == {"N6"}, added
        assert ext == {"N4": "93.8"}, ext


if __name__ == "__main__":
    for t in (test_orf_call_frame, test_filter_stop_events, test_start_codon_extension):
        t()
        print(f"{t.__name__} PASSED")
