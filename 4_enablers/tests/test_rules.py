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
  4. focal_frame_stops (2026-11-15): a 33-nt focal CDS and an ancestor lacking its first two bases. Stage 3 step 13
     reads the ancestor from focal position 3 and reports only a TGA at 12 (off the focal frame); the focal-frame stops
     are TAA at 7, TAG at 16 and the shared focal stop TAA at 31. The rebuilt table holds 7, 16 and 31, keeps the
     metadata of the step 13 row, and keeps the step 13 table as .stage3.
  5. filter_stop_events with a focal CDS that lacks its stop codon (as in the Brassicaceae Arabis alpina alignments):
     a stop at L-2 is a premature stop and is kept.
  6. focal_frame_stops edge cases: a stop triplet split by a row insertion (a base in a focal-gap column) is still the
     triplet of the three focal columns; a codon with a gap is no stop; a locus without step 13 rows takes Chromosome,
     Start and End from the file name (chromosome "scaffold_7" with an underscore) and strand NA; a ".p" AncSeq ID
     finds its alignment.
  7. focal_frame_stops against an independent codon walk on 300 random alignments (gaps, insertions, stop-rich).
  8. Step 5a then 5b, twice: the same final table both times, the start-codon row added once, and the step 5b backup
     made from the focal-frame table.
  9. Tree configs: data/tree_config.json excludes N0, N1 as ORF nodes and N0, N1 as enabler nodes, and equals
     make_tree_config.py --ingroup-crown N3 --allow-orf-node N2 --allow-orf-node N3 in both lists; the published preset
     equals --ingroup-crown N3 (N0-N3). Step 7c with the default nodes keeps a locus with the ORF in N3 and removes one
     with the ORF in N1.
 10. drop_unaligned_end_indels: in a 96-nt CDS, the table row of N6 at 8-20 (a deletion trimmed by the stage 3 mask)
     sits in a gap run 1-20 of its alignment row, an unaligned 5' end, and is removed; N4's 2-nt gap at 30-31 is
     internal and kept; Arabis alpina's row at 90-93 is part of a run reaching base 96 and is removed; N9's row at 50-51
     is not a gap in the alignment and is kept; N3's row at 44-50 lies in an internal gap run 40-50 (the remainder of
     a gap trimmed by the stage 3 mask) and is removed only with indel_mask_remainder_from_alignment.
 11. older_orf_pdng: DNG gA (ORF node N9; ORF also at N3, absent at N4 and N6) becomes PDNG; DNG gB (ORF only at N9)
     stays DNG; DNG gC (species-specific; ORF at N6) becomes PDNG; the PDNG gD is untouched; both summary files change.

Usage: python 4_enablers/tests/test_rules.py   (or pytest); exits non-zero on failure
"""
import csv
import json
import os
import random
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


def test_focal_frame_stops():
    focal = "ATGGCTCAAAGTGGAAAGCTTCATGGTGCATAA"
    anc = "--GGCTTAAAGTGAATAGCTTCATGGTGCATAA"
    with tempfile.TemporaryDirectory() as d:
        w(f"{d}/aln/{F}_g1_chr1_1_33.aln", f">{F}.chr1\n{focal}\n>N4.x\n{anc}\n")
        fields = ["file_name", "gff_transcriptID", "Chromosome", "Start", "End", "Species", "StopCodonStart_RefPos",
                  "StopCodonEnd_RefPos", "StopCodon", "gff_strand_ref", "cds_length", "reference_length"]
        old = [f"{F}_g1_chr1_1_33", "g1", "chr1", "1", "33", "N4", "12", "14", "TGA", "+", "31", "2"]
        w(f"{d}/stops/{F}_stops.tsv", "\t".join(fields) + "\n" + "\t".join(old) + "\n")
        w(f"{d}/anc.tsv", "geneID\tNode_AncSeq_ORF\tN4\ng1\tN9\t0\n")
        run("focal_frame_stops.py", f"{d}/stops/{F}_stops.tsv", f"{d}/aln", f"{d}/anc.tsv", "--focal-species", F, opts=f"{d}/none.json")
        rows = list(csv.DictReader(open(f"{d}/stops/{F}_stops.tsv"), delimiter="\t"))
        assert [(r["StopCodonStart_RefPos"], r["StopCodon"]) for r in rows] == [("7", "TAA"), ("16", "TAG"), ("31", "TAA")], rows
        assert all(r["gff_strand_ref"] == "+" and r["reference_length"] == "33" and r["cds_length"] == "31" for r in rows)
        assert open(f"{d}/stops/{F}_stops.tsv.stage3").read().count("TGA") == 1
        run("focal_frame_stops.py", f"{d}/stops/{F}_stops.tsv", f"{d}/aln", f"{d}/anc.tsv", "--focal-species", F, opts=f"{d}/none.json")
        assert list(csv.DictReader(open(f"{d}/stops/{F}_stops.tsv"), delimiter="\t")) == rows  # rerun gives the same table


def test_terminal_stop_needs_focal_stop():
    cds = CDS[:-3]  # 93 nt, no stop codon
    with tempfile.TemporaryDirectory() as d:
        opts = f"{d}/opts.json"
        w(opts, json.dumps({"drop_terminal_stop": True, "stops_focal_frame_only": True}))
        w(f"{d}/aln/{F}_g1_chr1_1_93.aln", f">{F}.chr1\n{cds}\n>N9.x\n{cds}\n")
        w(f"{d}/stops/all_stops_nodes_{F}.tsv", "geneID\tStopCodonStart_RefPos\tNodes_stops\ng1\t91\tN6|N4\n")
        run("filter_stop_events.py", f"{d}/stops", f"{d}/aln", "--focal-species", F, opts=opts)
        kept = [r["StopCodonStart_RefPos"] for r in csv.DictReader(open(f"{d}/stops/all_stops_nodes_{F}.tsv"), delimiter="\t")]
        assert kept == ["91"], kept


FIELDS = ["file_name", "gff_transcriptID", "Chromosome", "Start", "End", "Species", "StopCodonStart_RefPos",
          "StopCodonEnd_RefPos", "StopCodon", "gff_strand_ref", "cds_length", "reference_length"]


def test_focal_frame_stops_edge_cases():
    focal = "ATG-AAACAA" + "GCT" * 3 + "TAA"
    #        cols: ATG | - (ins col) | AAA | CAA ...; row: T in insertion column between T and AA -> focal cols read TAA
    row = "TAGCTAA-AA" + "GCT" * 3 + "TAA"
    with tempfile.TemporaryDirectory() as d:
        w(f"{d}/aln/{F}_g2_scaffold_7_100_120.aln", f">{F}.x\n{focal}\n>N6.x\n{row}\n")
        w(f"{d}/stops/{F}_stops.tsv", "\t".join(FIELDS) + "\n")
        w(f"{d}/anc.tsv", "geneID\tNode_AncSeq_ORF\tN6\ng2.p\tN9\t0\n")
        run("focal_frame_stops.py", f"{d}/stops/{F}_stops.tsv", f"{d}/aln", f"{d}/anc.tsv", "--focal-species", F, opts=f"{d}/none.json")
        rows = list(csv.DictReader(open(f"{d}/stops/{F}_stops.tsv"), delimiter="\t"))
        got = [(r["StopCodonStart_RefPos"], r["StopCodon"]) for r in rows]
        # codon 1 = T,A,G -> TAG (the insertion column C is skipped); codon 2 = T,A,A at cols 4,5,6 -> TAA;
        # codon 3 = cols 8,9 with a gap at col 7 -> '-AA', no stop; last codon TAA
        assert got == [("1", "TAG"), ("4", "TAA"), ("19", "TAA")], got
        r0 = rows[0]
        assert (r0["gff_transcriptID"], r0["Chromosome"], r0["Start"], r0["End"], r0["gff_strand_ref"]) == \
            ("g2", "scaffold_7", "100", "120", "NA"), r0


def _walk(ref, row):
    out, pos, cod, start = [], 0, [], None
    for c, fch in enumerate(ref):
        if fch == "-":
            continue
        pos += 1
        if (pos - 1) % 3 == 0:
            cod, start = [], pos
        cod.append(row[c])
        if len(cod) == 3 and "".join(cod) in ("TAA", "TAG", "TGA"):
            out.append((start, "".join(cod)))
    return out


def test_focal_frame_stops_random():
    sys.path.insert(0, S)
    from focal_frame_stops import focal_frame_stops
    rng = random.Random(20261115)
    n_stops = 0
    for _ in range(300):
        n = rng.randrange(9, 120)
        ref, row = [], []
        for _ in range(n):
            if rng.random() < 0.1:
                ref.append("-")
                row.append(rng.choice("ACGT"))
            ref.append(rng.choice("ACGT"))
            row.append("-" if rng.random() < 0.15 else rng.choice("AATTGC"))
        ref, row = "".join(ref), "".join(row)
        exp = _walk(ref, row)
        got = focal_frame_stops(ref, row)
        assert got == exp, (ref, row, got, exp)
        n_stops += len(exp)
    assert n_stops > 100, n_stops  # the generator must actually produce stops


def test_step5a_then_5b_twice():
    with tempfile.TemporaryDirectory() as d:
        opts = f"{d}/opts.json"
        w(opts, json.dumps({"start_codon_events": True, "start_codon_extension_max": 50, "orf_count_stop": True}))
        n6 = "ATA" + "GCT" * 2 + "TAA" + "GCT" * 27 + "TAA"
        w(f"{d}/aln/{F}_g1_chr1_1_96.aln", f">{F}.chr1\n{CDS}\n>N6.x\n{n6}\n")
        old = [f"{F}_g1_chr1_1_96", "g1", "chr1", "1", "96", "N6", "11", "13", "TAA", "+", "96", "3"]
        w(f"{d}/stops/{F}_stops.tsv", "\t".join(FIELDS) + "\n" + "\t".join(old) + "\n")
        w(f"{d}/anc.tsv", "geneID\tNode_AncSeq_ORF\tN6\ng1\tN9\t0\n")
        tab = f"{d}/stops/{F}_stops.tsv"
        finals = []
        for _ in range(2):
            run("focal_frame_stops.py", tab, f"{d}/aln", f"{d}/anc.tsv", "--focal-species", F, opts=opts)
            assert not os.path.exists(tab + ".pre_start")
            run("add_start_codon_events.py", f"{d}/stops", f"{d}/aln", f"{d}/anc.tsv", "--focal-species", F,
                "--table", tab, opts=opts)
            finals.append(open(tab).read())
            pre = [r["StopCodonStart_RefPos"] for r in csv.DictReader(open(tab + ".pre_start"), delimiter="\t")]
            assert pre == ["10", "94"], pre  # step 5b backup = the focal-frame table (TAA at codon 4 and the focal stop)
        assert finals[0] == finals[1]
        rows = [(r["StopCodonStart_RefPos"], r["StopCodon"]) for r in csv.DictReader(open(tab), delimiter="\t")]
        assert rows == [("10", "TAA"), ("94", "TAA"), ("1", "noATG_ATA")], rows
        assert "11\t13" in open(tab + ".stage3").read()


def test_tree_configs_and_step7c():
    data = os.path.join(os.path.dirname(HERE), "data")
    pres = os.path.join(os.path.dirname(HERE), "presets")
    d0 = json.load(open(os.path.join(data, "tree_config.json")))
    p0 = json.load(open(os.path.join(pres, "tree_config_brassicaceae_published.json")))
    assert d0["excluded_orf_nodes"] == ["N0", "N1"] and d0["excluded_enabler_nodes"] == ["N0", "N1"]
    assert p0["excluded_orf_nodes"] == ["N0", "N1", "N2", "N3"] and p0["excluded_enabler_nodes"] == ["N0", "N1"]
    assert {k: v for k, v in d0.items() if k != "excluded_orf_nodes"} == {k: v for k, v in p0.items() if k != "excluded_orf_nodes"}
    with tempfile.TemporaryDirectory() as d:
        for extra, ref in (([], p0), (["--allow-orf-node", "N2", "--allow-orf-node", "N3"], d0)):
            subprocess.run([sys.executable, os.path.join(S, "make_tree_config.py"), os.path.join(data, "species_tree_N.nwk"),
                            "--ingroup-crown", "N3", "--outdir", d] + extra, check=True, capture_output=True)
            m = json.load(open(f"{d}/tree_config.json"))
            assert m["excluded_orf_nodes"] == ref["excluded_orf_nodes"], (extra, m["excluded_orf_nodes"])
            assert m["excluded_enabler_nodes"] == ref["excluded_enabler_nodes"]
            assert m["node_clades"] == ref["node_clades"]
        nodes = [f"N{i}" for i in range(22)]
        w(f"{d}/nodes/sp_nodes_summary_filtered.tsv", "geneID\tGene_type\tORF_Node\ngA\tPDNG\tN3\ngB\tDNG\tN9\n")
        cells = {"gA": {"N3", "N4"}, "gB": {"N1", "N9"}}
        w(f"{d}/anc.tsv", "geneID\tNode_AncSeq_ORF\t" + "\t".join(nodes) + "\n" + "".join(
            g + "\tNA\t" + "\t".join("1" if n in c else "0" for n in nodes) + "\n" for g, c in cells.items()))
        env = dict(os.environ, LINGUA_CONFIG=os.path.join(data, "tree_config.json"))
        subprocess.run([sys.executable, os.path.join(S, "remove_conserved_orfs.py"), f"{d}/nodes", f"{d}/anc.tsv", "sp"]
                       + d0["excluded_orf_nodes"], check=True, capture_output=True, env=env)
        kept = [r["geneID"] for r in csv.DictReader(open(f"{d}/nodes/sp_nodes_summary_filtered.tsv"), delimiter="\t")]
        assert kept == ["gA"], kept


def test_drop_unaligned_end_indels():
    n6 = "-" * 20 + CDS[20:]
    n4 = CDS[:29] + "--" + CDS[31:]
    aal = CDS[:89] + "-" * 7
    n3 = CDS[:39] + "-" * 11 + CDS[50:]
    rows = [("N6", 8, 20, 13), ("N4", 30, 31, 2), ("Arabis_alpina", 90, 93, 4), ("N9", 50, 51, 2), ("N3", 44, 50, 7)]
    hdr = ["file_name", "gff_transcriptID", "Chromosome", "Start", "End", "Species", "IndelStart_RefPos",
           "IndelEnd_RefPos", "IndelLength", "gff_strand_ref", "cds_length", "reference_length"]
    with tempfile.TemporaryDirectory() as d:
        fn = f"{F}_g1_chr1_1_96"
        w(f"{d}/aln/{fn}.aln",
          f">{F}.c\n{CDS}\n>N6.x\n{n6}\n>N4.x\n{n4}\n>Arabis_alpina.x\n{aal}\n>N9.x\n{CDS}\n>N3.x\n{n3}\n")
        w(f"{d}/ind/{F}_indels_summary.tsv", "\t".join(hdr) + "\n" + "".join(
            "\t".join([fn, "g1", "chr1", "1", "96", sp, str(a), str(b), str(n), "+", "96", "96"]) + "\n" for sp, a, b, n in rows))
        tab = f"{d}/ind/{F}_indels_summary.tsv"
        w(f"{d}/ends.json", json.dumps({"indel_unaligned_end_from_alignment": True}))
        w(f"{d}/both.json", json.dumps({"indel_unaligned_end_from_alignment": True, "indel_mask_remainder_from_alignment": True}))
        # reruns read the .pre_unaligned original, so each run depends only on its own options
        for opts, exp_kept, exp_gone in ((f"{d}/ends.json", ["N4", "N9", "N3"], ["N6", "Arabis_alpina"]),
                                         (f"{d}/both.json", ["N4", "N9"], ["N6", "Arabis_alpina", "N3"]),
                                         (f"{d}/both.json", ["N4", "N9"], ["N6", "Arabis_alpina", "N3"])):
            run("drop_unaligned_end_indels.py", tab, f"{d}/aln", "--focal-species", F, opts=opts)
            kept = [r["Species"] for r in csv.DictReader(open(tab), delimiter="\t")]
            gone = [(r["Species"], r["reason"]) for r in csv.DictReader(open(tab + ".unaligned_end_rows"), delimiter="\t")]
            assert kept == exp_kept and [g for g, _ in gone] == exp_gone, (opts, kept, gone)
        assert dict(gone) == {"N6": "unaligned_end", "Arabis_alpina": "unaligned_end", "N3": "mask_remainder"}, gone
        assert len(list(csv.DictReader(open(tab + ".pre_unaligned"), delimiter="\t"))) == 5


def test_older_orf_pdng():
    nodes = [f"N{i}" for i in range(22)]
    cells = {"gA": {"N3", "N9"}, "gB": {"N9"}, "gC": {"N6"}, "gD": {"N3", "N9"}}
    calls = [("gA", "N9", "DNG"), ("gB", "N9", "DNG"), ("gC", f"{F}-specific", "DNG"), ("gD", "N9", "PDNG")]
    with tempfile.TemporaryDirectory() as d:
        nd = f"{d}/final/nodes_summary_by_species"
        body = "geneID\tORF_Node\tEnablers_Nodes\tGene_type\n" + "".join(f"{g}\t{o}\tN4|N6\t{t}\n" for g, o, t in calls)
        w(f"{nd}/athaliana_nodes_summary_filtered.tsv", body)
        w(f"{nd}/athaliana_nodes_summary.tsv", body)
        w(f"{d}/anc.tsv", "geneID\tNode_AncSeq_ORF\t" + "\t".join(nodes) + "\n" + "".join(
            g + "\tNA\t" + "\t".join("1" if n in c else "0" for n in nodes) + "\n" for g, c in cells.items()))
        w(f"{d}/inputs.tsv", "\t".join([F, "athaliana", "athaliana", "athaliana", "aln", "ind", "stp", f"{d}/anc.tsv"]) + "\n")
        run("older_orf_pdng.py", f"{d}/final", f"{d}/inputs.tsv", opts=f"{d}/none.json")
        for suf in ("_nodes_summary_filtered.tsv", "_nodes_summary.tsv"):
            got = {r["geneID"]: r["Gene_type"] for r in csv.DictReader(open(f"{nd}/athaliana{suf}"), delimiter="\t")}
            assert got == {"gA": "PDNG", "gB": "DNG", "gC": "PDNG", "gD": "PDNG"}, got
        log = [r["geneID"] + ":" + r["older_ancestors_with_ORF"] for r in csv.DictReader(open(f"{nd}/athaliana_older_orf_pdng.log"), delimiter="\t")]
        assert log == ["gA:N3", "gC:N6"], log


if __name__ == "__main__":
    for t in (test_orf_call_frame, test_filter_stop_events, test_start_codon_extension, test_focal_frame_stops,
              test_terminal_stop_needs_focal_stop, test_focal_frame_stops_edge_cases, test_focal_frame_stops_random,
              test_step5a_then_5b_twice, test_tree_configs_and_step7c, test_drop_unaligned_end_indels,
              test_older_orf_pdng):
        t()
        print(f"{t.__name__} PASSED")
