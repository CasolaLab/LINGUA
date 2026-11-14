"""
test_toy.py

Steps 1, 1a, 2, 4 (families) and 5 on a four-species toy tree whose answer was worked out by hand.

Tree ((A,B)N2,C)N1,D)N0; ingroup crown N1, ORFs may not originate at N0.
  a1 (A, ORF N2, DNG) and b1 (B, ORF N2, PDNG) are the same locus, found in A and B -> one DNG locus at N2, no deaths.
  c1 (C, ORF N1, PDNG), found only in C, not in any HOG -> PDNG birth at N1, death on N2 (A and B both lack it).
  a2, a3 (A-specific PDNGs, same HOG) -> one PDNG birth on A.
  a4 (A-specific DNG) and a5 (A-specific PDNG) in one HOG -> DNG birth on A; a5 dropped as a paralog of its DNG.
  a6, a7 (A-specific PDNGs, unclustered, one BLASTp hit) -> one family, a7 (longer) kept -> one PDNG birth on A.
Expected births: DNG N2 1, A 1; PDNG N1 1, A 2. Deaths: PDNG N2 1. Rates with N2 = 5: DNG birth rate 0.200.

Usage: python 5_birth_death/tests/test_toy.py   (or pytest); exits non-zero on failure
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

S = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")


def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def run(*args, fail=None):
    p = subprocess.run([sys.executable] + list(args), capture_output=True, text=True)
    if fail:
        assert p.returncode != 0 and fail in p.stderr, f"{os.path.basename(args[0])} should stop with '{fail}':\n{p.stderr}"
        return
    assert p.returncode == 0, f"{os.path.basename(args[0])} failed:\n{p.stderr}"


def build(d):
    if True:
        sp = ["A", "B", "C", "D"]
        cfg = {"source_tree": "toy", "ingroup_crown": {"label": "N1", "original": "N1"}, "species": sp,
               "ingroup_species": ["A", "B", "C"], "aliases": {s: s for s in sp},
               "node_original_names": {n: n for n in ("N0", "N1", "N2")},
               "node_order_young_to_old": ["N2", "N1", "N0"],
               "node_clades": {"N0": sp, "N1": ["A", "B", "C"], "N2": ["A", "B"]},
               "excluded_orf_nodes": ["N0"], "excluded_enabler_nodes": ["N0"]}
        w(f"{d}/tree_config.json", json.dumps(cfg))
        hdr = "file_name\tspecies_name\tchromosome\tstart\tend\n"
        w(f"{d}/coords/A_all_species_summary.tsv", hdr + "".join(
            f"A_{g}_c1_{s}_{e}.aln\t{t}\t{c}\t{s}\t{e}\n" for g, t, c, s, e in [
                ("a1", "A", "c1", 100, 200), ("a1", "B", "b1", 1000, 1100), ("a2", "A", "c1", 300, 400),
                ("a3", "A", "c1", 500, 600), ("a4", "A", "c1", 700, 800), ("a5", "A", "c1", 900, 950),
                ("a6", "A", "c2", 100, 200), ("a7", "A", "c2", 300, 400)]))
        w(f"{d}/coords/B_all_species_summary.tsv", hdr + "B_b1_b1_1000_1100.aln\tB\tb1\t1000\t1100\n"
          "B_b1_b1_1000_1100.aln\tA\tc1\t100\t200\n")
        w(f"{d}/coords/C_all_species_summary.tsv", hdr + "C_c1_x_5_50.aln\tC\tx\t5\t50\nC_c1_x_5_50.aln\tA\tc9\t5\t50\n")
        nh = "geneID\tORF_Node\tSynteny_Node\tEnablers_Nodes\tGene_type\n"
        w(f"{d}/nodes/A_nodes_summary_filtered.tsv", nh + "a1\tN2\tN1\tN1\tDNG\n" + "".join(
            f"{g}\tA-specific\tN2\tNA\t{t}\n" for g, t in
            [("a2", "PDNG"), ("a3", "PDNG"), ("a4", "DNG"), ("a5", "PDNG"), ("a6", "PDNG"), ("a7", "PDNG")]))
        w(f"{d}/nodes/B_nodes_summary_filtered.tsv", nh + "b1\tN2\tN1\tNA\tPDNG\n")
        w(f"{d}/nodes/C_nodes_summary_filtered.tsv", nh + "c1\tN1\tN0\tNA\tPDNG\n")
        th = "geneID\tNode_tblastn\tA\tB\tC\tD\n"
        w(f"{d}/tb/A_cdngs_tblastn_7_summary.tsv", th + "a1\tN2\t1\t1\t0\t0\n" + "".join(
            f"{g}\tA-specific\t1\t0\t0\t0\n" for g in ("a2", "a3", "a4", "a5", "a6", "a7")))
        w(f"{d}/tb/B_cdngs_tblastn_1_summary.tsv", th + "b1\tN2\t1\t1\t0\t0\n")
        w(f"{d}/tb/C_cdngs_tblastn_1_summary.tsv", th + "c1\tN1\t0\t0\t1\t0\n")
        w(f"{d}/N0.tsv", "HOG\tOG\tGene Tree Parent Clade\tA\tB\tC\tD\n"
          "N0.H1\tOG1\tn\ta1\tb1\t\t\nN0.H2\tOG2\tn\ta2, a3\t\t\t\nN0.H3\tOG3\tn\ta4, a5\t\t\t\n")
        w(f"{d}/prot/A.faa", ">a6__A\nMAAAA\n>a7__A\nMAAAAAAA\n")
        w(f"{d}/prot/C.faa", ">c1__C\nMCCC\n")
        w(f"{d}/len.tsv", "branch\tlength\nN1\t10\nN2\t5\nA\t2\nB\t2\nC\t7\nD\t20\n")
        return f"{d}/tree_config.json"


def test_toy():
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        run(f"{S}/step1_coordinates.py", "--tree-config", cf, "--coords", f"{d}/coords", "--nodes", f"{d}/nodes",
            "--tblastn", f"{d}/tb", "--out", f"{d}/out/coords.tsv")
        run(f"{S}/step2_loci.py", "--tree-config", cf, "--coords", f"{d}/out/coords.tsv", "--outdir", f"{d}/out")
        run(f"{S}/step4_unclustered_families.py", "select", "--tree-config", cf, "--loci", f"{d}/out/loci.tsv",
            "--hogs", f"{d}/N0.tsv", "--proteins", f"{d}/prot", "--outdir", f"{d}/out/unc")
        assert sorted(os.listdir(f"{d}/out/unc")) == ["A_unclustered.faa", "unclustered_counts.tsv"]
        w(f"{d}/out/unc/A_unclustered.blastp.tsv", "a6\ta6\t100\t5\t5\t5\t1e-10\t20\na6\ta7\t100\t5\t5\t8\t1e-6\t15\n")
        run(f"{S}/step4_unclustered_families.py", "families", "--outdir", f"{d}/out/unc")
        run(f"{S}/step5_rates.py", "--no-overlap", "--tree-config", cf, "--loci", f"{d}/out/loci.tsv", "--genes", f"{d}/out/genes.tsv",
            "--hogs", f"{d}/N0.tsv", "--branch-lengths", f"{d}/len.tsv", "--families",
            f"{d}/out/unc/unclustered_families.tsv", "--outdir", f"{d}/out")
        r = {x["Branch"]: x for x in csv.DictReader(open(f"{d}/out/rates_by_branch.tsv"), delimiter="\t")}
        assert list(r) == ["N1", "N2", "A", "B", "C", "D"], list(r)
        got = {b: tuple(int(r[b][c]) for c in ("DNG_births", "PDNG_births", "DNG_deaths", "PDNG_deaths")) for b in r}
        want = {"N1": (0, 1, 0, 0), "N2": (1, 0, 0, 1), "A": (1, 2, 0, 0), "B": (0, 0, 0, 0), "C": (0, 0, 0, 0),
                "D": (0, 0, 0, 0)}
        assert got == want, got
        assert r["N2"]["DNG_birth_rate"] == "0.200" and r["N2"]["PDNG_death_rate"] == "0.200"
        st = {x["geneID"]: x["status"] for x in csv.DictReader(open(f"{d}/out/gene_status.tsv"), delimiter="\t")}
        assert st["a5"].startswith("excluded: PDNG paralog") and st["a6"].startswith("excluded: duplicate"), st


def test_id_form_mismatch_stops():
    """HOG table written with 'db|ID|' identifiers: step 5 must stop instead of counting every gene as its own origin."""
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        w(f"{d}/N0.tsv", "HOG\tOG\tGene Tree Parent Clade\tA\tB\tC\tD\n"
          "N0.H1\tOG1\tn\tref|a1|\tref|b1|\t\t\nN0.H2\tOG2\tn\tref|a2|, ref|a3|\t\t\t\n")
        run(f"{S}/step1_coordinates.py", "--tree-config", cf, "--coords", f"{d}/coords", "--nodes", f"{d}/nodes",
            "--tblastn", f"{d}/tb", "--out", f"{d}/out/coords.tsv")
        run(f"{S}/step2_loci.py", "--tree-config", cf, "--coords", f"{d}/out/coords.tsv", "--outdir", f"{d}/out")
        run(f"{S}/step5_rates.py", "--no-overlap", "--tree-config", cf, "--loci", f"{d}/out/loci.tsv", "--genes", f"{d}/out/genes.tsv",
            "--hogs", f"{d}/N0.tsv", "--branch-lengths", f"{d}/len.tsv", "--outdir", f"{d}/out", fail="differ in form")
        w(f"{d}/nodes/A_nodes_summary_filtered.tsv", "geneID\tORF_Node\tSynteny_Node\tEnablers_Nodes\tGene_type\n"
          + "".join(f"A.{g}\tA-specific\tN2\tNA\tPDNG\n" for g in ("a1", "a2", "a3", "a4", "a5", "a6", "a7")))
        run(f"{S}/step1_coordinates.py", "--tree-config", cf, "--coords", f"{d}/coords", "--nodes", f"{d}/nodes",
            "--tblastn", f"{d}/tb", "--out", f"{d}/out2/coords.tsv", fail="probably differ in form")


def test_overlap_inside_long_gene():
    """A locus (own gene cand01_gene, minus strand) inside the intron of a 100 kb plus-strand gene g1, with small genes between
    the start of g1 and the locus: antisense_overlap (found although g1 starts far upstream). A lncRNA covering the
    locus is ignored."""
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        gff = ["##gff-version 3", "c1\tx\tgene\t1\t100000\t.\t+\t.\tID=g1;gene_biotype=protein_coding",
               "c1\tx\tmRNA\t1\t100000\t.\t+\t.\tID=g1.t;Parent=g1", "c1\tx\tCDS\t1\t300\t.\t+\t0\tID=c;Parent=g1.t",
               "c1\tx\tgene\t40000\t60000\t.\t-\t.\tID=nc1;gene_biotype=lncRNA"]
        for i, st in enumerate((5000, 9000, 20000, 30000)):
            gff.append(f"c1\tx\tgene\t{st}\t{st + 100}\t.\t+\t.\tID=s{i}")
        gff.append("c1\tx\tgene\t50001\t50100\t.\t-\t.\tID=cand01_gene")
        w(f"{d}/gff/A.gff", "\n".join(gff) + "\n")
        w(f"{d}/gffs.tsv", "species\tgff\nA\tgff/A.gff\n")
        w(f"{d}/c.tsv", "geneID\tfocal_species\tfile_name\tspecies_name\tchromosome\tstart\tend\tORF_Node\t"
          "Synteny_Node\tEnablers_Nodes\tGene_type\ncand01\tA\tf\tA\tc1\t50000\t50100\tA-specific\tN2\tNA\tPDNG\n")
        run(f"{S}/step3_overlap.py", "--tree-config", cf, "--coords", f"{d}/c.tsv", "--gff-table", f"{d}/gffs.tsv",
            "--out", f"{d}/ov.tsv")
        got = [(r["orientation"], r["other_gene"]) for r in csv.DictReader(open(f"{d}/ov.tsv"), delimiter="\t")]
        assert got == [("antisense_overlap", "g1")], got


def test_hog_lineage_split():
    """HOG H1 joins a PDNG locus at N2 (A, B) and a C-specific PDNG locus; C lies outside clade N2, so with the split
    they are two origins (N2 and C) and with merge one (N2). HOG H2 joins a PDNG locus at N1 (found in C) and an
    A-specific one; A lies inside clade N1, so one origin at N1 either way, carried by C and A, lost on B."""
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        w(f"{d}/l.tsv", "ClusterID\tNode\tGene_type\tA\tB\tC\tD\n1\tN2\tPDNG\tx1\tx1b\t\t\n2\tC-specific\tPDNG\t\t\tx2\t\n"
          "3\tN1\tPDNG\t\t\ty1\t\n4\tA-specific\tPDNG\ty2\t\t\t\n")
        w(f"{d}/g.tsv", "species\tgeneID\tORF_Node\tGene_type\tClusterID\tCluster_node\tpresent_in\n"
          "A\tx1\tN2\tPDNG\t1\tN2\tA;B\nB\tx1b\tN2\tPDNG\t1\tN2\tA;B\nC\tx2\tC-specific\tPDNG\t2\tC-specific\tC\n"
          "C\ty1\tN1\tPDNG\t3\tN1\tC\nA\ty2\tA-specific\tPDNG\t4\tA-specific\tA\n")
        w(f"{d}/h.tsv", "HOG\tOG\tGene Tree Parent Clade\tA\tB\tC\tD\nH1\tO1\tn\tx1\tx1b\tx2\t\nH2\tO2\tn\ty2\t\ty1\t\n")
        for mode, want in (("split", {"N1": (1, 0), "N2": (1, 0), "B": (0, 1), "C": (1, 0)}),
                           ("merge", {"N1": (1, 0), "N2": (1, 0), "B": (0, 1)})):
            run(f"{S}/step5_rates.py", "--no-overlap", "--tree-config", cf, "--loci", f"{d}/l.tsv", "--genes", f"{d}/g.tsv",
                "--hogs", f"{d}/h.tsv", "--branch-lengths", f"{d}/len.tsv", "--hog-span", mode,
                "--outdir", f"{d}/{mode}")
            r = csv.DictReader(open(f"{d}/{mode}/rates_by_branch.tsv"), delimiter="\t")
            got = {x["Branch"]: (int(x["PDNG_births"]), int(x["PDNG_deaths"])) for x in r
                   if int(x["PDNG_births"]) or int(x["PDNG_deaths"])}
            assert got == want, (mode, got)


def test_overlap_classes():
    """Own genes found through their transcript IDs (gene-<id> with mRNA <id>), chromosome named Chr1 in the GFF and chr1
    in the coordinates. All on the minus strand:
      loc1, CDS 1001-1100 and 1401-1500, inside same-strand host1 (900-2000) whose CDS 1200-1300 is in loc1's intron
        -> same_strand_overlap (counted)
      loc3, CDS 9001-9200, and same-strand ov3 with CDS 9150-9300 -> same_strand_cds_overlap (left out)
      loc8, CDS 12001-12200, and plus-strand anti8 (12100-12500) with CDS 12400-12500 -> antisense_overlap (counted)"""
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        g = ["##gff-version 3"]
        for gid, st, s0, e0, cds, own in (("loc1", "-", 1001, 1500, [(1001, 1100), (1401, 1500)], True),
                                          ("host1", "-", 900, 2000, [(1200, 1300)], False),
                                          ("loc3", "-", 9001, 9200, [(9001, 9200)], True),
                                          ("ov3", "-", 9100, 9400, [(9150, 9300)], False),
                                          ("loc8", "-", 12001, 12200, [(12001, 12200)], True),
                                          ("anti8", "+", 12100, 12500, [(12400, 12500)], False)):
            gene_id, tid = (f"gene-{gid}", gid) if own else (gid, f"{gid}.t")
            g += [f"Chr1\tx\tgene\t{s0}\t{e0}\t.\t{st}\t.\tID={gene_id};gene_biotype=protein_coding",
                  f"Chr1\tx\tmRNA\t{s0}\t{e0}\t.\t{st}\t.\tID={tid};Parent={gene_id}"]
            g += [f"Chr1\tx\tCDS\t{a}\t{b}\t.\t{st}\t0\tParent={tid}" for a, b in cds]
        w(f"{d}/gff/A.gff", "\n".join(g) + "\n")
        w(f"{d}/gffs.tsv", "species\tgff\nA\tgff/A.gff\n")
        w(f"{d}/c.tsv", "geneID\tfocal_species\tfile_name\tspecies_name\tchromosome\tstart\tend\tORF_Node\t"
          "Synteny_Node\tEnablers_Nodes\tGene_type\n" + "".join(
              f"{x}\tA\tf\tA\tchr1\t{s0}\t{e0}\tA-specific\tN2\tNA\tPDNG\n"
              for x, s0, e0 in (("loc1", 1000, 1500), ("loc3", 9000, 9200), ("loc8", 12000, 12200))))
        run(f"{S}/step3_overlap.py", "--tree-config", cf, "--coords", f"{d}/c.tsv", "--gff-table", f"{d}/gffs.tsv",
            "--out", f"{d}/ov.tsv")
        got = [(r["orientation"], r["other_gene"]) for r in csv.DictReader(open(f"{d}/ov.tsv"), delimiter="\t")]
        assert got == [("same_strand_overlap", "host1"), ("same_strand_cds_overlap", "ov3"),
                       ("antisense_overlap", "anti8")], got


def test_own_gene_rules():
    """Gene G1.g (minus) has transcripts Locus1.v1.0 (CDS 1001-1100, 1401-1500) and Sib.v1.0 (CDS 1201-1300). Same-strand
    gene ov2 has CDS 1220-1280, inside Locus1's intron but in Sib's CDS. Locus1 is found through Locus1.v1.0 and only its
    own CDS count -> same_strand_overlap. Locus zzzz9999 overlaps G1.g but no gene carries its ID ->
    own_gene_not_found; with 1 of 2 loci unmatched the default --max-unmatched stops the run, and with 1 of 3 loci
    unknown so does the default --max-unknown. A .fai with a GFF lacking ##sequence-region lines stops the run.
    Chromosomes: GFF CM1 and CM2 (both 50000 bp) and CM3 (70000 bp); .fai NC1, NC2 (50000) and NC3 (70000). Only CM3
    is renamed (NC3); a locus on NC1 is unknown, not compared with genes of CM1 or CM2."""
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        g = ["##gff-version 3", "##sequence-region CM1 1 50000", "##sequence-region CM2 1 50000",
             "##sequence-region CM3 1 70000",
             "CM3\tx\tgene\t1001\t1500\t.\t-\t.\tID=G1.g;gene_biotype=protein_coding",
             "CM3\tx\tmRNA\t1001\t1500\t.\t-\t.\tID=Locus1.v1.0;Parent=G1.g",
             "CM3\tx\tCDS\t1001\t1100\t.\t-\t0\tParent=Locus1.v1.0", "CM3\tx\tCDS\t1401\t1500\t.\t-\t0\tParent=Locus1.v1.0",
             "CM3\tx\tmRNA\t1001\t1500\t.\t-\t.\tID=Sib.v1.0;Parent=G1.g", "CM3\tx\tCDS\t1201\t1300\t.\t-\t0\tParent=Sib.v1.0",
             "CM3\tx\tgene\t1150\t1350\t.\t-\t.\tID=ov2", "CM3\tx\tmRNA\t1150\t1350\t.\t-\t.\tID=ov2.t;Parent=ov2",
             "CM3\tx\tCDS\t1220\t1280\t.\t-\t0\tParent=ov2.t",
             "CM1\tx\tgene\t101\t900\t.\t+\t.\tID=Loc9", "CM2\tx\tgene\t101\t900\t.\t+\t.\tID=other9"]
        w(f"{d}/gff/A.gff", "\n".join(g) + "\n")
        w(f"{d}/gff/A.fai", "NC1\t50000\t1\t60\t61\nNC2\t50000\t1\t60\t61\nNC3\t70000\t1\t60\t61\n")
        w(f"{d}/gffs.tsv", "species\tgff\tfai\nA\tgff/A.gff\tgff/A.fai\n")
        hdr = "geneID\tfocal_species\tfile_name\tspecies_name\tchromosome\tstart\tend\tORF_Node\tSynteny_Node\t" \
              "Enablers_Nodes\tGene_type\n"
        rows = (("Locus1", "NC3", 1000, 1500), ("zzzz9999", "NC3", 1050, 1090), ("Loc9", "NC1", 100, 900))
        w(f"{d}/c.tsv", hdr + "".join(f"{x}\tA\tf\tA\t{c}\t{s0}\t{e0}\tA-specific\tN2\tNA\tPDNG\n" for x, c, s0, e0 in rows))
        s3 = [f"{S}/step3_overlap.py", "--tree-config", cf, "--coords", f"{d}/c.tsv", "--gff-table", f"{d}/gffs.tsv",
              "--out", f"{d}/ov.tsv"]
        run(*s3, fail="probably belongs to another assembly")
        run(*s3, "--max-unknown", "0.5", fail="probably differ in form")
        w(f"{d}/gff/B.gff", "\n".join(x for x in g if not x.startswith("##sequence-region")) + "\n")
        w(f"{d}/gffs_b.tsv", "species\tgff\tfai\nA\tgff/B.gff\tgff/A.fai\n")
        run(*s3[:6], f"{d}/gffs_b.tsv", "--out", f"{d}/ov_b.tsv", fail="no ##sequence-region lines")
        run(*s3, "--max-unknown", "0.5", "--max-unmatched", "0.6")
        got = [(r["geneID"], r["orientation"], r["other_gene"]) for r in csv.DictReader(open(f"{d}/ov.tsv"), delimiter="\t")]
        assert got == [("Locus1", "same_strand_overlap", "ov2"), ("zzzz9999", "own_gene_not_found", ""),
                       ("Loc9", "unknown", "")], got


def test_id_suffix_rules():
    """HOG-table IDs with a bare '.p' (a1.p, b1.p) match the loci a1, b1: the run succeeds with a1/b1 in H1. HOG-table
    IDs whose trailing number is part of the gene name (evm.model.c1.1, evm.model.c1.2) would become one key: stop."""
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        run(f"{S}/step1_coordinates.py", "--tree-config", cf, "--coords", f"{d}/coords", "--nodes", f"{d}/nodes",
            "--tblastn", f"{d}/tb", "--out", f"{d}/out/coords.tsv")
        run(f"{S}/step2_loci.py", "--tree-config", cf, "--coords", f"{d}/out/coords.tsv", "--outdir", f"{d}/out")
        s5 = [f"{S}/step5_rates.py", "--no-overlap", "--tree-config", cf, "--loci", f"{d}/out/loci.tsv", "--genes", f"{d}/out/genes.tsv",
              "--branch-lengths", f"{d}/len.tsv"]
        w(f"{d}/p.tsv", "HOG\tOG\tGene Tree Parent Clade\tA\tB\tC\tD\n"
          "N0.H1\tOG1\tn\ta1.p\tb1.p\t\t\nN0.H2\tOG2\tn\ta2, a3\t\t\t\nN0.H3\tOG3\tn\ta4, a5\t\t\t\n")
        run(*s5, "--hogs", f"{d}/p.tsv", "--outdir", f"{d}/r1")
        st = {x["geneID"]: x["HOG_unit"] for x in csv.DictReader(open(f"{d}/r1/gene_status.tsv"), delimiter="\t")}
        assert st["a1"] == st["b1"] == "N0.H1", st
        assert len(st) == 9, st
        w(f"{d}/e.tsv", "HOG\tOG\tGene Tree Parent Clade\tA\tB\tC\tD\n"
          "N0.H1\tOG1\tn\ta1, evm.model.c1.1\tb1\t\t\nN0.H2\tOG2\tn\ta2, a3, evm.model.c1.2\t\t\t\n")
        run(*s5, "--hogs", f"{d}/e.tsv", "--outdir", f"{d}/r2", fail="would join different IDs")


def test_overlap_required():
    """The overlap check cannot be skipped silently. Step 5 stops without --overlap or --no-overlap. Step 3 stops when a
    species of the coordinates (B) has no GFF, runs when B is named in --no-gff, and stops when a species is both in the
    GFF table and in --no-gff."""
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        run(f"{S}/step1_coordinates.py", "--tree-config", cf, "--coords", f"{d}/coords", "--nodes", f"{d}/nodes",
            "--tblastn", f"{d}/tb", "--out", f"{d}/out/coords.tsv")
        run(f"{S}/step2_loci.py", "--tree-config", cf, "--coords", f"{d}/out/coords.tsv", "--outdir", f"{d}/out")
        run(f"{S}/step5_rates.py", "--tree-config", cf, "--loci", f"{d}/out/loci.tsv", "--genes", f"{d}/out/genes.tsv",
            "--hogs", f"{d}/N0.tsv", "--branch-lengths", f"{d}/len.tsv", "--outdir", f"{d}/r",
            fail="one of the arguments --overlap --no-overlap is required")
        w(f"{d}/gff/A.gff", "##gff-version 3\nchr1\tx\tgene\t1\t10\t.\t+\t.\tID=zzzz1\n")
        w(f"{d}/gffs.tsv", "species\tgff\nA\tgff/A.gff\n")
        hdr = "geneID\tfocal_species\tfile_name\tspecies_name\tchromosome\tstart\tend\tORF_Node\tSynteny_Node\t" \
              "Enablers_Nodes\tGene_type\n"
        w(f"{d}/c.tsv", hdr + "zzzz1\tA\tf\tA\tchr1\t1\t10\tA-specific\tN2\tNA\tPDNG\n"
          "yy\tB\tf\tB\tchr1\t1\t10\tB-specific\tN2\tNA\tPDNG\n")
        s3 = [f"{S}/step3_overlap.py", "--tree-config", cf, "--coords", f"{d}/c.tsv", "--gff-table", f"{d}/gffs.tsv"]
        run(*s3, "--out", f"{d}/o1.tsv", fail="name them in --no-gff")
        run(*s3, "--no-gff", "A,B", "--out", f"{d}/o2.tsv", fail="also listed in --no-gff")
        run(*s3, "--no-gff", "B", "--out", f"{d}/o3.tsv")
        got = [(r["species"], r["geneID"]) for r in csv.DictReader(open(f"{d}/o3.tsv"), delimiter="\t")]
        assert got == [("A", "zzzz1")], got


def test_unknown_genomes():
    """HOG H2 of test_hog_lineage_split: an origin at N1 carried by C (y1) and A (y2), lost on B. With B listed as
    unknown for y1 the death on B is not counted; listing D (outside the origin's clade) or A (a carrier) changes
    nothing, and no listing adds a death."""
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        w(f"{d}/l.tsv", "ClusterID\tNode\tGene_type\tA\tB\tC\tD\n3\tN1\tPDNG\t\t\ty1\t\n4\tA-specific\tPDNG\ty2\t\t\t\n")
        w(f"{d}/g.tsv", "species\tgeneID\tORF_Node\tGene_type\tClusterID\tCluster_node\tpresent_in\n"
          "C\ty1\tN1\tPDNG\t3\tN1\tC\nA\ty2\tA-specific\tPDNG\t4\tA-specific\tA\n")
        w(f"{d}/h.tsv", "HOG\tOG\tGene Tree Parent Clade\tA\tB\tC\tD\nH2\tO2\tn\ty2\t\ty1\t\n")
        uh = "species\tgeneID\tunknown_in\n"
        for name, unk, want_b in (("none", None, 1), ("B", uh + "C\ty1\tB\n", 0), ("D", uh + "C\ty1\tD\n", 1),
                                  ("A", uh + "A\ty2\tA;D\n", 1)):
            extra = []
            if unk:
                w(f"{d}/u_{name}.tsv", unk)
                extra = ["--unknown", f"{d}/u_{name}.tsv"]
            run(f"{S}/step5_rates.py", "--no-overlap", "--tree-config", cf, "--loci", f"{d}/l.tsv", "--genes", f"{d}/g.tsv",
                "--hogs", f"{d}/h.tsv", "--branch-lengths", f"{d}/len.tsv", *extra, "--outdir", f"{d}/{name}")
            r = {x["Branch"]: x for x in csv.DictReader(open(f"{d}/{name}/rates_by_branch.tsv"), delimiter="\t")}
            deaths = {b: int(r[b]["PDNG_deaths"]) for b in r if int(r[b]["PDNG_deaths"])}
            assert deaths == ({"B": 1} if want_b else {}), (name, deaths)
            assert int(r["N1"]["PDNG_births"]) == 1, name


def test_orf_presence():
    """Step 1a on one alignment of focal A (31-codon CDS, ATG ... TAA): B carries the same ORF (present); C has a stop
    at codon 10, so its ORF covers 30 of 96 nt, below 80% (absent); D has bases at 9 of 96 CDS columns
    (unknown). Locus a2 has no alignment: B, C and D unknown. The stage 4 options must use the length ORF test."""
    with tempfile.TemporaryDirectory() as d:
        cf = build(d)
        cds = "ATG" + "GCT" * 30 + "TAA"
        c_row = cds[:27] + "TAA" + cds[30:]
        d_row = cds[:9] + "-" * (len(cds) - 9)
        w(f"{d}/aln/A_stitched/A_a1_c1_100_200.aln", f">A\n{cds}\n>B\n{cds}\n>C\n{c_row}\n>D\n{d_row}\n>N2\n{cds}\n")
        w(f"{d}/tb1a/A_cdngs_tblastn_2_summary.tsv", "geneID\tNode_tblastn\tA\tB\tC\tD\na1\tN2\t1\t1\t1\t0\n"
          "a2\tA-specific\t1\t0\t0\t0\n")
        opts = {"orf_test": "length", "orf_min_len": 80, "orf_min_cov": 80, "orf_min_in_frame": 50, "orf_count_stop": True}
        w(f"{d}/opts.json", json.dumps(opts))
        step = [f"{S}/step1a_orf_presence.py", "--tree-config", cf, "--tblastn", f"{d}/tb1a", "--alignments", f"{d}/aln"]
        run(*step, "--stage4-options", f"{d}/opts.json", "--outdir", f"{d}/pres")
        got = {r["geneID"]: (r["A"], r["B"], r["C"], r["D"])
               for r in csv.DictReader(open(f"{d}/pres/A_cdngs_tblastn_2_summary.tsv"), delimiter="\t")}
        assert got == {"a1": ("1", "1", "0", "0"), "a2": ("1", "0", "0", "0")}, got
        unk = {(r["species"], r["geneID"]): r["unknown_in"]
               for r in csv.DictReader(open(f"{d}/pres/unknown_presence.tsv"), delimiter="\t")}
        assert unk == {("A", "a1"): "D", ("A", "a2"): "B;C;D"}, unk
        w(f"{d}/opts_tb.json", json.dumps({"orf_test": "tblastn"}))
        run(*step, "--stage4-options", f"{d}/opts_tb.json", "--outdir", f"{d}/pres2", fail="does not use the length ORF test")


if __name__ == "__main__":
    for t in (test_toy, test_id_form_mismatch_stops, test_overlap_inside_long_gene, test_hog_lineage_split,
              test_overlap_classes, test_own_gene_rules,
              test_id_suffix_rules, test_overlap_required, test_unknown_genomes, test_orf_presence):
        t()
        print(f"{t.__name__} PASSED")
