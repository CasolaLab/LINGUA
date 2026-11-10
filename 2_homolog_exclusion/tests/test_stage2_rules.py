"""
test_stage2_rules.py - checks for the options added in 0.2.0 (the Casola Lab's stage-2 re-screen rules, 2026-11) and for the
hog-rule module. Run from anywhere:   python3 tests/test_stage2_rules.py
No search tools needed. Expected answers were worked out by hand and are written next to each check.
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
sys.path.insert(0, os.path.join(ROOT, "bin"))
import domain_cdd as cdd  # noqa: E402
import domain_interpro as ipr  # noqa: E402
import hep_common as h  # noqa: E402
import hog_rule as hr  # noqa: E402
import homology_blast as hb  # noqa: E402
import homology_jackhmmer as hj  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + ((" | " + str(extra)) if (extra and not cond) else ""))
    if not cond:
        fails.append(name)


tmp = tempfile.mkdtemp()
M = h.TaxonMatcher(["Arabidopsis_thaliana"], ["Tarenaya"])

# ---- intervals
check("merge_intervals joins touching and overlapping", h.merge_intervals([(5, 10), (11, 12), (1, 3), (2, 4)]) == [(1, 12)])
check("overlap_length", h.overlap_length([(1, 10), (20, 30)], [(5, 25)]) == 6 + 6)
check("format/parse round trip", h.parse_intervals(h.format_intervals([(3, 9), (1, 2)])) == [(1, 9)])
check("genus_of", h.genus_of("Homo_sapiens") == "homo" and h.genus_of("Homo sapiens") == "homo" and h.genus_of(None) is None)


# ---- homology-blast: min_genera, filter_evalue, single_hsp
def brow(q, sid, e, qs, qe, qlen=100):
    return {"qid": q, "sid": sid, "pident": 50.0, "length": qe - qs + 1, "qstart": qs, "qend": qe, "sstart": 1, "send": qe - qs + 1,
            "evalue": e, "bitscore": 50.0, "qlen": qlen, "slen": 200}


cfg = h.resolve_config("homology-blast")[0]
cfg.update(min_qcov=0, min_genera=2, filter_evalue=1e-5)
rows = [brow("g1", "a__Homo_sapiens", 1e-10, 1, 50), brow("g1", "b__Homo_neanderthalensis", 1e-10, 1, 50),   # one genus
        brow("g2", "c__Homo_sapiens", 1e-10, 1, 50), brow("g2", "d__Mus_musculus", 1e-10, 1, 50),            # two genera
        brow("g3", "e__Homo_sapiens", 1e-10, 1, 50), brow("g3", "f__Mus_musculus", 1e-4, 1, 50),               # 2nd above 1e-5
        brow("g4", "g__Tarenaya_hassleriana", 1e-30, 1, 50), brow("g4", "h__Mus_musculus", 1e-10, 1, 50)]      # Tarenaya ignored
hits, counted, _ = hb.evaluate(rows, cfg, "blastp", M, None, "blast")
check("blast min_genera=2: one genus (g1) not counted, two genera (g2) counted", counted == {"g2"}, counted)
check("blast: demoted hits carry a note", any(r["gene_id"] == "g1" and "min_genera" in (r.get("note") or "") for r in hits))
check("blast: q_intervals written", any(r["gene_id"] == "g2" and r.get("q_intervals") == "1-50" for r in hits))
cfg.update(min_genera=1, filter_evalue=None, min_qcov=50, qcov_mode="single_hsp")
rows = [brow("t1", "x__TE_one", 1e-5, 1, 30), brow("t1", "x__TE_one", 1e-5, 40, 70),     # union 61%, longest HSP 31%
        brow("t2", "y__TE_two", 1e-5, 1, 60)]                                              # one HSP 60%
_, counted, _ = hb.evaluate(rows, cfg, "blastp", M, None, "te")
check("blast qcov_mode=single_hsp: only one HSP >= 50% counts", counted == {"t2"}, counted)
cfg["qcov_mode"] = "union"
_, counted, _ = hb.evaluate(rows, cfg, "blastp", M, None, "te")
check("blast qcov_mode=union (default) counts both", counted == {"t1", "t2"}, counted)
try:
    c2 = dict(cfg); c2["filter_evalue"] = 1e-1; hb.check_settings(c2); ok = False
except h.HepError:
    ok = True
check("blast: filter_evalue looser than evalue is refused", ok)


# ---- homology-jackhmmer: bias rule and min_genera
def jrow(q, tid, e, score, bias, qs=1, qe=40):
    return {"tid": tid, "tlen": 300, "qid": q, "qlen": 80, "evalue": e, "bitscore": score, "bias": bias, "ievalue": e,
            "qstart": qs, "qend": qe, "sstart": 1, "send": qe - qs + 1}


cfg = h.resolve_config("homology-jackhmmer")[0]
cfg.update(min_qcov=0, max_bias_frac=0.5, min_genera=2, filter_evalue=1e-5)
rows = [jrow("j1", "a__Homo_sapiens", 1e-8, 40, 25), jrow("j1", "b__Mus_musculus", 1e-8, 40, 5),      # 1st biased -> one genus
        jrow("j2", "c__Homo_sapiens", 1e-8, 40, 5), jrow("j2", "d__Mus_musculus", 1e-8, 40, 19.9)]     # both fine
hits, counted, _ = hj.evaluate(rows, cfg, M, None, "jack")
check("jackhmmer bias >= 0.5 x score not counted, so j1 falls below 2 genera", counted == {"j2"}, counted)
check("jackhmmer: bias note", any("composition bias" in (r.get("note") or "") for r in hits))
weak = jrow("j3", "e__Homo_sapiens", 1e-8, 40, 1, 5, 30); weak["ievalue"] = 0.1
hits, _, _ = hj.evaluate([weak], cfg, M, None, "jack")
check("jackhmmer: region falls back to all domains when none passes alone", hits[0].get("q_intervals") == "5-30", hits[0])

# ---- domain-cdd: ignore_domains
cfg = h.resolve_config("domain-cdd")[0]
cfg.update(min_domain_cov=0, ignore_domains="pfam02713")
crow = lambda q, acc: {"qid": q, "sid": "gnl|CDD|1", "pident": 40.0, "qstart": 1, "qend": 50, "sstart": 1, "send": 50,
                       "evalue": 1e-4, "bitscore": 40.0, "qlen": 100, "slen": 60, "stitle": "%s, DUF, desc" % acc}
hits, counted, ign = cdd.evaluate([crow("c1", "pfam02713"), crow("c2", "pfam00069")], cfg)
check("cdd ignore_domains: pfam02713 hit does not count, pfam00069 does", counted == {"c2"} and ign == {"c1"}, (counted, ign))

# ---- domain-interpro: PANTHER family rule, require_evalue, exemption
fam = os.path.join(tmp, "fam.tsv")
open(fam, "w").write("family\tspecies\tn_genes\nPTHR1\tARATH\t3\nPTHR1\tBRANA\t2\nPTHR2\tARATH\t1\nPTHR2\tHUMAN\t1\n")
cfg = h.resolve_config("domain-interpro")[0]
cfg.update(panther_family_species=fam, panther_exclude_species="ARATH,BRANA,BRACM,THECC", max_evalue=1e-3,
           max_evalue_exempt="PANTHER", require_evalue=True)
irow = lambda q, an, sig, v: {"qid": q, "length": 100, "analysis": an, "sig": sig, "desc": "", "start": 5, "stop": 60, "score": v,
                              "status": "T", "ipr": "-", "ipr_desc": "-"}
rows = [irow("p1", "PANTHER", "PTHR1:SF2", "1e-50"),     # family only in excluded species -> excluded
        irow("p2", "PANTHER", "PTHR2", "0.5"),           # family has HUMAN; E 0.5 but PANTHER exempt -> counts
        irow("p3", "Pfam", "PF00001", "-"),              # require_evalue -> not counted
        irow("p4", "Pfam", "PF00002", "1e-5"),           # counts
        irow("p5", "PANTHER", "PTHR9", "1e-50")]         # family absent from the table -> counted, with a warning
hits, counted, excl, _ = ipr.evaluate(rows, cfg)
check("interpro PANTHER rule, exemption and require_evalue", counted == {"p2", "p4", "p5"} and "p1" in excl, (counted, excl))
check("interpro q_intervals", any(r.get("q_intervals") == "5-60" for r in hits))

# ---- hog-rule
hogs = os.path.join(tmp, "N0.tsv")
open(hogs, "w").write("HOG\tOG\tGene Tree Parent Clade\tSpA\tSpB\nN0.HOG1\tOG1\tn1\tL1, L2\tM1\nN0.HOG2\tOG2\tn2\tL3\tM2, M3\n")
H = hr.read_hogs(hogs)
check("read_hogs", H == {"L1": "N0.HOG1", "L2": "N0.HOG1", "M1": "N0.HOG1", "L3": "N0.HOG2", "M2": "N0.HOG2", "M3": "N0.HOG2"}, H)
al = {"L1": {"M1": [(10, 40)]},          # overlaps M1 region 30-60 by 11 aa -> PASS at 20
      "L2": {"M1": [(25, 70)]},          # overlaps by 31 aa -> REMOVED_HOG
      "L3": {"M2": [(1, 5)], "M3": [(1, 100)]}}   # M2 region far away; M3 has no coordinates -> REGION_UNKNOWN
dec = hr.decide(al, {"M1": [(30, 60)], "M2": [(50, 80)]}, {"M3"}, 20)
check("hog-rule decide", dec == {"L1": ("PASS", "M1", 11), "L2": ("REMOVED_HOG", "M1", 31), "L3": ("REGION_UNKNOWN", None, 0)}, dec)
check("hog-rule unknown_region=pass", hr.decide(al, {"M1": [(30, 60)]}, {"M3"}, 20, "pass")["L3"][0] == "PASS")
bl = os.path.join(tmp, "hog_blastp.tsv")
open(bl, "w").write("# BLASTP 2.17.0+\n# Query: L2\n# Fields: x\n# 2 hits found\n"
                    "L2\tM1\t1e-20\t80\t1\t46\t25\t70\nL2\tM2\t1e-20\t80\t1\t46\t1\t40\n# BLAST processed 1 queries\n")
al2, conf = hr.parse_hog_blast(bl, {("L2", "M1")}, 1e-3)
check("parse_hog_blast keeps only same-HOG pairs", al2 == {"L2": {"M1": [(25, 70)]}} and conf == {"L2"}, al2)

# ---- hog-rule end to end (--parse-only), on a combine folder written by hand
import json  # noqa: E402
inp = os.path.join(tmp, "in"); os.makedirs(inp)
open(os.path.join(inp, "SpA_final.faa"), "w").write(">L1\nMKVAAAA\n>L2\nMKVAAAC\n>L3\nMKVAAAD\n>L4\nMKVAAAE\n")
open(os.path.join(inp, "SpB_final.faa"), "w").write(">M1\nMKVCCCC\n>M2\nMKVCCCD\n>M3\nMKVCCCE\n")
open(hogs, "w").write("HOG\tOG\tGene Tree Parent Clade\tSpA\tSpB\nN0.HOG1\tOG1\tn1\tL1, L2\tM1\n"
                      "N0.HOG2\tOG2\tn2\tL3\tM2, M3\nN0.HOG3\tOG3\tn3\tL4\t\n")
run = os.path.join(tmp, "blast"); os.makedirs(run)
open(os.path.join(run, "hits.tsv"), "w").write("\t".join(h.HITS_COLUMNS) + "\n" + "".join(
    "\t".join({"species": "SpB", "gene_id": g, "counts_as_hit": "yes", "q_intervals": iv}.get(c, "NA") for c in h.HITS_COLUMNS) + "\n"
    for g, iv in (("M1", "30-60"), ("M2", "50-80"), ("M3", "NA"))))
comb = os.path.join(tmp, "combine"); os.makedirs(comb)
json.dump({"runs": [{"label": "blast", "path": run}]}, open(os.path.join(comb, "run.json"), "w"))
open(os.path.join(comb, "gene_matrix.tsv"), "w").write("species\tgene_id\tblast\tverdict\tremoved_by\n" + "".join(
    "%s\t%s\tX\t%s\t%s\n" % r for r in (("SpA", "L1", "PASS", "NA"), ("SpA", "L2", "PASS", "NA"), ("SpA", "L3", "PASS", "NA"),
                                          ("SpA", "L4", "PASS", "NA"), ("SpB", "M1", "REMOVED", "blast"),
                                          ("SpB", "M2", "REMOVED", "blast"), ("SpB", "M3", "REMOVED", "blast"))))
bl = os.path.join(tmp, "hog_blastp2.tsv")
open(bl, "w").write("# BLASTP 2.17.0+\n# Query: L1\n# Query: L2\n# Query: L3\n"
                    "L1\tM1\t1e-20\t80\t1\t31\t10\t40\nL2\tM1\t1e-20\t80\t1\t46\t25\t70\n"
                    "L3\tM2\t1e-20\t80\t1\t5\t1\t5\nL3\tM3\t1e-20\t80\t1\t99\t1\t100\n# BLAST processed 3 queries\n")
o = os.path.join(tmp, "hog_out")
res = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "hog_rule.py"), "-i", inp, "-o", o, "--combine", comb, "--hogs", hogs,
                      "--parse-only", bl, "--set", "evidence_runs=blast", "--quiet"], capture_output=True, text=True)
got = {}
if os.path.isfile(os.path.join(o, "hog_rule.tsv")):
    for line in open(os.path.join(o, "hog_rule.tsv")).read().splitlines()[1:]:
        c = line.split("\t"); got[c[1]] = (c[4], c[6])
check("hog-rule end to end: verdicts", got == {"L1": ("PASS", "11"), "L2": ("REMOVED_HOG", "31"), "L3": ("REGION_UNKNOWN", "0"),
                                               "L4": ("PASS", "0")}, (got, res.stderr[-400:]))
check("hog-rule end to end: exit status 2 when a verdict is REGION_UNKNOWN", res.returncode == 2, res.returncode)
pas = open(os.path.join(o, "passed_ids.tsv")).read().split() if os.path.isfile(os.path.join(o, "passed_ids.tsv")) else []
bl3 = os.path.join(tmp, "hog_blastp3.tsv")          # L2 missing from the BLAST output
open(bl3, "w").write("# BLASTP 2.17.0+\n# Query: L1\n# Query: L3\nL1\tM1\t1e-20\t80\t1\t31\t10\t40\n"
                     "L3\tM2\t1e-20\t80\t1\t5\t1\t5\n# BLAST processed 2 queries\n")
o3 = os.path.join(tmp, "hog_out3")
res3 = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "hog_rule.py"), "-i", inp, "-o", o3, "--combine", comb, "--hogs", hogs,
                       "--parse-only", bl3, "--set", "evidence_runs=bl*", "--set", "unknown_region=pass", "--quiet"], capture_output=True, text=True)
got3 = {l.split("\t")[1]: l.split("\t")[4] for l in open(os.path.join(o3, "hog_rule.tsv")).read().splitlines()[1:]} \
    if os.path.isfile(os.path.join(o3, "hog_rule.tsv")) else res3.stderr[-300:]
check("hog-rule: candidate missing from the BLAST output is NOT_RUN (wildcard evidence_runs bl*)",
      got3 == {"L1": "PASS", "L2": "NOT_RUN", "L3": "PASS", "L4": "PASS"} and res3.returncode == 2, (got3, res3.returncode))
check("hog-rule end to end: passed list", pas == ["species", "gene_id", "SpA", "L1", "SpA", "L4"], pas)

# ---- paralog-rule
import paralog_rule as pr  # noqa: E402
# C1..C4 are candidates (C3 is a candidate paralog, not judged); O1, O2 are established genes of the same proteome
rows = [["C1", "C1", "100", "1e-50", "200", "1", "100", "100"],     # self hit: ignored
        ["C1", "O1", "60", "1e-20", "90", "1", "60", "100"],       # 60% identity over 60%: removes C1
        ["C2", "C3", "95", "1e-40", "150", "1", "100", "100"],     # paralog is a candidate: ignored
        ["C2", "O2", "45", "1e-10", "50", "1", "100", "100"],      # identity below 50: C2 passes, O2 reported
        ["C4", "O1", "80", "1e-12", "60", "1", "40", "100"],       # 40% of length: not enough
        ["C4", "O2", "90", "1e-3", "40", "1", "100", "100"]]       # E above 1e-5: ignored
al = pr.best_paralogs(rows, {"C1", "C2", "C3", "C4"}, 1e-5)
check("paralog-rule ignores self, candidate and weak-E hits", al == {"C1": [("O1", 60.0, 60.0, 1e-20)], "C2": [("O2", 45.0, 100.0, 1e-10)],
                                                                    "C4": [("O1", 80.0, 40.0, 1e-12)]}, al)
d = pr.decide(al, 50, 50)
check("paralog-rule decide", {k: v[:2] for k, v in d.items()} == {"C1": ("REMOVED_PARALOG", "O1"), "C2": ("PASS", "O2"),
                                                                    "C4": ("PASS", "O1")}, d)
check("paralog-rule: identity and coverage must hold in one alignment",
      pr.decide({"X": [("A", 90.0, 30.0, 1e-9), ("B", 40.0, 90.0, 1e-9)]}, 50, 50)["X"][0] == "PASS")
pin, pprot, pprev, pbl = (os.path.join(tmp, x) for x in ("par_in", "par_prot", "par_prev", "par_bl"))
for x in (pin, pprot, pprev, pbl):
    os.makedirs(x)
open(os.path.join(pin, "SpA_final.faa"), "w").write(">C1\nMKVA\n>C2\nMKVC\n>C3\nMKVD\n>C4\nMKVE\n>C5\nMKVF\n")
open(os.path.join(pprot, "SpA.faa"), "w").write(">C1\nMKVA\n>C2\nMKVC\n>C3\nMKVD\n>C4\nMKVE\n>C5\nMKVF\n>O1\nMKVG\n>O2\nMKVH\n")
open(os.path.join(pprev, "passed_ids.tsv"), "w").write("species\tgene_id\nSpA\tC1\nSpA\tC2\nSpA\tC4\nSpA\tC5\n")
open(os.path.join(pbl, "SpA.tsv"), "w").write("# BLASTP 2.17.0+\n# Query: C1\n# Query: C2\n# Query: C4\n"
                                              + "".join("\t".join(r) + "\n" for r in rows) + "# BLAST processed 3 queries\n")
po = os.path.join(tmp, "par_out")                                # C5 has no BLAST record: NOT_RUN, exit status 2
res = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "paralog_rule.py"), "-i", pin, "-o", po, "--passed", pprev,
                      "--proteomes", pprot, "--parse-only", pbl, "--quiet"], capture_output=True, text=True)
got = {l.split("\t")[1]: l.split("\t")[2] for l in open(os.path.join(po, "paralog_rule.tsv")).read().splitlines()[1:]} \
    if os.path.isfile(os.path.join(po, "paralog_rule.tsv")) else res.stderr[-300:]
check("paralog-rule end to end: verdicts", got == {"C1": "REMOVED_PARALOG", "C2": "PASS", "C4": "PASS", "C5": "NOT_RUN"}, got)
check("paralog-rule end to end: exit status 2 when a gene is NOT_RUN", res.returncode == 2, res.returncode)
pas = open(os.path.join(po, "passed_ids.tsv")).read().split() if os.path.isfile(os.path.join(po, "passed_ids.tsv")) else []
check("paralog-rule end to end: passed list", pas == ["species", "gene_id", "SpA", "C2", "SpA", "C4"], pas)

# ---- panther-families
pf =os.path.join(tmp, "PTHR19.0_test")
open(pf, "w").write("ARATH|TAIR=AT1G1|UniProtKB=Q1\tQ1\tAT1G1\tPTHR10:SF1\tX\n"
                    "ARATH|TAIR=AT1G2|UniProtKB=Q2\tQ2\tAT1G2\tPTHR10:SF2\tX\nHUMAN|HGNC=1|UniProtKB=P1\tP1\tG\tPTHR11\tY\n")
out = subprocess.run([sys.executable, os.path.join(ROOT, "bin", "panther_families.py"), pf], capture_output=True, text=True)
check("panther_families table", out.stdout == "family\tspecies\tn_genes\nPTHR10\tARATH\t2\nPTHR11\tHUMAN\t1\n", out.stdout)

# ---- 0.3.2: min_genera counted across several databases (pooled run)
import combine as cb  # noqa: E402

pool_root = os.path.join(tmp, "pool_out")
os.makedirs(pool_root)
pfa = os.path.join(tmp, "SpA_pool.faa")
open(pfa, "w").write("".join(">G%d\nMKV%s\n" % (i, "ACDEFG"[i - 1]) for i in range(1, 7)))


def comp_run(label, rows, statuses, min_genera=1):
    d = os.path.join(pool_root, label)
    with h.ModuleOutput(d, "homology-blast", label) as o:
        o.add_species("SpA", pfa, [dict(gene_id=g, target=t, target_species=s, source=label, counts_as_hit=c)
                                    for g, t, s, c in rows], statuses)
        o.close({"settings": {"min_genera": min_genera}}, quiet=True)
    return d


# db1: G1 Alpha; G2 Alpha; G3 Alpha; G5 only an in-phylogeny hit; G6 Alpha + Epsilon
d1 = comp_run("t-db1", [("G1", "a1", "Alpha_one", True), ("G2", "a2", "Alpha_one", True), ("G3", "a3", "Alpha_two", True),
                        ("G5", "t1", "Tarenaya_hassleriana", False), ("G6", "a6", "Alpha_one", True),
                        ("G6", "e6", "Epsilon_one", True)],
              {"G1": h.HIT, "G2": h.HIT, "G3": h.HIT, "G4": h.NO_HIT, "G5": h.IN_PHYLOGENY_ONLY, "G6": h.HIT})
# db2: G1 Beta; G2 Alpha (same genus as db1); G3 and G6 not searched; G4 Gamma + Delta; G5 nothing
d2 = comp_run("t-db2", [("G1", "b1", "Beta_one", True), ("G2", "a2b", "Alpha_three", True),
                        ("G4", "g4", "Gamma_one", True), ("G4", "d4", "Delta_one", True)],
              {"G1": h.HIT, "G2": h.HIT, "G3": h.NOT_RUN, "G4": h.HIT, "G5": h.NO_HIT, "G6": h.NOT_RUN})
pooled = os.path.join(pool_root, "t")
h.pool_database_runs("homology-blast", "t", pooled, [d1, d2], [("SpA", pfa)], 2, {"settings": {"min_genera": 2}}, quiet=True)
got = {r["gene_id"]: r["status"] for r in h._read_tsv(os.path.join(pooled, "gene_status.tsv"))}
# by hand: G1 Alpha(db1)+Beta(db2) = 2 genera -> HIT, though each database alone has 1; G2 Alpha in both = 1 genus -> NO_HIT;
# G3 1 genus and not searched in db2 -> NOT_RUN; G4 2 genera in db2 -> HIT; G5 in-phylogeny only -> IN_PHYLOGENY_ONLY;
# G6 2 genera in db1 -> HIT although db2 never searched it (a removal stands)
check("pooled min_genera across databases: statuses",
      got == {"G1": "HIT", "G2": "NO_HIT", "G3": "NOT_RUN", "G4": "HIT", "G5": "IN_PHYLOGENY_ONLY", "G6": "HIT"}, got)
pas = sorted(r["gene_id"] for r in h._read_tsv(os.path.join(pooled, "passed_ids.tsv")))
check("pooled run: passed list", pas == ["G2", "G5"], pas)
notes = {(r["gene_id"], r["target"]): r["note"] for r in h._read_tsv(os.path.join(pooled, "hits.tsv"))}
check("pooled run: demoted hit says why", "below min_genera=2" in notes.get(("G2", "a2"), ""), notes.get(("G2", "a2")))
check("pooled run: components marked", all(h.read_marker(d).get(h.COMPONENT_KEY) == "t" for d in (d1, d2)))
pj = json.load(open(os.path.join(pooled, "run.json")))
check("pooled run: run.json records scope and parts",
      pj.get("min_genera_scope") == "all_databases" and [x["label"] for x in pj.get("pooled_from", [])] == ["t-db1", "t-db2"])
check("combine expands to the pooled run only", cb.expand_runs([pool_root]) == [pooled], cb.expand_runs([pool_root]))
try:
    cb.Run(d1)
    check("combine refuses a component run", False)
except h.HepError:
    check("combine refuses a component run", True)
d3 = comp_run("t-db3", [], {g: h.NO_HIT for g in ("G1", "G2", "G3", "G4", "G5", "G6")}, min_genera=2)
try:
    h.pool_database_runs("homology-blast", "t3", os.path.join(pool_root, "t3"), [d3], [("SpA", pfa)], 2, {}, quiet=True)
    check("pooling refuses a component evaluated with min_genera>1", False)
except h.HepError:
    check("pooling refuses a component evaluated with min_genera>1", True)
check("pooling_wanted: side by side, default scope", h.pooling_wanted({"min_genera": 2}, 2, False, False) is True)
check("pooling_wanted: per_database", h.pooling_wanted({"min_genera_scope": "per_database"}, 2, False, False) is False)
check("pooling_wanted: sequential and single database", not h.pooling_wanted({}, 2, True, False)
      and not h.pooling_wanted({}, 1, False, False))
try:
    h.pooling_wanted({"min_genera_scope": "everywhere"}, 2, False, False)
    check("pooling_wanted rejects an unknown scope", False)
except h.HepError:
    check("pooling_wanted rejects an unknown scope", True)

print("ALL PASSED" if not fails else "%d FAILED: %s" % (len(fails), ", ".join(fails)))
sys.exit(1 if fails else 0)
