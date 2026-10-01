"""
test_combine.py - checks for bin/combine.py. Run from anywhere:   python3 tests/test_combine.py
No search tools needed: the test writes run folders exactly as the modules do (gene_status.tsv, passed_ids.tsv, hits.tsv,
run.json, passed/<species>.faa, .hep_run) and the expected answers were worked out by hand.

Three runs on one species, Sp, in this order (the chain order):
  blast  (homology-blast, all 10 genes)         g1 HIT  g2 NO_HIT  g3 NO_HIT  g4 IN_PHYLOGENY_ONLY  g5 NOT_RUN  g6 HIT
                                                g7 NO_HIT  g8 NO_HIT  g9 HIT  g10 NO_HIT
  cdd    (domain-cdd, all-at-once, g1..g9)      g1 NO_HIT  g2 HIT  g3 NO_HIT  g4 NO_HIT  g5 NO_HIT  g6 NO_HIT
                                                g7 SPURIOUS (flag only: passes)  g8 NOT_RUN  g9 HIT      (g10 not given)
  jack   (homology-jackhmmer, chained after blast: receives blast's survivors g2 g3 g4 g7 g8 g10)
                                                g2 NO_HIT  g3 HIT  g4 NO_HIT  g7 NO_HIT  g8 NO_HIT  g10 NO_HIT
Expected verdicts:
  g1  REMOVED (blast); jack NOT_SEARCHED, dropped_by blast
  g2  REMOVED (cdd)
  g3  REMOVED (jack)
  g4  PASS
  g5  INCOMPLETE (blast NOT_RUN; jack NOT_SEARCHED because blast stopped it)
  g6  REMOVED (blast); jack NOT_SEARCHED
  g7  PASS (SPURIOUS is only a flag)
  g8  INCOMPLETE (cdd NOT_RUN, nobody removed it)
  g9  REMOVED (blast + cdd); jack NOT_SEARCHED
  g10 INCOMPLETE (cdd never received it and no earlier run explains why)
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
SCRIPT = os.path.join(ROOT, "bin", "combine.py")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import combine as c  # noqa: E402
import hep_common as h  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + ((" | " + str(extra)) if (extra and not cond) else ""))
    if not cond:
        fails.append(name)


tmp = tempfile.mkdtemp()
SEQ = {"g%d" % i: "MKV" + "A" * i for i in range(1, 11)}


def make_run(path, module, label, statuses, hits=(), position=None, spurious_pass=(), species="Sp", info_extra=None):
    os.makedirs(os.path.join(path, "passed"), exist_ok=True)
    open(os.path.join(path, ".hep_run"), "w").write("module=%s\nlabel=%s\n" % (module, label))
    passed = [g for g, s in statuses.items() if s in (h.NO_HIT, h.EXCLUDED_ONLY, h.IN_PHYLOGENY_ONLY) or g in spurious_pass]
    with open(os.path.join(path, "gene_status.tsv"), "w") as fh:
        fh.write("species\tgene_id\tstatus\tn_counted_hits\n")
        for g, s in statuses.items():
            fh.write("%s\t%s\t%s\t%d\n" % (species, g, s, 1 if s == h.HIT else 0))
    with open(os.path.join(path, "passed_ids.tsv"), "w") as fh:
        fh.write("species\tgene_id\n")
        for g in passed:
            fh.write("%s\t%s\n" % (species, g))
    with open(os.path.join(path, "passed", species + ".faa"), "w") as fh:
        for g in passed:
            fh.write(">%s original header text\n%s\n" % (g, SEQ[g]))
    with open(os.path.join(path, "hits.tsv"), "w") as fh:
        fh.write("\t".join(h.HITS_COLUMNS) + "\n")
        for g, target, tsp, counts in hits:
            fh.write("\t".join([species, g, target, tsp, label, "1e-10", "100", "90.0", "50.0", "NA", "yes" if counts else "no", "NA"]) + "\n")
    info = {"module": module, "label": label, "chain": ({"sequential": True, "position": position, "of": 3} if position else None),
            "tools": [{"name": module + " tool", "version": "1.0"}], "date": "2026-09-26"}
    info.update(info_extra or {})
    json.dump(info, open(os.path.join(path, "run.json"), "w"))
    return path


def table(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return [dict(zip(rows[0], r)) for r in rows[1:]]


blast_st = {"g1": "HIT", "g2": "NO_HIT", "g3": "NO_HIT", "g4": "IN_PHYLOGENY_ONLY", "g5": "NOT_RUN", "g6": "HIT",
            "g7": "NO_HIT", "g8": "NO_HIT", "g9": "HIT", "g10": "NO_HIT"}
cdd_st = {"g1": "NO_HIT", "g2": "HIT", "g3": "NO_HIT", "g4": "NO_HIT", "g5": "NO_HIT", "g6": "NO_HIT", "g7": "SPURIOUS",
          "g8": "NOT_RUN", "g9": "HIT"}
jack_st = {"g2": "NO_HIT", "g3": "HIT", "g4": "NO_HIT", "g7": "NO_HIT", "g8": "NO_HIT", "g10": "NO_HIT"}
blast_hits = [("g1", "t1__Ceratodon_purpureus", "Ceratodon_purpureus", True), ("g1", "t2__Blasia_pusilla", "Blasia_pusilla", True),
              ("g4", "t3__Sp", "Sp", False),
              ("g6", "t4__Ceratodon_purpureus_R40", "Ceratodon_purpureus_R40", True),
              ("g6", "t5__Ceratodon_purpureus", "Ceratodon_purpureus", True),
              ("g9", "t6__Blasia_pusilla", "Blasia_pusilla", True)]
cdd_hits = [("g2", "cd00001", "NA", True), ("g9", "cd00002", "NA", True)]
jack_hits = [("g3", "t7__Fossombronia_cristula", "Fossombronia_cristula", True)]
R = os.path.join(tmp, "runs")
os.makedirs(R)
p_blast = make_run(os.path.join(R, "a_blast"), "homology-blast", "blast", blast_st, blast_hits, position=1)
p_cdd = make_run(os.path.join(R, "b_cdd"), "domain-cdd", "cdd", cdd_st, cdd_hits, spurious_pass=("g7",))
p_jack = make_run(os.path.join(R, "c_jack"), "homology-jackhmmer", "jack", jack_st, jack_hits, position=3)
def rc_c(rc, out):
    """combine exits 2 when some gene is INCOMPLETE (its pass cannot be claimed), otherwise 0."""
    inc = any(r["verdict"] == "INCOMPLETE" for r in table(os.path.join(out, "gene_matrix.tsv")))
    return rc == (2 if inc else 0)


counter = [0]


def run(args, expect_ok=True):
    counter[0] += 1
    out = os.path.join(tmp, "out%d" % counter[0])
    p = subprocess.run([sys.executable, SCRIPT, "-o", out, "--quiet"] + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       universal_newlines=True)
    return p, out


p, out = run(["--runs", p_blast, p_cdd, p_jack])
check("combine finishes with exit status 2, because 3 genes are INCOMPLETE (the tables are still written)", p.returncode == 2 and rc_c(p.returncode, out), p.stderr)
gm = {r["gene_id"]: r for r in table(os.path.join(out, "gene_matrix.tsv"))}
check("verdicts: g1 REMOVED g2 REMOVED g3 REMOVED g4 PASS g5 INCOMPLETE g6 REMOVED g7 PASS g8 INCOMPLETE g9 REMOVED g10 INCOMPLETE",
      {g: r["verdict"] for g, r in gm.items()} == {"g1": "REMOVED", "g2": "REMOVED", "g3": "REMOVED", "g4": "PASS",
                                                    "g5": "INCOMPLETE", "g6": "REMOVED", "g7": "PASS", "g8": "INCOMPLETE",
                                                    "g9": "REMOVED", "g10": "INCOMPLETE"}, {g: r["verdict"] for g, r in gm.items()})
check("a gene an earlier run removed is NOT_SEARCHED in the later chained run, never 0 or blank (g1, g6, g9)",
      all(gm[g]["jack"] == "NOT_SEARCHED" for g in ("g1", "g6", "g9")) and gm["g1"]["blast"] == "HIT")
check("dropped_by names the run that stopped the gene: g1 blast; g2 cdd; g5 blast (NOT_RUN also stops a chain); g4 NA",
      gm["g1"]["dropped_by"] == "blast" and gm["g2"]["dropped_by"] == "cdd" and gm["g5"]["dropped_by"] == "blast"
      and gm["g4"]["dropped_by"] == "NA", {g: gm[g]["dropped_by"] for g in gm})
check("removed_by lists every run that removed it (g9: blast;cdd) and NA for none",
      gm["g9"]["removed_by"] == "blast;cdd" and gm["g3"]["removed_by"] == "jack" and gm["g4"]["removed_by"] == "NA")
check("g5: blast NOT_RUN and jack NOT_SEARCHED, so the pass cannot be claimed", gm["g5"]["blast"] == "NOT_RUN" and gm["g5"]["verdict"] == "INCOMPLETE")
check("g7: SPURIOUS is only a flag, so the gene passes", gm["g7"]["cdd"] == "SPURIOUS" and gm["g7"]["verdict"] == "PASS")
check("g10: never given to cdd and no earlier run explains why: NOT_SEARCHED with dropped_by NA, INCOMPLETE, and a warning",
      gm["g10"]["cdd"] == "NOT_SEARCHED" and gm["g10"]["dropped_by"] == "NA" and gm["g10"]["verdict"] == "INCOMPLETE"
      and "g10" in p.stderr, (gm["g10"], p.stderr))
check("breadth: g1 hits 2 species and 2 genera; g6 hits 2 names of ONE species (strain) and 1 genus; g9 1 and 1",
      (gm["g1"]["n_hit_species"], gm["g1"]["n_hit_genera"]) == ("2", "2")
      and (gm["g6"]["n_hit_species"], gm["g6"]["n_hit_genera"]) == ("1", "1")
      and (gm["g9"]["n_hit_species"], gm["g9"]["n_hit_genera"]) == ("1", "1"), {g: (gm[g]["n_hit_species"], gm[g]["n_hit_genera"]) for g in gm})
check("breadth: 0 where a homology run searched the gene and found nothing that counts (g4: its only hit is in-phylogeny, g2)",
      (gm["g4"]["n_hit_species"], gm["g2"]["n_hit_species"]) == ("0", "0"))
check("breadth: a domain hit is not a species (g2 hit cdd only); g5 was searched by no homology run, so NA, not 0",
      gm["g5"]["n_hit_species"] == "NA" and gm["g5"]["n_hit_genera"] == "NA")
check("n_runs_searched counts the runs that really decided: g1 2 (blast, cdd), g10 2 (blast, jack), g5 1 (blast NOT_RUN)", gm["g1"]["n_runs_searched"] == "2" and gm["g10"]["n_runs_searched"] == "2"
      and gm["g5"]["n_runs_searched"] == "1", {g: gm[g]["n_runs_searched"] for g in gm})
check("no table has an empty cell",
      all(v != "" for f in os.listdir(out) if f.endswith(".tsv") for r in table(os.path.join(out, f)) for v in r.values()))
rs = {r["run"]: r for r in table(os.path.join(out, "run_summary.tsv"))}
check("run_summary: blast received 10, removed 3 (g1 g6 g9), passed 6, NOT_RUN 1 (g5); all 3 removed here first",
      (rs["blast"]["genes_received"], rs["blast"]["removed"], rs["blast"]["passed"], rs["blast"]["not_run"],
       rs["blast"]["removed_first"]) == ("10", "3", "6", "1", "3"), rs["blast"])
check("run_summary: cdd removed g2 and g9; only g2 is new (g9 was already removed by blast): removed_first 1",
      (rs["cdd"]["removed"], rs["cdd"]["removed_first"], rs["cdd"]["not_run"], rs["cdd"]["not_searched"]) == ("2", "1", "1", "1"), rs["cdd"])
check("run_summary: jack received 6, removed g3 (new), and 4 genes were never given to it",
      (rs["jack"]["genes_received"], rs["jack"]["removed"], rs["jack"]["removed_first"], rs["jack"]["not_searched"]) == ("6", "1", "1", "4"), rs["jack"])
ov = {r["removed_by"]: r["n_genes"] for r in table(os.path.join(out, "overlaps.tsv"))}
check("overlaps: blast alone 2 (g1 g6), blast+cdd 1 (g9), cdd alone 1 (g2), jack alone 1 (g3)",
      ov == {"blast": "2", "blast+cdd": "1", "cdd": "1", "jack": "1"}, ov)
sm = {r["species"]: r for r in table(os.path.join(out, "summary.tsv"))}
check("summary: 10 genes = 2 PASS + 5 REMOVED + 3 INCOMPLETE",
      (sm["ALL"]["genes"], sm["ALL"]["PASS"], sm["ALL"]["REMOVED"], sm["ALL"]["INCOMPLETE"]) == ("10", "2", "5", "3"), sm["ALL"])
ids = [r["gene_id"] for r in table(os.path.join(out, "passed_all_ids.tsv"))]
check("passed_all_ids lists exactly the PASS genes", sorted(ids) == ["g4", "g7"], ids)
fa = open(os.path.join(out, "passed_all", "Sp.faa")).read()
check("passed_all FASTA has those two genes with their original headers and sequences",
      fa.count(">") == 2 and ">g4 original header text\n" + SEQ["g4"] in fa and ">g7 original header text" in fa, fa)
rj = json.load(open(os.path.join(out, "run.json")))
check("run.json records the run order, the tools and the counts", rj["run_order"] == ["blast", "cdd", "jack"]
      and rj["counts"]["PASS"] == 2 and len(rj["runs"]) == 3, rj)
check("the end-of-run note is written and names the tools used",
      "blast tool" in open(os.path.join(out, "citations.txt")).read().replace("homology-blast tool", "blast tool"))

# --- binary tables
p, out = run(["--runs", p_blast, p_cdd, p_jack, "--binary"])
bm = {r["gene_id"]: r for r in table(os.path.join(out, "gene_matrix_binary.tsv"))}
check("binary: 1 = the run removed it, 0 = it searched and did not, NA = not searched or not finished (never 0)",
      (bm["g1"]["blast"], bm["g1"]["cdd"], bm["g1"]["jack"]) == ("1", "0", "NA")
      and (bm["g5"]["blast"], bm["g5"]["jack"]) == ("NA", "NA") and (bm["g8"]["cdd"], bm["g8"]["jack"]) == ("NA", "0")
      and (bm["g7"]["cdd"], bm["g4"]["blast"]) == ("0", "0"), bm)
p, out = run(["--runs", p_blast, p_cdd, p_jack, "--binary", "--na", "-1"])
bm = {r["gene_id"]: r for r in table(os.path.join(out, "gene_matrix_binary.tsv"))}
check("--na changes the placeholder (here -1), for scripts that want a number", bm["g1"]["jack"] == "-1" and bm["g5"]["blast"] == "-1", bm["g1"])
p, out = run(["--runs", p_blast, "--na", " "])
check("a blank --na is refused", p.returncode != 0 and "blank" in p.stderr, p.stderr)
p, out = run(["--runs", p_blast, p_cdd, p_jack])
check("the binary table is not written unless asked", not os.path.exists(os.path.join(out, "gene_matrix_binary.tsv")))

# --- a folder of runs is read in chain order
p, out = run(["--runs", R])
gm2 = table(os.path.join(out, "gene_matrix.tsv"))
check("a folder of run folders is expanded and put in chain order (position 1, then the unnumbered run, then position 3 last)",
      rc_c(p.returncode, out) and list(gm2[0].keys())[2:5] == ["blast", "jack", "cdd"], list(gm2[0].keys()))

# runs that were chained by hand (no --sequential): the folder's alphabetical order must NOT decide; the date does
D = os.path.join(tmp, "byhand")
os.makedirs(D)
make_run(os.path.join(D, "z_first"), "homology-blast", "first", blast_st, blast_hits, info_extra={"date": "2026-09-26T10:00:00"})
make_run(os.path.join(D, "a_second"), "homology-jackhmmer", "second", jack_st, jack_hits, info_extra={"date": "2026-09-26T11:00:00"})
p, out = run(["--runs", D])
gmh = table(os.path.join(out, "gene_matrix.tsv"))
check("runs chained by hand are read in the order they were made, not alphabetically (z_first was made before a_second)",
      rc_c(p.returncode, out) and list(gmh[0].keys())[2:4] == ["first", "second"] and "missing from run" not in p.stderr, (list(gmh[0].keys()), p.stderr))
REAL_DIR = os.path.join(ROOT, "example_data", "combine", "real_sample")
p, out = run(["--runs", REAL_DIR])
check("the real sample folder read as a folder: chain order found from the dates, no false warnings, same result",
      p.returncode == 0 and "missing from run" not in p.stderr
      and {r["species"]: r["PASS"] for r in table(os.path.join(out, "summary.tsv"))}["ALL"] == "12", p.stderr)

# --- what is refused
p, out = run(["--runs", os.path.join(tmp, "nowhere")])
check("a missing folder is refused", p.returncode != 0 and "not a folder" in p.stderr, p.stderr)
os.makedirs(os.path.join(tmp, "plain"))
p, out = run(["--runs", os.path.join(tmp, "plain")])
check("a folder that is not a run is refused", p.returncode != 0 and "no run folders" in p.stderr, p.stderr)
dup = make_run(os.path.join(tmp, "dup_blast"), "homology-blast", "blast", blast_st, blast_hits)
p, out = run(["--runs", p_blast, dup])
check("two runs with the same label are refused", p.returncode != 0 and "same label" in p.stderr, p.stderr)
bad = make_run(os.path.join(tmp, "badlabel"), "homology-blast", "verdict", blast_st, blast_hits)
p, out = run(["--runs", bad])
check("a run label equal to a column name is refused", p.returncode != 0 and "column name" in p.stderr, p.stderr)
incons = make_run(os.path.join(tmp, "incons"), "homology-blast", "inc", {"g1": "NO_HIT", "g2": "NO_HIT"}, [])
with open(os.path.join(incons, "passed_ids.tsv"), "w") as fh:
    fh.write("species\tgene_id\nSp\tg2\n")
p, out = run(["--runs", incons])
check("run files that disagree (NO_HIT but not in passed_ids) are refused, not guessed at",
      p.returncode != 0 and "disagree" in p.stderr, p.stderr)
half = make_run(os.path.join(tmp, "half"), "homology-blast", "half", {"g1": "NO_HIT"}, [])
os.remove(os.path.join(half, "run.json"))
p, out = run(["--runs", half])
check("a run folder with a missing file is refused", p.returncode != 0 and "missing" in p.stderr, p.stderr)
p, out1 = run(["--runs", p_blast])
p = subprocess.run([sys.executable, SCRIPT, "--runs", p_cdd, "-o", out1, "--quiet"], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                   universal_newlines=True)
check("running again into the same output folder is refused", p.returncode != 0 and "already holds results" in p.stderr, p.stderr)
p = subprocess.run([sys.executable, SCRIPT, "--runs", p_cdd, "-o", out1, "--quiet", "--force"], stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, universal_newlines=True)
check("--force overwrites it", rc_c(p.returncode, out1), p.stderr)
p, out2 = run(["--runs", out1])
check("the output of combine itself is not accepted as a run", p.returncode != 0, p.stderr)
lost = make_run(os.path.join(tmp, "lost"), "homology-blast", "lost", {"g1": "NO_HIT"}, [])
os.remove(os.path.join(lost, "passed", "Sp.faa"))
p, out = run(["--runs", lost])
check("a passed gene whose sequence cannot be found is an error, never silently dropped",
      p.returncode != 0 and "sequences" in p.stderr, p.stderr)

# --- real run folders (example_data/combine/real_sample: 200 real moss proteins; blastp against 7 real proteomes, then jackhmmer
# on BLAST's 13 survivors (chained), and jackhmmer on all 200 (all at once); hits.tsv trimmed to 5 rows per gene)
REAL = os.path.join(ROOT, "example_data", "combine", "real_sample")
p, out = run(["--runs", os.path.join(REAL, "homology-blast"), os.path.join(REAL, "chain_jack"), os.path.join(REAL, "homology-jackhmmer")])
rg = {r["gene_id"]: r for r in table(os.path.join(out, "gene_matrix.tsv"))}
rsm = {r["species"]: r for r in table(os.path.join(out, "summary.tsv"))}["ALL"]
check("real runs: 200 genes = 12 PASS + 188 REMOVED + 0 INCOMPLETE", p.returncode == 0
      and (rsm["genes"], rsm["PASS"], rsm["REMOVED"], rsm["INCOMPLETE"]) == ("200", "12", "188", "0"), (p.stderr, rsm))
nsr = [g for g, r in rg.items() if r["chain_jack"] == "NOT_SEARCHED"]
check("real runs: the 187 genes BLAST removed are NOT_SEARCHED in the chained run, dropped_by homology-blast",
      len(nsr) == 187 and all(rg[g]["dropped_by"] == "homology-blast" and rg[g]["homology-blast"] == "HIT" for g in nsr))
rrs = {r["run"]: r for r in table(os.path.join(out, "run_summary.tsv"))}
check("real runs: the chained run received 13 genes and removed 1; the all-at-once jackhmmer removed 188 (1 new after BLAST)",
      (rrs["chain_jack"]["genes_received"], rrs["chain_jack"]["removed"], rrs["homology-jackhmmer"]["removed"]) == ("13", "1", "188"), rrs)
check("real runs: no local user name or path leaks into the example data",
      not any("ludtson" in open(os.path.join(r, f), errors="replace").read() for r, _, fs in os.walk(REAL) for f in fs))

# --- other pieces
check("a single run works: every verdict is PASS or REMOVED or INCOMPLETE by that run alone",
      set(r["verdict"] for r in table(os.path.join(run(["--runs", p_blast])[1], "gene_matrix.tsv"))) == {"REMOVED", "PASS", "INCOMPLETE"})
two_sp = make_run(os.path.join(tmp, "sp2"), "homology-blast", "two", {"g1": "NO_HIT"}, [], species="Other")
sp2 = c.Run(two_sp)
check("the Run reader keeps species apart: the same gene ID in two species is two genes", ("Other", "g1") in sp2.status)

print()
print("FAILED: %d" % len(fails) if fails else "ALL PASSED")
sys.exit(1 if fails else 0)
