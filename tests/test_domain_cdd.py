"""
test_domain_cdd.py - checks for bin/domain_cdd.py using the mock data in example_data/domain-cdd/.

Run from anywhere:   python3 tests/test_domain_cdd.py
Needs no external tools and no databases (it uses --parse-only on a mock RPS-BLAST result).
The expected answers below were worked out by hand from the mock data, not copied from the program's output:

  gene_a1  Pfam kinase, covers 250 of 264 model positions (94.7%), E 1e-40          -> counts
  gene_a2  COG, covers 241 of 400 (60.25%), E 1e-8                                  -> counts
  gene_a3  SMART, covers 50 of 250 (20%)                                            -> too little of the domain
  gene_a4  searched, no hits                                                        -> no hit
  gene_a5  missing from the result file                                             -> NOT_RUN
  gene_a6  CD, E 0.05 (above the 0.01 threshold)                                    -> no hit
  gene_a7  Pfam, two alignments cover 45 and 46 of 100 alone, 91 together           -> counts (union)
  gene_a8  PRK (header carries the __SpeciesA suffix), covers 200 of 300 (66.7%)    -> counts
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
SCRIPT = os.path.join(ROOT, "bin", "domain_cdd.py")
DATA = os.path.join(ROOT, "example_data", "domain-cdd")
INPUT = os.path.join(DATA, "input")
MOCK = os.path.join(DATA, "mock_rpsblast")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import domain_cdd as d  # noqa: E402
import hep_common as h  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + ((" | " + str(extra)) if (extra and not cond) else ""))
    if not cond:
        fails.append(name)


tmp = tempfile.mkdtemp()
counter = [0]


def run(extra_args, mock=MOCK, inp=INPUT):
    """Run domain_cdd.py --parse-only; return (exit code, stderr, output dir)."""
    counter[0] += 1
    out = os.path.join(tmp, "run%d" % counter[0])
    cmd = [sys.executable, SCRIPT, "-i", inp, "-o", out, "--parse-only", mock, "--quiet"] + extra_args
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return p.returncode, p.stderr, out


def statuses(out):
    rows = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(out, "gene_status.tsv"))][1:]
    return {r[1]: r[2] for r in rows}

def rc_ok(rc, out):
    """The exit-status rule: 0 when every gene was searched, 2 when some gene is NOT_RUN (the results are still written)."""
    return rc == (2 if "NOT_RUN" in statuses(out).values() else 0)



def table(path):
    lines = [ln.rstrip("\n").split("\t") for ln in open(path)]
    return [dict(zip(lines[0], r)) for r in lines[1:]]


BASE = {"gene_a1": "HIT", "gene_a2": "HIT", "gene_a3": "NO_HIT", "gene_a4": "NO_HIT", "gene_a5": "NOT_RUN",
        "gene_a6": "NO_HIT", "gene_a7": "HIT", "gene_a8": "HIT"}

# --- default rules
rc, err, out = run([])
check("default run succeeds (exit 2 because the mock has a NOT_RUN gene; results written)", rc_ok(rc, out) and rc == 2, err)
check("default statuses", statuses(out) == BASE, statuses(out))


def blank_cells(folder):
    """(file, line) of every table row with an empty cell. The tool never writes blanks: missing values are NA."""
    bad = []
    for name in ("hits.tsv", "gene_status.tsv", "summary.tsv", "passed_ids.tsv"):
        for n, line in enumerate(open(os.path.join(folder, name)).read().splitlines(), 1):
            if "" in line.split("\t"):
                bad.append((name, n))
    return bad


check("no output table has an empty cell (missing values are written as NA, so scripts never see a shifted column)",
      blank_cells(out) == [], blank_cells(out))
check("passed genes are exactly those without a counted hit",
      open(os.path.join(out, "passed_ids.tsv")).read().split() == ["species", "gene_id", "SpeciesA", "gene_a3",
                                                                      "SpeciesA", "gene_a4", "SpeciesA", "gene_a6"])
passed_fa = open(os.path.join(out, "passed", "SpeciesA.faa")).read()
check("passed FASTA has only the passers", passed_fa.count(">") == 3 and ">gene_a3 " in passed_fa)
hits = {r["gene_id"]: r for r in table(os.path.join(out, "hits.tsv"))}
check("a hit reports CDD accession and source",
      hits["gene_a1"]["target"] == "pfam00069" and hits["gene_a1"]["source"] == "PFAM"
      and hits["gene_a2"]["source"] == "COG" and hits["gene_a3"]["source"] == "SMART"
      and hits["gene_a6"]["source"] == "CD" and hits["gene_a8"]["source"] == "PRK", hits)
check("two alignments to one model are merged into one hit with union coverage",
      len([1 for r in table(os.path.join(out, "hits.tsv")) if r["gene_id"] == "gene_a7"]) == 1
      and hits["gene_a7"]["tcov"] == "91.0" and hits["gene_a7"]["qcov"] == "70.5", hits["gene_a7"])
check("coverage of a normal hit", hits["gene_a1"]["tcov"] == "94.7" and hits["gene_a1"]["qcov"] == "87.0", hits["gene_a1"])
check("a rejected hit says why (domain coverage)", "domain coverage 20.0%" in hits["gene_a3"]["note"]
      and hits["gene_a3"]["counts_as_hit"] == "no", hits["gene_a3"])
check("a rejected hit says why (E-value)", "E-value above 0.01" in hits["gene_a6"]["note"], hits["gene_a6"])
check("gene ID with __Species suffix is matched by its canonical ID", "gene_a8" in statuses(out))
check("citations.txt says results were parsed, not searched",
      "parsed from an earlier run" in open(os.path.join(out, "citations.txt")).read())

sm = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(out, "summary.tsv"))]
check("summary.tsv counts match the answers worked out by hand (8 genes: 4 hit, 3 no hit, 1 not run, 3 passed)",
      sm[1] == ["SpeciesA", "8", "3", "4", "3", "0", "0", "0", "1"] and sm[-1][0] == "ALL", sm)
again = [sys.executable, SCRIPT, "-i", INPUT, "-o", out, "--parse-only", MOCK, "--quiet"]
p = subprocess.run(again, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
check("running again into the same output folder is refused, with the way out named",
      p.returncode != 0 and "already holds results" in p.stderr and "--force" in p.stderr, p.stderr)
p = subprocess.run(again + ["--force"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
check("--force overwrites the results of the same run", rc_ok(p.returncode, out) and statuses(out) == BASE, p.stderr)

# --- settings change the answers as intended
rc, err, out = run(["--ignore-sources", "COG"])
s = statuses(out)
check("ignore_sources=COG: gene_a2 becomes EXCLUDED_ONLY and passes",
      s["gene_a2"] == "EXCLUDED_ONLY" and s["gene_a1"] == "HIT" and s["gene_a8"] == "HIT", s)
check("ignored hit is still reported, with the reason",
      {r["gene_id"]: r for r in table(os.path.join(out, "hits.tsv"))}["gene_a2"]["note"] == "ignored source COG")
rc, err, out = run(["--ignore-sources", "cog,Prk"])
s = statuses(out)
check("ignore_sources is case-insensitive and takes several sources",
      s["gene_a2"] == "EXCLUDED_ONLY" and s["gene_a8"] == "EXCLUDED_ONLY", s)
rc, err, out = run(["--min-domain-cov", "90"])
check("min_domain_cov=90", statuses(out) == dict(BASE, gene_a2="NO_HIT", gene_a8="NO_HIT"), statuses(out))
rc, err, out = run(["--min-qcov", "80"])
check("min_qcov=80 (gene_a7 covers 70.5% of the query, gene_a8 covers 80.4%)",
      statuses(out) == dict(BASE, gene_a2="NO_HIT", gene_a7="NO_HIT"), statuses(out))
rc, err, out = run(["--evalue", "1e-10"])
check("evalue=1e-10 is applied again when parsing", statuses(out) == dict(BASE, gene_a2="NO_HIT", gene_a7="NO_HIT"),
      statuses(out))
rj = __import__("json").load(open(os.path.join(out, "run.json")))
check("run.json records the changed setting", rj["non_default_settings"]["evalue"]["value"] == 1e-10, rj["non_default_settings"])

# --- an unfinished result must never read as "no hit"
mock_lines = open(os.path.join(MOCK, "SpeciesA.rpsblast.tsv")).read().splitlines()
d1 = os.path.join(tmp, "trunc_all")
os.makedirs(d1)
open(os.path.join(d1, "SpeciesA.rpsblast.tsv"), "w").write("\n".join(mock_lines[:-1]) + "\n")   # footer removed
rc, err, out = run([], mock=d1)
check("result without its completion footer: every gene is NOT_RUN",
      rc == 2 and set(statuses(out).values()) == {"NOT_RUN"}, statuses(out))
check("...and a warning says so", "unfinished" in err, err)

# Real rpsblast writes '# RPSBLAST', '# Query:' before every query and ONE '# BLAST processed N queries' footer at the
# end of a run. Several finished runs joined into one file (one per chunk) therefore have several footers.
# split point: the header line of gene_a6's query. Queries before it are a1-a4 (4), after it a6-a8 (3).
k = mock_lines.index("# Query: gene_a6 mock candidate gene") - 1
first, second = mock_lines[:k], mock_lines[k:-1]      # the mock's own footer (last line) is dropped


def joined(name, lines):
    d = os.path.join(tmp, name)
    os.makedirs(d)
    open(os.path.join(d, "SpeciesA.rpsblast.tsv"), "w").write("\n".join(lines) + "\n")
    return d


rc, err, out = run([], mock=joined("two_ok", first + ["# BLAST processed 4 queries"] + second
                                   + ["# BLAST processed 3 queries"]))
check("two finished runs joined in one file are both trusted", rc_ok(rc, out) and statuses(out) == BASE, statuses(out))

rc, err, out = run([], mock=joined("second_cut", first + ["# BLAST processed 4 queries"] + second))
s = statuses(out)
check("the second run was cut short: the first run's answers stand, the second run's genes are NOT_RUN",
      s["gene_a1"] == "HIT" and s["gene_a3"] == "NO_HIT" and s["gene_a4"] == "NO_HIT" and s["gene_a5"] == "NOT_RUN"
      and s["gene_a6"] == s["gene_a7"] == s["gene_a8"] == "NOT_RUN", s)

rc, err, out = run([], mock=joined("first_cut", first + second + ["# BLAST processed 3 queries"]))
check("the first run was cut short and a finished run follows: the counts disagree, so no gene is trusted",
      set(statuses(out).values()) == {"NOT_RUN"}, statuses(out))

rc, err, out = run([], mock=joined("wrong_count", first + second + ["# BLAST processed 9 queries"]))
check("a footer whose query count does not match: no gene is trusted",
      set(statuses(out).values()) == {"NOT_RUN"}, statuses(out))

# --- REAL rpsblast 2.17.0+ output (from a run against CDD 3.21), kept as a small regression test.
# Real output has '# RPSBLAST'/'# Query:' before every query and one footer at the end; the first version of the
# parser assumed otherwise and marked every gene NOT_RUN. Answers below were worked out by hand from the rows:
#   XP_024356364.1, XP_024356385.1  no hits found                                              -> NO_HIT
#   XP_024356422.1  pfam08626, covers 841 of 1221 model positions (68.9%), E 3.25e-53         -> HIT
#   XP_024356430.2  pfam07766, E 0.005 passes, but covers 68 of 229 model positions (29.7%)   -> NO_HIT (too little domain)
REAL = os.path.join(DATA, "real_sample")
rc, err, out = run([], mock=os.path.join(REAL, "rpsblast"), inp=os.path.join(REAL, "input"))
check("real output: succeeds and no gene is NOT_RUN", rc == 0 and "NOT_RUN" not in statuses(out).values(), err)
check("real output: statuses", statuses(out) == {"XP_024356364.1": "NO_HIT", "XP_024356385.1": "NO_HIT",
                                                 "XP_024356422.1": "HIT", "XP_024356430.2": "NO_HIT"}, statuses(out))
rh = {r["gene_id"]: r for r in table(os.path.join(out, "hits.tsv"))}
check("real output: accession, source and coverage",
      rh["XP_024356422.1"]["target"] == "pfam08626" and rh["XP_024356422.1"]["source"] == "PFAM"
      and rh["XP_024356422.1"]["tcov"] == "68.9" and rh["XP_024356422.1"]["qcov"] == "59.5"
      and rh["XP_024356430.2"]["tcov"] == "29.7" and "domain coverage 29.7%" in rh["XP_024356430.2"]["note"], rh)
rc, err, out = run(["--min-domain-cov", "25"], mock=os.path.join(REAL, "rpsblast"), inp=os.path.join(REAL, "input"))
check("real output: lowering min_domain_cov to 25 makes the LETM1 hit count", statuses(out)["XP_024356430.2"] == "HIT",
      statuses(out))

# --- bad input is refused, not guessed at
d3 = os.path.join(tmp, "badcols")
os.makedirs(d3)
bad = [ln if ln.startswith("#") else "\t".join(ln.split("\t")[:10]) for ln in mock_lines]
open(os.path.join(d3, "SpeciesA.rpsblast.tsv"), "w").write("\n".join(bad) + "\n")
rc, err, out = run([], mock=d3)
check("wrong column count stops the run and says which format is expected",
      rc != 0 and "expected 15" in err and "outfmt" in err, err)
d4 = os.path.join(tmp, "empty")
os.makedirs(d4)
rc, err, out = run([], mock=d4)
check("no result file for a species: all its genes are NOT_RUN, and the exit status is 2", rc == 2 and set(statuses(out).values()) == {"NOT_RUN"})
p = subprocess.run([sys.executable, SCRIPT, "-i", INPUT, "-o", os.path.join(tmp, "nodb")], stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, universal_newlines=True)
check("--db is required unless --parse-only", p.returncode != 0 and "--db is required" in p.stderr, p.stderr)

# --- pieces
check("source_of", d.source_of("pfam00069") == "PFAM" and d.source_of("cd00001") == "CD"
      and d.source_of("COG0515") == "COG" and d.source_of("") == "NA")
check("union of overlapping and adjacent intervals",
      d._union_len([(1, 10), (5, 20), (21, 30), (40, 50)]) == 30 + 11)
check("union is order- and direction-independent", d._union_len([(20, 5), (1, 10)]) == 20)
cmd = d.rps_command("q.faa", "/db/Cdd", "o.tsv", {"evalue": 0.01, "comp_based_stats": 1}, 4)
check("rpsblast command line", cmd[:3] == ["rpsblast", "-query", "q.faa"] and "-evalue" in cmd
      and cmd[cmd.index("-evalue") + 1] == "0.01" and cmd[cmd.index("-outfmt") + 1].startswith("7 qseqid")
      and cmd[cmd.index("-num_threads") + 1] == "4" and cmd[cmd.index("-comp_based_stats") + 1] == "1", cmd)

# --- real rpsblast + real CDD, only if the environment provides them
db = os.environ.get("HEP_CDD_DB")
if db and h.find_tool("rpsblast"):
    print("(real-data check: rpsblast and CDD found; see tests/real_cdd_check.md)")
else:
    print("SKIP real rpsblast run: set HEP_CDD_DB=/path/to/Cdd and put rpsblast on PATH to include it")

print()
print("FAILED: %d" % len(fails) if fails else "ALL PASSED")
sys.exit(1 if fails else 0)
