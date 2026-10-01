"""
test_domain_interpro.py - checks for bin/domain_interpro.py.

Run from anywhere:   python3 tests/test_domain_interpro.py
Needs no InterProScan and no databases: it uses mock InterProScan results, a small REAL InterProScan output kept in
example_data/domain-interpro/real_sample/, and (Linux/WSL only) a stand-in interproscan.sh for the run path.
The expected answers were worked out by hand:

Mock data (example_data/domain-interpro/), default settings:
  gene_b1   Pfam                                   -> counts
  gene_b2   MobiDB-Lite + Coils only               -> EXCLUDED_ONLY (not domain evidence)
  gene_b3   AntiFam only                           -> SPURIOUS (a spurious gene model, not a domain)
  gene_b4   AntiFam + Pfam                         -> counts (a real domain outweighs the AntiFam flag)
  gene_b5   no rows                                -> NO_HIT (the run is proven complete)
  gene_b6   Hamap, column 9 is a SCORE (32.0)      -> counts
  gene_b7   SMART E 0.05, no InterPro accession    -> counts (no cut-off by default)
  gene_b8   PROSITE pattern only                   -> EXCLUDED_ONLY
  gene_b9   Pfam (header carries __SpeciesA)       -> counts
  gene_b10  CDD, no InterPro accession             -> counts (real domain hits often have none)
  gene_b11  Pfam PF05212                           -> counts (EXCLUDED_ONLY if PF05212 is in ignore_domains)
Real output (real_sample/, InterProScan 5.78, 6 proteins): five have real domain matches; UPI0010EB9177 has only AntiFam
and MobiDB-Lite rows -> SPURIOUS; UPI0000167F57 has 22 MobiDB-Lite rows and ONE real Pfam domain (PF10744) -> counts.
"""
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
SCRIPT = os.path.join(ROOT, "bin", "domain_interpro.py")
DATA = os.path.join(ROOT, "example_data", "domain-interpro")
INPUT, MOCK = os.path.join(DATA, "input"), os.path.join(DATA, "mock_interproscan")
REAL_IN = os.path.join(DATA, "real_sample", "input")
REAL_TSV = os.path.join(DATA, "real_sample", "interproscan")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import domain_interpro as d  # noqa: E402
import hep_common as h  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + ((" | " + str(extra)) if (extra and not cond) else ""))
    if not cond:
        fails.append(name)


tmp = tempfile.mkdtemp()
counter = [0]


def run(extra, mock=MOCK, inp=INPUT, env=None):
    counter[0] += 1
    out = os.path.join(tmp, "run%d" % counter[0])
    cmd = [sys.executable, SCRIPT, "-i", inp, "-o", out, "--parse-only", mock, "--quiet"] + extra
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
    return p.returncode, p.stderr, out


def table(path):
    lines = [ln.rstrip("\n").split("\t") for ln in open(path)]
    return [dict(zip(lines[0], r)) for r in lines[1:]]


def statuses(out):
    return {r["gene_id"]: r["status"] for r in table(os.path.join(out, "gene_status.tsv"))}

def rc_ok(rc, out):
    """The exit-status rule: 0 when every gene was searched, 2 when some gene is NOT_RUN (the results are still written)."""
    return rc == (2 if "NOT_RUN" in statuses(out).values() else 0)



def passed(out):
    return sorted(ln.split("\t")[1] for ln in open(os.path.join(out, "passed_ids.tsv")).read().splitlines()[1:])


BASE = {"gene_b1": "HIT", "gene_b2": "EXCLUDED_ONLY", "gene_b3": "SPURIOUS", "gene_b4": "HIT", "gene_b5": "NO_HIT",
        "gene_b6": "HIT", "gene_b7": "HIT", "gene_b8": "EXCLUDED_ONLY", "gene_b9": "HIT", "gene_b10": "HIT",
        "gene_b11": "HIT"}

# --- default rules on the mock data
rc, err, out = run([])
check("default run succeeds (exit 2 only if the mock has a NOT_RUN gene)", rc_ok(rc, out), err)
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
check("genes that pass: no counted domain (AntiFam-only passes under 'flag')",
      passed(out) == ["gene_b2", "gene_b3", "gene_b5", "gene_b8"], passed(out))
hits = {}
for r in table(os.path.join(out, "hits.tsv")):
    hits.setdefault(r["gene_id"], []).append(r)
check("a counted match reports analysis, signature, E-value, InterPro accession and query coverage",
      hits["gene_b1"][0]["source"] == "Pfam" and hits["gene_b1"][0]["target"] == "PF00001"
      and hits["gene_b1"][0]["evalue"] == "1e-20" and hits["gene_b1"][0]["qcov"] == "86.0"
      and hits["gene_b1"][0]["note"] == "InterPro IPR000001", hits["gene_b1"])
check("Hamap's column 9 is a score, so it is reported as a bit score, not an E-value",
      hits["gene_b6"][0]["bitscore"] == "32.021942" and hits["gene_b6"][0]["evalue"] == "NA", hits["gene_b6"])
check("a domain match with no InterPro accession still counts, and says so",
      hits["gene_b10"][0]["counts_as_hit"] == "yes" and hits["gene_b10"][0]["note"] == "no InterPro entry", hits["gene_b10"])
check("MobiDB-Lite and Coils rows are listed but do not count, with the reason",
      all(r["counts_as_hit"] == "no" for r in hits["gene_b2"]) and hits["gene_b2"][0]["note"] == "not domain evidence (MobiDBLite)",
      hits["gene_b2"])
check("an AntiFam row is listed as a spurious-model signal",
      hits["gene_b3"][0]["counts_as_hit"] == "no" and hits["gene_b3"][0]["note"].startswith("spurious gene-model signal (AntiFam)"))
check("the AntiFam row on a gene with a real domain does not stop the domain from counting",
      [r["counts_as_hit"] for r in hits["gene_b4"]] == ["no", "yes"], hits["gene_b4"])
check("gene ID with __Species suffix is matched by its canonical ID", "gene_b9" in statuses(out))
sm = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(out, "summary.tsv"))]
check("summary.tsv counts (11 genes: 7 hit, 1 no hit, 2 excluded-only, 1 spurious; 4 passed)",
      sm[1] == ["SpeciesA", "11", "4", "7", "1", "2", "1", "0", "0"], sm)
rj = json.load(open(os.path.join(out, "run.json")))
check("run.json records that the lookup was off and the results were parsed, not searched",
      "disabled" in rj["lookup"] and rj["assumed_complete"] is False
      and rj["tools"][0]["name"].startswith("InterProScan, results parsed"), rj["lookup"])

# --- settings change the answers as intended
rc, err, out = run(["--antifam-action", "remove"])
check("antifam_action=remove: the spurious gene no longer passes",
      statuses(out) == BASE and passed(out) == ["gene_b2", "gene_b5", "gene_b8"], passed(out))
rc, err, out = run(["--max-evalue", "1e-2"])
s = statuses(out)
check("max_evalue applies to E-value analyses (gene_b7 E 0.05 stops counting)...", s["gene_b7"] == "NO_HIT", s)
check("...but never to score-type analyses (Hamap gene_b6 still counts)", s["gene_b6"] == "HIT" and s["gene_b10"] == "HIT", s)
rc, err, out = run(["--ignore-domains", "PF05212"])
check("ignore_domains by Pfam accession", statuses(out) == dict(BASE, gene_b11="EXCLUDED_ONLY"), statuses(out))
rc, err, out = run(["--ignore-domains", "ipr006000"])
check("ignore_domains by InterPro accession (any case)", statuses(out)["gene_b11"] == "EXCLUDED_ONLY", statuses(out))
rc, err, out = run(["--exclude-analyses", ""])
s = statuses(out)
check("exclude_analyses emptied: MobiDB-Lite, Coils and PROSITE patterns now count",
      s["gene_b2"] == "HIT" and s["gene_b8"] == "HIT" and s["gene_b3"] == "SPURIOUS", s)

# --- InterProScan's TSV lists only proteins with matches: proof of completion is needed for NO_HIT
plain = os.path.join(tmp, "tsv_only")
os.makedirs(plain)
open(os.path.join(plain, "SpeciesA.interproscan.tsv"), "w").write(open(os.path.join(MOCK, "SpeciesA.interproscan.tsv")).read())
rc, err, out = run([], mock=plain)
s = statuses(out)
check("no evidence the run finished: only genes with a counted domain are decided; the rest are NOT_RUN",
      rc == 2 and {g: v for g, v in s.items() if v == "HIT"} == {g: "HIT" for g, v in BASE.items() if v == "HIT"}
      and all(s[g] == "NOT_RUN" for g in ("gene_b2", "gene_b3", "gene_b5", "gene_b8")), s)
check("...and a warning says why", "nothing shows that the earlier run finished" in err, err)
rc, err, out = run(["--assume-complete"], mock=plain)
check("--assume-complete states the run finished: full answers", statuses(out) == BASE, statuses(out))
check("...and run.json records that this was assumed", json.load(open(os.path.join(out, "run.json")))["assumed_complete"] is True)
open(os.path.join(plain, "SpeciesA.interproscan.status"), "w").write("hep_status=complete\n")
rc, err, out = run([], mock=plain)
check("this module's own status file counts as proof of completion", statuses(out) == BASE, statuses(out))
open(os.path.join(plain, "SpeciesA.interproscan.status"), "w").write("hep_status=incomplete\n")
rc, err, out = run([], mock=plain)
check("a status file saying 'incomplete' is not proof", statuses(out)["gene_b5"] == "NOT_RUN", statuses(out))
empty = os.path.join(tmp, "nothing")
os.makedirs(empty)
rc, err, out = run(["--assume-complete"], mock=empty)
check("no result file for the species: all genes NOT_RUN, even with --assume-complete",
      rc == 2 and set(statuses(out).values()) == {"NOT_RUN"}, statuses(out))

# --- bad input is refused, not guessed at
bad = os.path.join(tmp, "badcols")
os.makedirs(bad)
open(os.path.join(bad, "SpeciesA.interproscan.tsv"), "w").write("gene_b1\tPfam\tPF00001\n")
rc, err, out = run([], mock=bad)
check("wrong column count stops the run and says what is expected", rc != 0 and "expected 15" in err and "no header" in err, err)
hdr = os.path.join(tmp, "withheader")
os.makedirs(hdr)
lines = open(os.path.join(MOCK, "SpeciesA.interproscan.tsv")).read().splitlines()
open(os.path.join(hdr, "SpeciesA.interproscan.tsv"), "w").write("\t".join("col%d" % i for i in range(15)) + "\n" + "\n".join(lines) + "\n")
rc, err, out = run([], mock=hdr)
check("a header line is refused rather than misread", rc != 0 and "cannot read a number" in err, err)
extra = os.path.join(tmp, "extra")
os.makedirs(extra)
open(os.path.join(extra, "SpeciesA.interproscan.tsv"), "w").write(
    open(os.path.join(MOCK, "SpeciesA.interproscan.tsv")).read()
    + "\t".join(["stranger", "0" * 32, "50", "Pfam", "PF00009", "x", "1", "40", "1E-9", "T", "26-09-2026", "IPR1", "x", "-", "-"]) + "\n")
open(os.path.join(extra, "SpeciesA.interproscan.log"), "w").write("100% done: InterProScan analyses completed\n")
rc, err, out = run([], mock=extra)
check("a protein in the results but not in the input is ignored with a warning",
      rc == 0 and "stranger" in err and statuses(out) == BASE, err)

# --- REAL InterProScan 5.78 output
REAL = {"O31533": "HIT", "UPI0010EB9177": "SPURIOUS", "UPI0004FABBC5": "HIT", "UPI00043D6473": "HIT",
        "UPI0002E0D40B": "HIT", "UPI0000167F57": "HIT"}
rc, err, out = run(["--assume-complete"], mock=REAL_TSV, inp=REAL_IN)
check("real output parses (78 rows, 15 columns, no header)", rc == 0 and statuses(out) == REAL, (err, statuses(out)))
check("real output: the only protein that passes is the AntiFam/MobiDB-only one", passed(out) == ["UPI0010EB9177"], passed(out))
rows = table(os.path.join(out, "hits.tsv"))
check("real output: all 78 rows are listed", len(rows) == 78, len(rows))
r = [x for x in rows if x["gene_id"] == "UPI0010EB9177" and x["target"] == "ANF00057"][0]
check("real output: query coverage of a real AntiFam row (96 of 243 = 39.5%)", r["qcov"] == "39.5" and r["evalue"] != "NA", r)
r = [x for x in rows if x["target"] == "PS50270"][0]
check("real output: a ProSiteProfiles row has a score, not an E-value", r["bitscore"] == "55.149357" and r["evalue"] == "NA", r)
mob = [x for x in rows if x["gene_id"] == "UPI0000167F57" and x["source"] == "MobiDBLite"]
pf = [x for x in rows if x["gene_id"] == "UPI0000167F57" and x["source"] == "Pfam"]
check("real output: 22 MobiDB-Lite rows do not count, the one Pfam domain does",
      len(mob) == 22 and not any(x["counts_as_hit"] == "yes" for x in mob) and len(pf) == 1 and pf[0]["counts_as_hit"] == "yes")
no_ipr = [x for x in rows if x["counts_as_hit"] == "yes" and x["note"] == "no InterPro entry"]
check("real output: domain matches without an InterPro accession still count (CDD, SFLD, ProSiteProfiles, PIRSF)",
      {x["source"] for x in no_ipr} >= {"CDD", "SFLD", "ProSiteProfiles", "PIRSF"}, {x["source"] for x in no_ipr})
rc, err, out = run(["--assume-complete", "--ignore-domains", "PF10744"], mock=REAL_TSV, inp=REAL_IN)
check("real output: ignoring the one Pfam domain leaves the MobiDB-Lite protein EXCLUDED_ONLY and passing",
      statuses(out)["UPI0000167F57"] == "EXCLUDED_ONLY" and "UPI0000167F57" in passed(out), statuses(out))
rc, err, out = run(["--assume-complete", "--exclude-analyses", ""], mock=REAL_TSV, inp=REAL_IN)
check("real output: if MobiDB-Lite were counted, the AntiFam protein would wrongly fail",
      statuses(out)["UPI0010EB9177"] == "HIT", statuses(out))
rc, err, out = run([], mock=REAL_TSV, inp=REAL_IN)
check("real output without proof of completion: the AntiFam-only protein is NOT_RUN, real domains still HIT",
      statuses(out)["UPI0010EB9177"] == "NOT_RUN" and statuses(out)["O31533"] == "HIT", statuses(out))

# --- pieces
check("analysis names compare ignoring case and punctuation",
      d.norm("MobiDB-Lite") == d.norm("MobiDBLite") == d.norm("mobidblite") and d.norm("ProSitePatterns") == d.norm("PrositePatterns"))
cmd = d.ips_command("ips.sh", "q.faa", "o.tsv", "/tmp/t", {"applications": "Pfam,CDD"}, 4)
check("the command line always has -dp (no lookup: local and offline)", "-dp" in cmd, cmd)
check("command line: format, output, temp folder, cpu, analyses",
      cmd[cmd.index("-f") + 1] == "tsv" and cmd[cmd.index("-o") + 1] == "o.tsv" and cmd[cmd.index("-T") + 1] == "/tmp/t"
      and cmd[cmd.index("-cpu") + 1] == "4" and cmd[cmd.index("-appl") + 1] == "Pfam,CDD", cmd)
defaults = h.load_defaults("domain-interpro")
check("no setting can switch -dp off", not any("dp" == k or "lookup" in k or "precalc" in k for k in defaults), sorted(defaults))
apps = h.as_list(defaults["applications"])
check("default analyses: the lab's 12 databases plus the four ignored signals",
      len(apps) == 16 and {"Pfam", "CDD", "MobiDBLite", "Coils", "ProSitePatterns", "AntiFam"} <= set(apps), apps)
check("log parsing: analyses and versions",
      d.analyses_from_log("x\nRunning the following analyses:\n[Pfam-38.2,CDD-3.21]\ny") == ["Pfam-38.2", "CDD-3.21"])
check("strip_stops removes '*' and counts sequences",
      (lambda p: (open(p, "w").write(">a\nMK*V\n>b\nMKV\n"), d.strip_stops(p), open(p).read())[1:])(os.path.join(tmp, "s.faa")) is not None)

# --- the run path, with a stand-in interproscan.sh (Linux and WSL only)
if os.name == "nt":
    print("SKIP run-path checks with a stand-in interproscan.sh: needs Linux or WSL")
else:
    stub = os.path.join(tmp, "fake_interproscan.sh")
    open(stub, "w").write(r'''#!/usr/bin/env bash
in=""; out=""; tmp=""
args="$*"
while [[ $# -gt 0 ]]; do
  case "$1" in -i) in="$2"; shift 2;; -o) out="$2"; shift 2;; -T) tmp="$2"; shift 2;;
    -cpu|-f|-appl) shift 2;; -dp) shift;; --version) echo "InterProScan version 5.99-test"; exit 0;; *) shift;; esac
done
echo "$args" >> "$FAKE_CALLS"
cp "$in" "$FAKE_SEEN/$(basename "$in")"
if [[ -n "${FAKE_INVALID:-}" ]]; then
  echo "Invalid input specified for -appl/--applications parameter:"; echo "Analysis Foo does not exist or is deactivated."; exit 1
fi
if grep -q '^>bad' "$in"; then echo "boom"; exit 1; fi
awk -v OFS='\t' '/^>/{split(substr($0,2),a," "); if (a[1] !~ /^nohit/) print a[1],"00000000000000000000000000000000",50,"Pfam","PF00001","Mock",1,40,"1.0E-20","T","26-09-2026","IPR000001","Mock","-","-"}' "$in" > "$out"
if [[ -z "${FAKE_NODONE:-}" ]]; then
  echo "Running the following analyses:"; echo "[Pfam-38.2,CDD-3.21]"; echo "100% done:  InterProScan analyses completed"
fi
exit 0
''')
    os.chmod(stub, 0o755)

    def make_input(name, seqs):
        p = os.path.join(tmp, name)
        os.makedirs(p)
        with open(os.path.join(p, "Sp_final.faa"), "w") as fh:
            for gid, seq in seqs:
                fh.write(">%s\n%s\n" % (gid, seq))
        return p

    def run_real(name, seqs, extra=(), env_extra=None):
        inp = make_input(name + "_in", seqs)
        out = os.path.join(tmp, name + "_out")
        calls, seen = os.path.join(tmp, name + "_calls.txt"), os.path.join(tmp, name + "_seen")
        os.makedirs(seen)
        env = dict(os.environ, FAKE_CALLS=calls, FAKE_SEEN=seen)
        env.update(env_extra or {})
        p = subprocess.run([sys.executable, SCRIPT, "-i", inp, "-o", out, "--iprscan-bin", stub, "--quiet",
                            "--set", "chunk_size=2", "--threads", "3"] + list(extra),
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
        return p, out, calls, seen

    seqs = [("g1", "MKVLAA"), ("nohit2", "MKVL"), ("g3", "MKVLAAGGG"), ("g4", "MK"), ("g5", "MKVLA")]
    p, out, calls, seen = run_real("ok", seqs)
    check("run path succeeds", p.returncode == 0, p.stderr)
    check("run path: genes with matches HIT, the gene with none NO_HIT (its chunk finished)",
          statuses(out) == {"g1": "HIT", "nohit2": "NO_HIT", "g3": "HIT", "g4": "HIT", "g5": "HIT"}, statuses(out))
    call_lines = open(calls).read().splitlines()
    check("run path: 5 genes in chunks of 2 make 3 InterProScan calls", len(call_lines) == 3, call_lines)
    check("run path: every call has -dp", all(" -dp" in " " + c for c in call_lines), call_lines)
    tdirs = [c.split("-T ")[1].split()[0] for c in call_lines]
    check("run path: every call has its own temp folder", len(set(tdirs)) == 3, tdirs)
    check("run path: every call carries the analyses list and the thread count",
          all("-appl Pfam," in c and "-cpu 3" in c and "AntiFam" in c for c in call_lines), call_lines)
    check("run path: the raw results and a 'complete' status are kept",
          os.path.isfile(os.path.join(out, "raw", "Sp.interproscan.tsv"))
          and "hep_status=complete" in open(os.path.join(out, "raw", "Sp.interproscan.status")).read())
    rj = json.load(open(os.path.join(out, "run.json")))
    names = [t["name"] for t in rj["tools"]]
    check("run path: run.json records InterProScan's version and the analyses it ran",
          rj["tools"][0] == {"name": "InterProScan", "version": "5.99-test"} and "Pfam" in names and "CDD" in names
          and rj["analyses_run"] == ["Pfam-38.2", "CDD-3.21"], rj["tools"])
    # a re-run of the raw results through --parse-only needs no proof beyond the status file this module wrote
    rc, err, out2 = run([], mock=os.path.join(out, "raw"), inp=os.path.join(tmp, "ok_in"))
    check("run path: the saved raw results can be re-read with --parse-only", rc == 0 and statuses(out2) == statuses(out), err)

    # Sequences are chunked longest first, so the chunks are [g1 (6), g5 (5)], [g3 (4), bad2 (3)], [g4 (2)].
    # The stand-in fails on the chunk that contains 'bad2', i.e. the middle one.
    p, out, calls, seen = run_real("fail", [("g1", "MKVLAA"), ("bad2", "MKV"), ("g3", "MKVL"), ("g4", "MK"), ("g5", "MKVLA")])
    s = statuses(out)
    check("a chunk that fails: its genes (the chunk holding bad2 and g3) are NOT_RUN; the other chunks keep their answers",
          p.returncode == 2 and s["bad2"] == "NOT_RUN" and s["g3"] == "NOT_RUN"
          and s["g1"] == "HIT" and s["g5"] == "HIT" and s["g4"] == "HIT", (p.stderr, s))
    check("the longest sequences were in the first chunk (longest-first chunking reaches this module)",
          open(os.path.join(seen, "chunk_00000.faa")).read().startswith(">g1"))
    check("...with a warning naming the chunk", "failed" in p.stderr, p.stderr)
    check("...and the status file says incomplete", "hep_status=incomplete" in open(os.path.join(out, "raw", "Sp.interproscan.status")).read())

    p, out, calls, seen = run_real("nodone", seqs, env_extra={"FAKE_NODONE": "1"})
    check("exit code 0 without InterProScan's '100% done' line is NOT trusted: all genes NOT_RUN",
          p.returncode == 2 and set(statuses(out).values()) == {"NOT_RUN"}, (p.stderr, statuses(out)))

    p, out, calls, seen = run_real("invalid", seqs, env_extra={"FAKE_INVALID": "1"})
    check("an analysis name InterProScan rejects stops the run with InterProScan's own message",
          p.returncode != 0 and "rejected the analyses list" in p.stderr and "does not exist or is deactivated" in p.stderr, p.stderr)

    p, out, calls, seen = run_real("stops", [("g1", "MK*VL"), ("g2", "MKV*")])
    seen_text = "".join(open(os.path.join(seen, f)).read() for f in os.listdir(seen))
    check("'*' stop symbols never reach InterProScan, and the user is told",
          "*" not in seen_text and "stop symbols were removed from 2" in p.stderr, (seen_text, p.stderr))

    p = subprocess.run([sys.executable, SCRIPT, "-i", make_input("nobin_in", seqs), "-o", os.path.join(tmp, "nobin_out"),
                        "--iprscan-bin", os.path.join(tmp, "does_not_exist.sh"), "--quiet"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    check("a missing interproscan.sh is reported", p.returncode != 0 and "interproscan.sh not found" in p.stderr, p.stderr)

print()
print("FAILED: %d" % len(fails) if fails else "ALL PASSED")
sys.exit(1 if fails else 0)
