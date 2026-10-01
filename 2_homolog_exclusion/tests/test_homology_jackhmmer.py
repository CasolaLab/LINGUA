"""
test_homology_jackhmmer.py - checks for bin/homology_jackhmmer.py.

Run from anywhere:   python3 tests/test_homology_jackhmmer.py
Needs no HMMER and no databases: it uses mock jackhmmer results in the exact --domtblout format real HMMER 3.4 writes and
(Linux/WSL only) a stand-in jackhmmer for the run path. The expected answers were worked out by hand.

Mock search (example_data/homology-jackhmmer/mock_jackhmmer). Query species Physcomitrium_patens; the database wrongly
also contains that species. Default settings: E-value 1e-3, query coverage >= 50%; a domain adds to coverage only if its
own E-value also passes.
  gene_j1   Ceratodon, one domain, covers 181 of 200 = 90.5%                        -> counts
  gene_j2   hit only to the analysis species itself                                -> IN_PHYLOGENY_ONLY
  gene_j3   hit to itself AND to Blasia (160 of 200 = 80%)                         -> counts (the outside hit)
  gene_j4   covers 52 of 200 = 26%                                                 -> too little of the query -> NO_HIT
  gene_j5   searched, no rows                                                      -> NO_HIT
  gene_j6   not in the searched list                                               -> NOT_RUN
  gene_j7   hit to a target with no species tag                                    -> counts (unknown species is counted)
  gene_j8   hit to Ceratodon_purpureus_R40 (a strain)                              -> counts
  gene_j9   two domains 60 + 51 = 111 of 200 = 55.5%                               -> counts (union)
  gene_j10  whole-sequence E-value 0.05, above 1e-3                                -> NO_HIT
  gene_j11  domain 1 (E 1e-9) covers 40, domain 2 (E 0.5) covers 120: only 40 count -> 20% -> NO_HIT
  gene_j12  covers 90% of the query but 10% of a 1000-residue target               -> counts (no target-coverage cut-off)
  gene_j13  whole-sequence E-value 0.05 (above 1e-3) but one strong domain         -> NO_HIT (the sequence E-value decides)
  gene_j14  two overlapping domains 1-60 and 10-70: merged 70 of 200 = 35%         -> NO_HIT (added up they would be 60.5%)
"""
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
SCRIPT = os.path.join(ROOT, "bin", "homology_jackhmmer.py")
DATA = os.path.join(ROOT, "example_data", "homology-jackhmmer")
INPUT, MOCK = os.path.join(DATA, "input"), os.path.join(DATA, "mock_jackhmmer")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import homology_jackhmmer as j  # noqa: E402
import hep_common as h  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + ((" | " + str(extra)) if (extra and not cond) else ""))
    if not cond:
        fails.append(name)


tmp = tempfile.mkdtemp()
counter = [0]


def run(extra, mock=MOCK, inp=INPUT):
    counter[0] += 1
    out = os.path.join(tmp, "run%d" % counter[0])
    cmd = [sys.executable, SCRIPT, "-i", inp, "-o", out, "--parse-only", mock, "--quiet"] + extra
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
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
    ids = []
    for f in os.listdir(os.path.join(out, "passed")):
        ids += [l[1:].split()[0] for l in open(os.path.join(out, "passed", f)) if l.startswith(">")]
    return sorted(ids)


# --- defaults on the mock results
rc, err, out = run([])
s = statuses(out)
want = {"gene_j1": "HIT", "gene_j2": "IN_PHYLOGENY_ONLY", "gene_j3": "HIT", "gene_j4": "NO_HIT", "gene_j5": "NO_HIT",
        "gene_j6": "NOT_RUN", "gene_j7": "HIT", "gene_j8": "HIT", "gene_j9": "HIT", "gene_j10": "NO_HIT",
        "gene_j11": "NO_HIT", "gene_j12": "HIT", "gene_j13": "NO_HIT", "gene_j14": "NO_HIT"}
check("default settings give the hand-worked status for all 14 genes (exit 2: gene_j6 is NOT_RUN)", rc == 2 and s == want, (err, s))
hits = {(r["gene_id"], r["target"]): r for r in table(os.path.join(out, "hits.tsv"))}
check("query coverage: 90.5, 80, 26, 55.5 (union), 20 (weak domain not counted), 35 (overlap merged, not added)",
      [hits[k]["qcov"] for k in (("gene_j1", "c1__Ceratodon_purpureus"), ("gene_j3", "b1__Blasia_pusilla"),
                                 ("gene_j4", "c2__Ceratodon_purpureus"), ("gene_j9", "f1__Fossombronia_cristula"),
                                 ("gene_j11", "f3__Fossombronia_cristula"),
                                 ("gene_j14", "f6__Fossombronia_cristula"))] == ["90.5", "80.0", "26.0", "55.5", "20.0", "35.0"],
      {k: v["qcov"] for k, v in hits.items()})
check("target coverage is reported: 10.0 for the long target of gene_j12", hits[("gene_j12", "f4__Fossombronia_cristula")]["tcov"] == "10.0")
check("the target species comes from the ID tag; an untagged target is NA and is counted with a note",
      hits[("gene_j1", "c1__Ceratodon_purpureus")]["target_species"] == "Ceratodon_purpureus"
      and hits[("gene_j7", "u1")]["target_species"] == "NA" and "unknown" in hits[("gene_j7", "u1")]["note"])
check("no table has an empty cell",
      all(v != "" for f in ("hits.tsv", "gene_status.tsv", "summary.tsv") for r in table(os.path.join(out, f)) for v in r.values()))
check("passed genes are exactly the ones that pass: j2, j4, j5, j10, j11, j13, j14 (NOT_RUN j6 never passes)",
      passed(out) == ["gene_j10", "gene_j11", "gene_j13", "gene_j14", "gene_j2", "gene_j4", "gene_j5"], passed(out))
rj = json.load(open(os.path.join(out, "run.json")))
check("run.json records that results were parsed and the analysis species", rj["analysis_species"] == ["Physcomitrium_patens"]
      and "parsed" in rj["tools"][0]["name"], rj.get("tools"))

# --- settings
rc, err, out = run(["--set", "min_qcov=0"])
check("min_qcov=0: low-coverage hits count (j4, j11, j14), but a failing sequence E-value still excludes (j13)",
      statuses(out)["gene_j4"] == "HIT" and statuses(out)["gene_j11"] == "HIT" and statuses(out)["gene_j14"] == "HIT"
      and statuses(out)["gene_j13"] == "NO_HIT", err)
rc, err, out = run(["--set", "min_tcov=50"])
check("min_tcov=50: gene_j12 (10% of its target) no longer counts", statuses(out)["gene_j12"] == "NO_HIT" and statuses(out)["gene_j1"] == "HIT", err)
rc, err, out = run(["--set", "ignore_in_phylogeny=no"])
check("ignore_in_phylogeny=no: the analysis species' own hit counts (j2)", statuses(out)["gene_j2"] == "HIT", err)
rc, err, out = run(["--set", "evalue=1e-40"])
check("a stricter E-value drops hits (j1 has 1e-35); j2 (whole-sequence 1e-40 but domain 1e-38) has no passing domain, so no coverage",
      statuses(out)["gene_j1"] == "NO_HIT" and statuses(out)["gene_j2"] == "NO_HIT", err)
rc, err, out = run(["--set", "ignore_taxa=Fossombronia"])
s = statuses(out)
check("ignore_taxa with one word is a genus: Fossombronia hits (j9, j12) are exempt", s["gene_j9"] == "IN_PHYLOGENY_ONLY"
      and s["gene_j12"] == "IN_PHYLOGENY_ONLY" and s["gene_j1"] == "HIT", (err, s))
rc, err, out = run(["--set", "ignore_taxa=Ceratodon_purpureus"])
s = statuses(out)
check("ignore_taxa at species level: a strain (R40) is the same species, exempt; j1 too", s["gene_j8"] == "IN_PHYLOGENY_ONLY"
      and s["gene_j1"] == "IN_PHYLOGENY_ONLY" and s["gene_j9"] == "HIT", (err, s))
rc, err, out = run(["--set", "ignore_taxa=Ceratodon_purpureus", "--set", "match_level=exact"])
check("match_level=exact: the strain is a different name, so it counts (j8)", statuses(out)["gene_j8"] == "HIT"
      and statuses(out)["gene_j1"] == "IN_PHYLOGENY_ONLY", err)
rc, err, out = run(["--set", "iterations=0"])
check("iterations=0 is refused", rc != 0 and "iterations" in err, err)
rc, err, out = run(["--set", "species_source=map"])
check("species_source=map without a taxon map is refused", rc != 0 and "taxon_map" in err, err)

# --- completeness
rc, err, out = run([], mock=os.path.join(MOCK, "truncated.domtbl"))
check("a result file without '# [ok]' is NOT trusted: every gene is NOT_RUN, even those with hits",
      rc == 2 and set(statuses(out).values()) == {"NOT_RUN"} and "# [ok]" in err, (err, statuses(out)))
rc, err, out = run([], mock=os.path.join(tmp, "nowhere"))
check("a missing result file gives NOT_RUN, never NO_HIT", set(statuses(out).values()) == {"NOT_RUN"}, (err, statuses(out)))
# a finished file made by hand (no record of what was searched) covers the whole input
hand = os.path.join(tmp, "hand.domtbl")
open(hand, "w").write(open(os.path.join(MOCK, "Physcomitrium_patens.jackhmmer.domtbl")).read())
rc, err, out = run([], mock=hand)
check("a finished file with no searched-list covers every input gene: gene_j6 becomes NO_HIT, the rest as before",
      statuses(out)["gene_j6"] == "NO_HIT" and statuses(out)["gene_j1"] == "HIT", (err, statuses(out)))

ok, rows = j.parse_domtbl(os.path.join(MOCK, "Physcomitrium_patens.jackhmmer.domtbl"))
check("parse_domtbl: finished file, 16 domain rows, coordinates read as query (hmm) and target (ali)",
      ok and len(rows) == 16 and rows[0]["qid"] == "gene_j1" and (rows[0]["qstart"], rows[0]["qend"]) == (1, 181)
      and (rows[0]["sstart"], rows[0]["send"]) == (5, 185) and rows[0]["qlen"] == 200 and rows[0]["tlen"] == 300, rows[:1])
check("parse_domtbl: truncated file is not complete", j.parse_domtbl(os.path.join(MOCK, "truncated.domtbl"))[0] is False)
bad = os.path.join(tmp, "bad.domtbl")
open(bad, "w").write("# header\nonly three columns\n# [ok]\n")
try:
    j.parse_domtbl(bad)
    check("a malformed domain line is an error, not a silent skip", False)
except h.HepError:
    check("a malformed domain line is an error, not a silent skip", True)

# --- real HMMER 3.4 output (example_data/homology-jackhmmer/real_sample: 3 real proteins and 2 shuffled ones searched with
# jackhmmer against 7 moss/liverwort proteomes; the file is the real --domtblout with at most 8 rows kept per query)
REAL = os.path.join(DATA, "real_sample")
rc, err, out = run([], mock=os.path.join(REAL, "jackhmmer"), inp=os.path.join(REAL, "input"))
check("real HMMER output: two proteins with hits are HIT; the finished file proves the other three, which have no rows, were searched",
      rc == 0 and statuses(out) == {"XP_024396482.1": "HIT", "XP_024397620.1": "HIT", "XP_073386763.1": "NO_HIT",
                                   "SHUF_NP_001365144.1": "NO_HIT", "SHUF_NP_001365145.1": "NO_HIT"}, (err, statuses(out)))
okr, real_rows = j.parse_domtbl(os.path.join(REAL, "jackhmmer", "Physcomitrium_patens.jackhmmer.domtbl"))
check("real HMMER output: read as finished, 16 domain rows, tagged targets", okr and len(real_rows) == 16
      and all(h.split_id(r["tid"])[1] for r in real_rows), len(real_rows))

# --- pieces
cmd = j.jackhmmer_command("q.faa", "db.faa", "o.domtbl", {"iterations": 3, "evalue": 0.001, "incE": 0.001, "dbsize": None})
check("jackhmmer command: iterations, E-values, no alignments, one CPU, big -o output discarded, table output, query then database",
      cmd[:1] == ["jackhmmer"] and cmd[cmd.index("-N") + 1] == "3" and cmd[cmd.index("-E") + 1] == "0.001"
      and cmd[cmd.index("--incE") + 1] == "0.001" and "--noali" in cmd and cmd[cmd.index("--cpu") + 1] == "1"
      and cmd[cmd.index("-o") + 1] == os.devnull and cmd[cmd.index("--domtblout") + 1] == "o.domtbl"
      and cmd[-2:] == ["q.faa", "db.faa"] and "-Z" not in cmd, cmd)
cmd = j.jackhmmer_command("q.faa", "db.faa", "o.domtbl", {"iterations": 3, "evalue": 0.001, "incE": 0.001, "dbsize": 250000})
check("-Z is passed when dbsize is set", cmd[cmd.index("-Z") + 1] == "250000", cmd)
cmd = j.jackhmmer_command("q.faa", "db.faa", "o.domtbl", {"iterations": 3, "evalue": 0.001, "incE": 0.001, "dbsize": None}, threads=3)
check("--threads is honoured like in the other modules: it becomes --cpu for each jackhmmer job", cmd[cmd.index("--cpu") + 1] == "3", cmd)
check("the network is never involved: no remote option exists in the defaults",
      not any("remote" in k for k in h.load_defaults("homology-jackhmmer")))
fa = os.path.join(tmp, "inv.faa")
open(fa, "w").write(">a__Sp_one x\nMK\n>b__Sp_one\nMK\n>c__Sp_two\nMK\n>d\nMK\n")
inv = j.fasta_inventory(fa)
check("FASTA inventory counts IDs and species tags", inv["ids"] == 4 and inv["tagged"] == 3
      and inv["species"] == {"Sp_one": 2, "Sp_two": 1}, inv)
for name, path, text in (("a BLAST prefix (not a file)", os.path.join(tmp, "prefix"), None),
                         ("a file that is not FASTA", os.path.join(tmp, "notfasta.txt"), "hello\n")):
    if text is not None:
        open(path, "w").write(text)
    try:
        j.inspect_db(path)
        check("database check refuses " + name, False)
    except j.DbProblem:
        check("database check refuses " + name, True)
check("a FASTA file is accepted as the database", j.inspect_db(fa)["type"] == "fasta")

p = subprocess.run([sys.executable, SCRIPT, "-i", INPUT, "-o", os.path.join(tmp, "nodb")], stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, universal_newlines=True)
check("--db is required unless --parse-only", p.returncode != 0 and "--db is required" in p.stderr, p.stderr)

# --- the run path, with a stand-in jackhmmer (Linux and WSL only)
if os.name == "nt":
    print("SKIP run-path checks with a stand-in jackhmmer: need Linux or WSL")
else:
    bindir = os.path.join(tmp, "fakebin")
    os.makedirs(bindir)
    stub = os.path.join(bindir, "jackhmmer")
    open(stub, "w").write(r'''#!/usr/bin/env bash
if [[ "$1" == "-h" ]]; then echo "# jackhmmer :: iteratively search a protein sequence against a protein database"; echo "# HMMER 9.9 (Jan 2099); http://hmmer.org/"; exit 0; fi
args="$*"; tbl=""; pos=()
while [[ $# -gt 0 ]]; do case "$1" in --domtblout) tbl="$2"; shift 2;; -N|-E|--incE|--cpu|-o|-Z) shift 2;; --noali) shift;; *) pos+=("$1"); shift;; esac; done
q="${pos[0]}"; db="${pos[1]}"
echo "$args" >> "$FAKE_CALLS"
cp "$q" "$FAKE_SEEN/$(basename "$q")"
if grep -q '^>bad' "$q"; then echo "boom" >&2; exit 1; fi
dbbase="$(basename "$db")"; var="HITS_${dbbase//./_}"; mode=all; list=""
if [[ -n "${!var+x}" ]]; then mode=list; list="${!var}"; fi
{
echo "# target name  accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias from to from to from to acc description of target"
awk -v MODE="$mode" -v L=" $list " '
function flush() { if (id == "") return;
  if ((MODE == "all" && id !~ /^nohit/) || (MODE == "list" && index(L, " " id " ") > 0))
    printf "T__Ceratodon_purpureus - %d %s - %d 1e-30 100.0 1.0 1 1 1e-30 1e-30 90.0 0.5 1 %d 1 %d 1 %d 0.9 -\n", len, id, len, len, len, len }
/^>/ { flush(); id = substr($1, 2); len = 0; next } { len += length($0) }
END { flush() }' "$q"
if [[ -z "${FAKE_NOFOOTER:-}" ]]; then echo "#"; echo "# Program:         jackhmmer"; echo "# [ok]"; fi
} > "$tbl"
exit 0
''')
    os.chmod(stub, 0o755)

    def make_input(name, seqs):
        p_ = os.path.join(tmp, name)
        os.makedirs(p_)
        with open(os.path.join(p_, "Sp_final.faa"), "w") as fh:
            for gid, seq in seqs:
                fh.write(">%s\n%s\n" % (gid, seq))
        return p_

    def make_db(name, ids=None):
        path = os.path.join(tmp, name)
        with open(path, "w") as fh:
            for i in (ids or ["s%d__Ceratodon_purpureus" % k for k in range(1, 6)]):
                fh.write(">%s\nMKVLAAGG\n" % i)
        return path

    def run_real(name, seqs, dbs=("db1.faa",), extra=(), env_extra=None):
        inp = make_input(name + "_in", seqs)
        out = os.path.join(tmp, name + "_out")
        calls, seen = os.path.join(tmp, name + "_calls.txt"), os.path.join(tmp, name + "_seen")
        os.makedirs(seen)
        env = dict(os.environ, PATH=bindir + os.pathsep + os.environ["PATH"], FAKE_CALLS=calls, FAKE_SEEN=seen)
        env.update(env_extra or {})
        cmd = [sys.executable, SCRIPT, "-i", inp, "-o", out, "--quiet", "--set", "chunk_size=2"]
        for d in dbs:
            cmd += ["--db", d if os.path.isabs(d) else make_db(name + "_" + d)]
        p = subprocess.run(cmd + list(extra), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
        return p, out, calls, seen

    seqs = [("g1", "MKVLAA"), ("nohit2", "MKVL"), ("g3", "MKVLAAGGG"), ("g4", "MK"), ("g5", "MKVLA")]
    p, out, calls, seen = run_real("ok", seqs)
    check("run path succeeds", p.returncode == 0, p.stderr)
    check("run path: genes with hits HIT, the gene with no rows NO_HIT (the finished file proves it was searched)",
          statuses(out) == {"g1": "HIT", "nohit2": "NO_HIT", "g3": "HIT", "g4": "HIT", "g5": "HIT"}, statuses(out))
    call_lines = open(calls).read().splitlines()
    check("run path: 5 genes in chunks of 2 make 3 jackhmmer calls", len(call_lines) == 3, call_lines)
    check("run path: every call has the iterations, E-values, --noali, one CPU and the table output",
          all("-N 3" in c and "-E 0.001" in c and "--incE 0.001" in c and "--noali" in c and "--cpu 1" in c
              and "--domtblout" in c for c in call_lines), call_lines)
    rj = json.load(open(os.path.join(out, "run.json")))
    check("run path: run.json records HMMER's version, the database and its species inventory",
          rj["tools"][0] == {"name": "jackhmmer (HMMER)", "version": "HMMER 9.9"} and rj["database_info"]["type"] == "fasta"
          and rj["database_inventory"]["ids"] == 5 and rj["database_inventory"]["species"] == {"Ceratodon_purpureus": 5}, rj)
    check("run path: the raw results and the searched-list are kept",
          os.path.isfile(os.path.join(out, "raw", "Sp.jackhmmer.domtbl")) and os.path.isfile(os.path.join(out, "raw", "Sp.jackhmmer.searched")))
    rc, err, out2 = run([], mock=os.path.join(out, "raw"), inp=os.path.join(tmp, "ok_in"))
    check("run path: the saved raw results can be re-read with --parse-only", rc == 0 and statuses(out2) == statuses(out), err)

    # chunks are longest first: [g1 (6), g5 (5)], [g3 (4), bad2 (3)], [g4 (2)]
    p, out, calls, seen = run_real("fail", [("g1", "MKVLAA"), ("bad2", "MKV"), ("g3", "MKVL"), ("g4", "MK"), ("g5", "MKVLA")])
    s = statuses(out)
    check("a chunk that fails: its genes are NOT_RUN; the other chunks keep their answers",
          p.returncode == 2 and s["bad2"] == "NOT_RUN" and s["g3"] == "NOT_RUN" and s["g1"] == "HIT" and s["g5"] == "HIT"
          and s["g4"] == "HIT", (p.stderr, s))
    searched = open(os.path.join(out, "raw", "Sp.jackhmmer.searched")).read().split()
    check("a failed chunk's genes are NOT in the searched-list, so re-reading the raw files keeps them NOT_RUN",
          "bad2" not in searched and "g3" not in searched and "g1" in searched, searched)
    rc, err, out2 = run([], mock=os.path.join(out, "raw"), inp=os.path.join(tmp, "fail_in"))
    check("re-reading a partly failed run gives the same statuses", statuses(out2) == s, (err, statuses(out2)))
    p, out, calls, seen = run_real("nofoot", seqs, env_extra={"FAKE_NOFOOTER": "1"})
    check("output without HMMER's '# [ok]' is NOT trusted, even with exit code 0: all genes NOT_RUN",
          p.returncode == 2 and set(statuses(out).values()) == {"NOT_RUN"}, (p.stderr, statuses(out)))
    p, out, calls, seen = run_real("nobin", seqs, env_extra={"PATH": "/usr/bin:/bin"})
    check("jackhmmer missing from PATH is reported", p.returncode != 0 and "jackhmmer not found" in p.stderr, p.stderr)
    p, out, calls, seen = run_real("blastdb", seqs, dbs=(os.path.join(tmp, "some_blast_prefix"),))
    check("a BLAST database prefix is refused with a clear message, and nothing is left behind",
          p.returncode != 0 and "not a file" in p.stderr and not os.path.exists(out), p.stderr)
    dbself = make_db("self.faa", ["a__Sp", "b__Sp", "c__Ceratodon_purpureus"])
    ana = os.path.join(tmp, "selfin_in")
    os.makedirs(ana)
    open(os.path.join(ana, "Ceratodon_purpureus_final.faa"), "w").write(">g1\nMKVLAA\n")
    outs = os.path.join(tmp, "selfin_out")
    env = dict(os.environ, PATH=bindir + os.pathsep + os.environ["PATH"], FAKE_CALLS=os.path.join(tmp, "selfin_calls"),
              FAKE_SEEN=tmp)
    p = subprocess.run([sys.executable, SCRIPT, "-i", ana, "-o", outs, "--db", dbself, "--quiet"], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True, env=env)
    check("an analysis species INSIDE the database is named in a warning, and its hits are not counted",
          p.returncode == 0 and "INSIDE this database" in p.stderr and "Ceratodon_purpureus" in p.stderr
          and statuses(outs) == {"g1": "IN_PHYLOGENY_ONLY"}, (p.stderr, statuses(outs)))

    # --- several databases
    seq5 = [("g1", "MKVLAA"), ("g2", "MKVL"), ("g3", "MKVLAAGGG"), ("g4", "MK"), ("g5", "MKVLA")]
    hits_env = {"HITS_" + "d1_faa": "g1 g2", "HITS_" + "d2_faa": "g3", "HITS_" + "d3_faa": "g4"}
    dbs3 = [os.path.join(tmp, "multi", n) for n in ("d1.faa", "d2.faa", "d3.faa")]
    os.makedirs(os.path.join(tmp, "multi"))
    for d in dbs3:
        make_db(os.path.join("multi", os.path.basename(d)))
    p, out, calls, seen = run_real("multi", seq5, dbs=dbs3, env_extra=hits_env)
    check("WITHOUT --sequential every database searches every gene (three labelled runs)",
          p.returncode == 0 and all(os.path.isdir(os.path.join(out, l)) for l in ("d1", "d2", "d3"))
          and len(open(calls).read().splitlines()) == 9 and statuses(os.path.join(out, "d1"))["g1"] == "HIT"
          and statuses(os.path.join(out, "d1"))["g3"] == "NO_HIT", p.stderr)
    p, out, calls, seen = run_real("seq", seq5, dbs=dbs3, extra=["--sequential"], env_extra=hits_env)
    check("--sequential succeeds", p.returncode == 0, p.stderr)
    check("--sequential: the second database searches only the 3 survivors of the first, the third only the 2 that remain",
          sorted(l for l in os.listdir(out) if os.path.isdir(os.path.join(out, l))) == ["d1", "d2", "d3"]
          and passed(os.path.join(out, "d1")) == ["g3", "g4", "g5"] and passed(os.path.join(out, "d2")) == ["g4", "g5"]
          and passed(os.path.join(out, "d3")) == ["g5"], (p.stderr, [passed(os.path.join(out, x)) for x in ("d1", "d2", "d3")]))
    p, out, calls, seen = run_real("seqfail", seq5, dbs=[dbs3[0], os.path.join(tmp, "missing.faa"), dbs3[2]],
                                   extra=["--sequential"], env_extra=hits_env)
    check("--sequential stops when a database cannot be used, and says which were not searched",
          p.returncode != 0 and "NOT searched" in p.stderr and "d3" in p.stderr and not os.path.isdir(os.path.join(out, "d3")), p.stderr)

print()
print("FAILED: %d" % len(fails) if fails else "ALL PASSED")
sys.exit(1 if fails else 0)
