"""
test_homology_blast.py - checks for bin/homology_blast.py.

Run from anywhere:   python3 tests/test_homology_blast.py
Needs no BLAST+ and no databases: it uses mock BLAST results in the exact format real BLAST+ 2.17.0 writes and (Linux/WSL
only) stand-in blastp/tblastn/blastdbcmd programs for the run path. The expected answers were worked out by hand.

Mock protein search (example_data/homology-blast/mock_blast). Query species Physcomitrium_patens; the database wrongly
also contains that species (targets tagged __Physcomitrium_patens). Default settings: E-value 1e-3, query coverage >= 50%.
  gene_c1   one hit, Ceratodon, covers 90.5% of the query                          -> counts
  gene_c2   hit only to the analysis species itself                                -> IN_PHYLOGENY_ONLY
  gene_c3   hit to itself AND to Blasia (covers 80%)                               -> counts (the outside hit)
  gene_c4   hit covers only 26% of the query                                       -> too little of the query
  gene_c5   searched, no hits                                                      -> NO_HIT
  gene_c6   missing from the result file                                           -> NOT_RUN
  gene_c7   hit to a target with no species tag                                    -> counts (unknown species is counted)
  gene_c8   hit to Ceratodon_purpureus_R40 (a strain)                              -> counts
  gene_c9   hit to Marchantia_polymorpha_subsp._ruderalis                          -> counts (ignore_taxa can exempt it)
  gene_c10  query ID carries __Physcomitrium_patens; hit to Blasia                 -> counts
  gene_c11  E-value 0.05, above 1e-3                                               -> NO_HIT
  gene_c12  two alignments of 60 and 51 residues: 111 of 200 together = 55.5%      -> counts (union)
  gene_c13  identity 25%                                                           -> counts (no identity cut-off by default)
  gene_c14  covers 90% of the query but 10% of a 1000-residue target               -> counts (no target-coverage cut-off)
"""
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
SCRIPT = os.path.join(ROOT, "bin", "homology_blast.py")
DATA = os.path.join(ROOT, "example_data", "homology-blast")
INPUT, MOCK = os.path.join(DATA, "input"), os.path.join(DATA, "mock_blast")
ORG_IN, ORG_MOCK = os.path.join(DATA, "input_organelle"), os.path.join(DATA, "mock_tblastn")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import homology_blast as b  # noqa: E402
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


BASE = {"gene_c1": "HIT", "gene_c2": "IN_PHYLOGENY_ONLY", "gene_c3": "HIT", "gene_c4": "NO_HIT", "gene_c5": "NO_HIT",
        "gene_c6": "NOT_RUN", "gene_c7": "HIT", "gene_c8": "HIT", "gene_c9": "HIT", "gene_c10": "HIT",
        "gene_c11": "NO_HIT", "gene_c12": "HIT", "gene_c13": "HIT", "gene_c14": "HIT"}

# --- default rules on the mock data
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
check("genes that pass: no counted hit (the hit-to-itself gene passes: it is IN_PHYLOGENY_ONLY)",
      passed(out) == ["gene_c11", "gene_c2", "gene_c4", "gene_c5"], passed(out))   # text order: c11 sorts before c2
hits = {}
for r in table(os.path.join(out, "hits.tsv")):
    hits.setdefault(r["gene_id"], []).append(r)
r = hits["gene_c1"][0]
check("a counted hit reports target, species, E-value, bit score, identity and both coverages",
      r["target"] == "T1__Ceratodon_purpureus" and r["target_species"] == "Ceratodon_purpureus" and r["evalue"] == "1e-30"
      and r["bitscore"] == "120.0" and r["pident"] == "60.0" and r["qcov"] == "90.5" and r["tcov"] == "95.3"
      and r["counts_as_hit"] == "yes" and r["source"] == "homology-blast", r)
check("a hit to the analysis species is listed but not counted, with the reason",
      hits["gene_c2"][0]["counts_as_hit"] == "no" and hits["gene_c2"][0]["note"] == "in-phylogeny hit (Physcomitrium_patens), not counted",
      hits["gene_c2"])
check("an outside hit counts even when a hit to the analysis species is also present",
      sorted((x["target"], x["counts_as_hit"]) for x in hits["gene_c3"]) == [("S2__Physcomitrium_patens", "no"), ("T2__Blasia_pusilla", "yes")],
      hits["gene_c3"])
check("a hit that covers too little of the query says so", hits["gene_c4"][0]["note"] == "query coverage 26.0% below 50%", hits["gene_c4"])
check("a hit to a target with no species tag counts, and says the species is unknown",
      hits["gene_c7"][0]["counts_as_hit"] == "yes" and hits["gene_c7"][0]["target_species"] == "NA"
      and hits["gene_c7"][0]["note"] == "target species unknown, counted", hits["gene_c7"])
check("a hit above the E-value threshold is listed with the reason", hits["gene_c11"][0]["note"] == "E-value above 0.001", hits["gene_c11"])
check("two alignments to one target are merged; coverage is the union (111 of 200 = 55.5%)",
      len(hits["gene_c12"]) == 1 and hits["gene_c12"][0]["qcov"] == "55.5" and hits["gene_c12"][0]["tcov"] == "55.5"
      and float(hits["gene_c12"][0]["evalue"]) == 1e-8, hits["gene_c12"])
check("a query ID with a __Species suffix is matched by its canonical ID", "gene_c10" in statuses(out))
sm = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(out, "summary.tsv"))]
check("summary.tsv (14 genes: 9 hit, 3 no hit, 1 in-phylogeny only, 1 not run; 4 passed)",
      sm[1] == ["Physcomitrium_patens", "14", "4", "9", "3", "0", "0", "1", "1"], sm)
rj = json.load(open(os.path.join(out, "run.json")))
check("run.json records the analysis species, the match level and the in-phylogeny rule",
      rj["analysis_species"] == ["Physcomitrium_patens"] and rj["match_level"] == "species"
      and rj["ignore_in_phylogeny"] is True and rj["program"] == "blastp", rj)

# --- settings change the answers as intended
rc, err, out = run(["--ignore-in-phylogeny", "no"])
check("ignore_in_phylogeny=no: a hit to the analysis species now counts", statuses(out) == dict(BASE, gene_c2="HIT"), statuses(out))
rc, err, out = run(["--ignore-taxa", "Marchantia polymorpha"])
check("ignore_taxa: a subspecies is treated as the species (species level)",
      statuses(out) == dict(BASE, gene_c9="IN_PHYLOGENY_ONLY"), statuses(out))
rc, err, out = run(["--ignore-taxa", "Marchantia polymorpha", "--match-level", "exact"])
check("match_level=exact: the subspecies is NOT the same, so its hit counts again",
      statuses(out) == BASE and statuses(out)["gene_c2"] == "IN_PHYLOGENY_ONLY", statuses(out))
rc, err, out = run(["--ignore-taxa", "Ceratodon"])
s = statuses(out)
check("ignore_taxa with a genus: every Ceratodon hit is exempt, strains included",
      s["gene_c1"] == "IN_PHYLOGENY_ONLY" and s["gene_c8"] == "IN_PHYLOGENY_ONLY" and s["gene_c3"] == "HIT", s)
rc, err, out = run(["--match-level", "genus"])
check("match_level=genus still recognises the analysis species", statuses(out) == BASE, statuses(out))
rc, err, out = run(["--min-qcov", "0"])
check("min_qcov=0: the 26% hit now counts", statuses(out) == dict(BASE, gene_c4="HIT"), statuses(out))
rc, err, out = run(["--min-pident", "40"])
check("min_pident=40: only the 25% identity hit stops counting", statuses(out) == dict(BASE, gene_c13="NO_HIT"), statuses(out))
rc, err, out = run(["--min-tcov", "50"])
check("min_tcov=50: only the hit covering 10% of its target stops counting", statuses(out) == dict(BASE, gene_c14="NO_HIT"), statuses(out))
rc, err, out = run(["--evalue", "1e-12"])
check("evalue=1e-12 is applied again when parsing",
      statuses(out) == dict(BASE, gene_c3="IN_PHYLOGENY_ONLY", gene_c10="NO_HIT", gene_c12="NO_HIT", gene_c13="NO_HIT"), statuses(out))
rj = json.load(open(os.path.join(out, "run.json")))
check("run.json records the changed setting", rj["non_default_settings"]["evalue"]["value"] == 1e-12, rj["non_default_settings"])

# --- species from a taxon map instead of the ID tag
tmap = os.path.join(MOCK, "taxon_map.tsv")
rc, err, out = run(["--species-source", "map", "--taxon-map", tmap])
s = statuses(out)
check("species_source=map: species come from the map (with or without the __Species part of the ID)",
      s["gene_c2"] == "IN_PHYLOGENY_ONLY" and s["gene_c7"] == "IN_PHYLOGENY_ONLY", s)
check("species_source=map: hits whose IDs are not in the map are counted, never exempted",
      s["gene_c3"] == "HIT" and s["gene_c1"] == "HIT", s)
rc, err, out = run(["--species-source", "map"])
check("species_source=map without a taxon map is refused", rc != 0 and "taxon_map" in err, err)

# --- organelle search (tblastn): the analysis species' own organelle genome must count
rc, err, out = run(["--program", "tblastn"], mock=ORG_MOCK, inp=ORG_IN)
s = statuses(out)
oh = {r["gene_id"]: r for r in table(os.path.join(out, "hits.tsv"))}
check("tblastn with defaults: strong hit counts, weak hit covers too little, own-species organelle is exempt",
      s == {"gene_o1": "HIT", "gene_o2": "NO_HIT", "gene_o3": "IN_PHYLOGENY_ONLY"}, s)
check("tblastn does not report target coverage (the target coordinates are nucleotides)", oh["gene_o1"]["tcov"] == "NA", oh["gene_o1"])
preset = os.path.join(ROOT, "presets", "organelle_tblastn.conf")
rc, err, out = run(["--config", preset], mock=ORG_MOCK, inp=ORG_IN)
check("the organelle preset: any hit counts, including the analysis species' own organelle genome",
      statuses(out) == {"gene_o1": "HIT", "gene_o2": "HIT", "gene_o3": "HIT"}, (err, statuses(out)))
rc, err, out = run(["--program", "tblastn", "--min-tcov", "50"], mock=ORG_MOCK, inp=ORG_IN)
check("min_tcov with tblastn is refused, with the reason", rc != 0 and "nucleotides" in err, err)
rc, err, out = run(["--program", "tblastn", "--ungapped", "yes"], mock=ORG_MOCK, inp=ORG_IN)
check("ungapped without comp_based_stats=0 is refused before running", rc != 0 and "comp_based_stats=0" in err, err)
rc, err, out = run(["--program", "megablast"])
check("an unknown program is refused", rc != 0 and "blastp or tblastn" in err, err)

# --- REAL BLAST+ 2.17.0 output, kept as small regression tests (real proteins, real lines). Answers worked out by hand:
# blastp against a database that (wrongly) contains the query species:
#   XP_073393566.1  one hit, to itself (Physcomitrium_patens), 100%                    -> IN_PHYLOGENY_ONLY
#   XP_024403912.1  one hit, to itself                                                 -> IN_PHYLOGENY_ONLY
#   YP_539006.1     hits to itself AND two to Sphagnum_fallax (84% over 175 aa; 55%)   -> HIT (the outside hits count)
REAL = os.path.join(DATA, "real_sample")
rc, err, out = run([], mock=os.path.join(REAL, "blast"), inp=os.path.join(REAL, "input"))
check("real blastp output: succeeds and no gene is NOT_RUN", rc == 0 and "NOT_RUN" not in statuses(out).values(), err)
check("real blastp output: statuses",
      statuses(out) == {"XP_073393566.1": "IN_PHYLOGENY_ONLY", "XP_024403912.1": "IN_PHYLOGENY_ONLY", "YP_539006.1": "HIT"},
      statuses(out))
rh = table(os.path.join(out, "hits.tsv"))
counted = sorted((r["target"], r["target_species"]) for r in rh if r["counts_as_hit"] == "yes")
check("real blastp output: the two Sphagnum hits count and the self hit does not",
      counted == [("KAH8932845.1__Sphagnum_fallax", "Sphagnum_fallax"), ("KAH8952199.1__Sphagnum_fallax", "Sphagnum_fallax")], counted)
check("real blastp output: 5 alignment rows become 5 (gene, target) hits with real coverage",
      len(rh) == 5 and {r["target"]: r["qcov"] for r in rh}["KAH8952199.1__Sphagnum_fallax"] == "57.1", [(r["target"], r["qcov"]) for r in rh])
rc, err, out = run(["--ignore-in-phylogeny", "no"], mock=os.path.join(REAL, "blast"), inp=os.path.join(REAL, "input"))
check("real blastp output: with the exemption off, a gene whose only hit is itself now fails",
      set(statuses(out).values()) == {"HIT"}, statuses(out))
# tblastn against the real Marchantia chloroplast genome (a real organelle genome):
#   Plastid_encoded: NP_904194.1 (93% over all 475 aa, E 0) and NP_904193.1 (85% over all 312 aa) -> HIT; NP_904188.1 no hit -> NO_HIT
#   Physcomitrium_patens: YP_539019.1 has five alignments to the genome, best E 9.13e-52. Their query ranges 264-354,
#   71-190, 213-246, 119-201 and 292-350 overlap; together they cover 256 of its 489 residues (52.4%), just above the
#   default 50% rule, so the gene has a counted hit.
RT = os.path.join(DATA, "real_sample_tblastn")
rc, err, out = run(["--program", "tblastn"], mock=os.path.join(RT, "blast"), inp=os.path.join(RT, "input"))
s = {(r["species"], r["gene_id"]): r["status"] for r in table(os.path.join(out, "gene_status.tsv"))}
check("real tblastn output: strong plastid hits count, the protein with no hit passes",
      s[("Plastid_encoded", "NP_904194.1")] == "HIT" and s[("Plastid_encoded", "NP_904193.1")] == "HIT"
      and s[("Plastid_encoded", "NP_904188.1")] == "NO_HIT", s)
oh = [r for r in table(os.path.join(out, "hits.tsv")) if r["gene_id"] == "YP_539019.1"]
check("real tblastn output: several alignments of one protein to one genome merge into one hit; no target coverage is reported",
      len(oh) == 1 and oh[0]["tcov"] == "NA" and oh[0]["evalue"] == "9.13e-52", oh)
check("real tblastn output: the five overlapping alignments cover 256 of 489 residues (52.4%), so the default 50% rule is met",
      oh[0]["qcov"] == "52.4" and s[("Physcomitrium_patens", "YP_539019.1")] == "HIT", (oh, s))
rc, err, out = run(["--config", os.path.join(ROOT, "presets", "organelle_tblastn.conf")], mock=os.path.join(RT, "blast"), inp=os.path.join(RT, "input"))
s = {(r["species"], r["gene_id"]): r["status"] for r in table(os.path.join(out, "gene_status.tsv"))}
check("real tblastn output with the organelle preset: the protein YP_539019.1 has a hit and counts",
      s[("Physcomitrium_patens", "YP_539019.1")] == "HIT" and s[("Plastid_encoded", "NP_904188.1")] == "NO_HIT", s)

# --- a result that never finished must never read as "no hit"
mock_lines = open(os.path.join(MOCK, "Physcomitrium_patens.blast.tsv")).read().splitlines()


def joined(name, lines):
    d = os.path.join(tmp, name)
    os.makedirs(d)
    open(os.path.join(d, "Physcomitrium_patens.blast.tsv"), "w").write("\n".join(lines) + "\n")
    return d


rc, err, out = run([], mock=joined("nofooter", mock_lines[:-1]))
check("result without BLAST's completion footer: every gene is NOT_RUN and a warning says so",
      set(statuses(out).values()) == {"NOT_RUN"} and "unfinished" in err, (err, statuses(out)))
k = mock_lines.index("# Query: gene_c7 mock candidate gene") - 1     # gene_c7's header line: 5 queries before it, 8 after
first, second = mock_lines[:k], mock_lines[k:-1]
n1 = sum(1 for ln in first if ln.startswith("# Query:"))
n2 = sum(1 for ln in second if ln.startswith("# Query:"))
rc, err, out = run([], mock=joined("two_ok", first + ["# BLAST processed %d queries" % n1] + second + ["# BLAST processed %d queries" % n2]))
check("two finished runs joined in one file are both trusted", statuses(out) == BASE, statuses(out))
rc, err, out = run([], mock=joined("second_cut", first + ["# BLAST processed %d queries" % n1] + second))
s = statuses(out)
check("the second run was cut short: the first run's answers stand, the second run's genes are NOT_RUN",
      s["gene_c1"] == "HIT" and s["gene_c2"] == "IN_PHYLOGENY_ONLY" and s["gene_c5"] == "NO_HIT"
      and all(s[g] == "NOT_RUN" for g in ("gene_c7", "gene_c8", "gene_c9", "gene_c14")), s)
rc, err, out = run([], mock=joined("first_cut", first + second + ["# BLAST processed %d queries" % n2]))
check("the first run was cut short and a finished run follows: no gene is trusted", set(statuses(out).values()) == {"NOT_RUN"}, statuses(out))
rc, err, out = run([], mock=joined("wrong_count", mock_lines[:-1] + ["# BLAST processed 99 queries"]))
check("a footer whose query count does not match: no gene is trusted", set(statuses(out).values()) == {"NOT_RUN"}, statuses(out))
bad = [ln if ln.startswith("#") else "\t".join(ln.split("\t")[:9]) for ln in mock_lines]
rc, err, out = run([], mock=joined("badcols", bad))
check("wrong column count stops the run and says which format is expected", rc != 0 and "expected 14" in err and "outfmt" in err, err)
empty = os.path.join(tmp, "nothing")
os.makedirs(empty)
rc, err, out = run([], mock=empty)
check("no result file for a species: all its genes are NOT_RUN, and the exit status is 2", rc == 2 and set(statuses(out).values()) == {"NOT_RUN"}, statuses(out))
p = subprocess.run([sys.executable, SCRIPT, "-i", INPUT, "-o", os.path.join(tmp, "nodb")], stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, universal_newlines=True)
check("--db is required unless --parse-only", p.returncode != 0 and "--db is required" in p.stderr, p.stderr)

# --- output-folder protection reaches this module
again = [sys.executable, SCRIPT, "-i", INPUT, "-o", run([])[2], "--parse-only", MOCK, "--quiet"]
rc, err, out = run([])
p = subprocess.run([sys.executable, SCRIPT, "-i", INPUT, "-o", out, "--parse-only", MOCK, "--quiet"],
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
check("running again into the same output folder is refused", p.returncode != 0 and "already holds results" in p.stderr, p.stderr)

# --- pieces
check("taxon map: comments and both separators", (lambda p: (open(p, "w").write("# c\nA\tSp one\nB,Sp two\n"), b.read_taxon_map(p))[1])
      (os.path.join(tmp, "m.tsv")) == {"A": "Sp one", "B": "Sp two"})
check("db_label strips the FASTA extension and odd characters", b.db_label("/x/my db.fasta") == "my_db" and b.db_label("/x/prot") == "prot")
cmd = b.blast_command("blastp", "q.faa", "/db/p", "o.tsv", {"evalue": 0.001, "max_target_seqs": 5000, "seg": True,
                                                            "max_hsps": None, "ungapped": False, "comp_based_stats": None,
                                                            "dbsize": None}, 4)
check("blast command line: program, database, E-value, format, no truncation, low-complexity masking on",
      cmd[0] == "blastp" and cmd[cmd.index("-db") + 1] == "/db/p" and cmd[cmd.index("-evalue") + 1] == "0.001"
      and cmd[cmd.index("-outfmt") + 1].startswith("7 qseqid") and cmd[cmd.index("-max_target_seqs") + 1] == "5000"
      and cmd[cmd.index("-seg") + 1] == "yes" and cmd[cmd.index("-num_threads") + 1] == "4"
      and "-ungapped" not in cmd and "-max_hsps" not in cmd, cmd)
check("the network is never involved: no -remote option can be produced",
      "-remote" not in cmd and not any("remote" in k for k in h.load_defaults("homology-blast")))

# --- the run path, with stand-in blastp / tblastn / blastdbcmd (Linux and WSL only)
if os.name == "nt":
    print("SKIP run-path checks with stand-in BLAST+ programs: need Linux or WSL")
else:
    bindir = os.path.join(tmp, "fakebin")
    os.makedirs(bindir)
    blast_stub = r'''#!/usr/bin/env bash
q=""; out=""; db=""; args="$*"
while [[ $# -gt 0 ]]; do case "$1" in -query) q="$2"; shift 2;; -out) out="$2"; shift 2;; -db) db="$2"; shift 2;; -version) echo "%NAME%: 2.99.0+"; exit 0;; -evalue|-outfmt|-max_target_seqs|-num_threads|-seg|-max_hsps|-comp_based_stats|-dbsize) shift 2;; *) shift;; esac; done
echo "%NAME% $args" >> "$FAKE_CALLS"
cp "$q" "$FAKE_SEEN/$(basename "$q")"
if grep -q '^>bad' "$q"; then echo "boom" >&2; exit 1; fi
# which genes hit this database: by default every gene except those named nohit*; HITS_<dbname>="g1 g2" overrides that
dbbase="$(basename "$db")"; var="HITS_${dbbase}"; mode=all; list=""
if [[ -n "${!var+x}" ]]; then mode=list; list="${!var}"; fi
awk -v P="%UPPER%" -v MODE="$mode" -v L=" $list " '
function flush() { if (id == "") return; nq++;
  print "# " P " 2.99.0+"; print "# Query: " id; print "# Database: fake";
  if ((MODE == "all" && id ~ /^nohit/) || (MODE == "list" && index(L, " " id " ") == 0)) { print "# 0 hits found" } else {
    print "# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, query length, subject length";
    print "# 1 hits found"; printf "%s\t%s\t60.0\t%d\t0\t0\t1\t%d\t1\t%d\t1e-30\t100\t%d\t%d\n", id, "T__Ceratodon_purpureus", len, len, len, len, len } }
/^>/ { flush(); id = substr($1, 2); len = 0; next } { len += length($0) }
END { flush(); if (ENVIRON["FAKE_NOFOOTER"] == "") print "# BLAST processed " nq " queries" }' "$q" > "$out"
exit 0
'''
    for name in ("blastp", "tblastn"):
        p_ = os.path.join(bindir, name)
        open(p_, "w").write(blast_stub.replace("%NAME%", name).replace("%UPPER%", name.upper()))
        os.chmod(p_, 0o755)
    dbcmd = os.path.join(bindir, "blastdbcmd")
    open(dbcmd, "w").write(r'''#!/usr/bin/env bash
db=""; dbtype=""; mode=""
while [[ $# -gt 0 ]]; do case "$1" in -db) db="$2"; shift 2;; -dbtype) dbtype="$2"; shift 2;; -info) mode=info; shift;; -entry) mode=entry; shift 2;; *) shift;; esac; done
[[ "$db" == *broken* ]] && { echo "BLAST Database error: No alias or index file found"; exit 1; }
want="${FAKE_DBTYPE:-prot}"
if [[ "$mode" == info ]]; then
  [[ "$dbtype" == "$want" ]] || { echo "BLAST Database error: no $dbtype database"; exit 1; }
  printf 'Database: fake db\n\t1,234 sequences; 99,999 total residues\n'; exit 0
fi
case "${FAKE_IDS:-tagged}" in
  tagged) for i in $(seq 1 30); do echo "s${i}__Ceratodon_purpureus"; done;;
  untagged) for i in $(seq 1 30); do echo "s${i}"; done;;
  ord) for i in $(seq 1 30); do echo "gnl|BL_ORD_ID|${i}"; done;;
  self) for i in $(seq 1 10); do echo "s${i}__Sp"; done; for i in $(seq 1 5); do echo "c${i}__Ceratodon_purpureus"; done;;
  mixed) for i in $(seq 1 6); do echo "s${i}__Ceratodon_purpureus"; done; for i in $(seq 1 4); do echo "u${i}"; done;;
esac
''')
    os.chmod(dbcmd, 0o755)
    mkdb = os.path.join(bindir, "makeblastdb")
    open(mkdb, "w").write(r'''#!/usr/bin/env bash
in=""; out=""; dbtype=""; allargs="$*"
while [[ $# -gt 0 ]]; do case "$1" in -in) in="$2"; shift 2;; -out) out="$2"; shift 2;; -dbtype) dbtype="$2"; shift 2;; -title) shift 2;; *) shift;; esac; done
echo "makeblastdb $allargs" >> "$FAKE_CALLS"
if [[ -n "${FAKE_MAKEDB_FAIL:-}" ]]; then echo "BLAST Database creation error: bad sequence" >&2; exit 1; fi
if [[ "$dbtype" == prot ]]; then : > "${out}.pin"; : > "${out}.psq"; : > "${out}.phr"; else : > "${out}.nin"; : > "${out}.nsq"; : > "${out}.nhr"; fi
exit 0
''')
    os.chmod(mkdb, 0o755)
    # a second folder of stand-ins WITHOUT makeblastdb, to test its absence
    bindir2 = os.path.join(tmp, "fakebin2")
    os.makedirs(bindir2)
    for n_ in ("blastp", "tblastn", "blastdbcmd"):
        os.symlink(os.path.join(bindir, n_), os.path.join(bindir2, n_))

    def make_input(name, seqs):
        p_ = os.path.join(tmp, name)
        os.makedirs(p_)
        with open(os.path.join(p_, "Sp_final.faa"), "w") as fh:
            for gid, seq in seqs:
                fh.write(">%s\n%s\n" % (gid, seq))
        return p_

    def run_real(name, seqs, dbs=("fake_db",), extra=(), env_extra=None):
        inp = make_input(name + "_in", seqs)
        out = os.path.join(tmp, name + "_out")
        calls, seen = os.path.join(tmp, name + "_calls.txt"), os.path.join(tmp, name + "_seen")
        os.makedirs(seen)
        env = dict(os.environ, PATH=bindir + os.pathsep + os.environ["PATH"], FAKE_CALLS=calls, FAKE_SEEN=seen)
        env.update(env_extra or {})
        cmd = [sys.executable, SCRIPT, "-i", inp, "-o", out, "--quiet", "--set", "chunk_size=2", "--threads", "3"]
        for d in dbs:
            cmd += ["--db", d]
        p = subprocess.run(cmd + list(extra), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
        return p, out, calls, seen

    seqs = [("g1", "MKVLAA"), ("nohit2", "MKVL"), ("g3", "MKVLAAGGG"), ("g4", "MK"), ("g5", "MKVLA")]
    p, out, calls, seen = run_real("ok", seqs)
    check("run path succeeds", p.returncode == 0, p.stderr)
    check("run path: genes with hits HIT (the hit is to another species), the gene with none NO_HIT",
          statuses(out) == {"g1": "HIT", "nohit2": "NO_HIT", "g3": "HIT", "g4": "HIT", "g5": "HIT"}, statuses(out))
    call_lines = open(calls).read().splitlines()
    check("run path: 5 genes in chunks of 2 make 3 BLAST calls", len(call_lines) == 3, call_lines)
    check("run path: every call has the E-value, no truncation, masking on, the fixed output format and the thread count",
          all("-evalue 0.001" in c and "-max_target_seqs 5000" in c and "-seg yes" in c and "-outfmt 7 qseqid" in c
              and "-num_threads 3" in c for c in call_lines), call_lines)
    rj = json.load(open(os.path.join(out, "run.json")))
    check("run path: run.json records BLAST+'s version, the database, and its species inventory (every ID counted)",
          rj["tools"][0] == {"name": "BLASTP (BLAST+)", "version": "2.99.0+"} and rj["database_info"]["type"] == "prot"
          and rj["database_info"]["sequences"] == 1234 and rj["database_inventory"]["ids"] == 30
          and rj["database_inventory"]["tagged"] == 30 and rj["database_inventory"]["species"] == {"Ceratodon_purpureus": 30}, rj)
    check("run path: the raw results are kept so the run can be re-read", os.path.isfile(os.path.join(out, "raw", "Sp.blast.tsv")))
    rc, err, out2 = run([], mock=os.path.join(out, "raw"), inp=os.path.join(tmp, "ok_in"))
    check("run path: the saved raw results can be re-read with --parse-only", rc == 0 and statuses(out2) == statuses(out), err)

    # chunks are longest first: [g1 (6), g5 (5)], [g3 (4), bad2 (3)], [g4 (2)]
    p, out, calls, seen = run_real("fail", [("g1", "MKVLAA"), ("bad2", "MKV"), ("g3", "MKVL"), ("g4", "MK"), ("g5", "MKVLA")])
    s = statuses(out)
    check("a chunk that fails: its genes are NOT_RUN; the other chunks keep their answers",
          p.returncode == 2 and s["bad2"] == "NOT_RUN" and s["g3"] == "NOT_RUN" and s["g1"] == "HIT" and s["g5"] == "HIT"
          and s["g4"] == "HIT", (p.stderr, s))
    p, out, calls, seen = run_real("nofoot", seqs, env_extra={"FAKE_NOFOOTER": "1"})
    check("output without BLAST's completion footer is NOT trusted, even with exit code 0: all genes NOT_RUN",
          p.returncode == 2 and set(statuses(out).values()) == {"NOT_RUN"}, (p.stderr, statuses(out)))

    p, out, calls, seen = run_real("nucl", seqs, env_extra={"FAKE_DBTYPE": "nucl"})
    check("a nucleotide database with blastp is refused, and tblastn is suggested",
          p.returncode != 0 and "needs a protein database" in p.stderr and "tblastn" in p.stderr, p.stderr)
    p, out, calls, seen = run_real("tbn", seqs, extra=["--program", "tblastn"], env_extra={"FAKE_DBTYPE": "nucl"})
    check("tblastn against a nucleotide database runs", p.returncode == 0 and "tblastn" in open(calls).read()
          and json.load(open(os.path.join(out, "run.json")))["program"] == "tblastn", p.stderr)
    p, out, calls, seen = run_real("untag", seqs, env_extra={"FAKE_IDS": "untagged"})
    check("a database whose IDs carry no species tag: a clear warning that such hits will all be counted",
          p.returncode == 0 and "none of the 30 sequence IDs in this database carries a '__Species_Name' tag" in p.stderr, p.stderr)
    p, out, calls, seen = run_real("mixed", seqs, env_extra={"FAKE_IDS": "mixed"})
    check("a database where only some IDs carry a species tag: the untagged share is reported",
          p.returncode == 0 and "4 of 10 sequence IDs carry no species tag" in p.stderr, p.stderr)
    # the analysis species (input file Sp_final.faa is species 'Sp') sitting inside the database: named in a warning
    p, out, calls, seen = run_real("inside", seqs, env_extra={"FAKE_IDS": "self"})
    check("the analysis species inside the database is named, with how many sequences, and the consequence is stated",
          p.returncode == 0 and "INSIDE this database: Sp (10 sequences)" in p.stderr
          and "reported but not counted (IN_PHYLOGENY_ONLY)" in p.stderr, p.stderr)
    p, out, calls, seen = run_real("inside2", seqs, extra=["--ignore-in-phylogeny", "no"], env_extra={"FAKE_IDS": "self"})
    check("...and when the exemption is off, the warning says those hits WILL be counted", "WILL be counted" in p.stderr, p.stderr)
    p, out, calls, seen = run_real("skipinv", seqs, extra=["--skip-db-inventory"], env_extra={"FAKE_IDS": "self"})
    check("--skip-db-inventory reads no IDs: no warning, and run.json says the inventory was skipped",
          p.returncode == 0 and "INSIDE" not in p.stderr and json.load(open(os.path.join(out, "run.json")))["database_inventory"] is None, p.stderr)
    p, out, calls, seen = run_real("badfolder", seqs, env_extra={"FAKE_DBTYPE": "nucl"})
    check("a database that cannot be used leaves no output folder behind (no header-only tables)",
          p.returncode != 0 and not os.path.exists(out), (p.stderr, os.path.exists(out)))
    p, out, calls, seen = run_real("ord", seqs, env_extra={"FAKE_IDS": "ord"})
    check("a database built without -parse_seqids is refused, with the fix named",
          p.returncode != 0 and "-parse_seqids" in p.stderr, p.stderr)
    p, out, calls, seen = run_real("flags", seqs, extra=["--set", "ungapped=yes", "--set", "comp_based_stats=0",
                                                         "--set", "max_hsps=5", "--set", "dbsize=1000000", "--seg", "no"])
    c = open(calls).read()
    check("optional BLAST+ settings are passed through", p.returncode == 0 and "-ungapped" in c and "-comp_based_stats 0" in c
          and "-max_hsps 5" in c and "-dbsize 1000000" in c and "-seg no" in c, (p.stderr, c))

    # several databases in one command: each its own labelled run; one broken database does not stop the other
    p, out, calls, seen = run_real("multi", seqs, dbs=("/x/alpha", "/x/broken_db", "/x/beta"))
    check("several databases: a broken one is reported, the others still complete, and the exit code says something failed",
          p.returncode != 0 and "1 of 3 database run(s) did not complete: broken_db" in p.stderr
          and os.path.isfile(os.path.join(out, "alpha", "gene_status.tsv")) and os.path.isfile(os.path.join(out, "beta", "gene_status.tsv"))
          and not os.path.exists(os.path.join(out, "broken_db", "gene_status.tsv")), p.stderr)
    check("several databases: each run is labelled by its database and its hits say which database they came from",
          json.load(open(os.path.join(out, "alpha", "run.json")))["label"] == "alpha"
          and {r["source"] for r in table(os.path.join(out, "beta", "hits.tsv"))} == {"beta"})
    p, out, calls, seen = run_real("multilabel", seqs, dbs=("/x/alpha", "/x/beta"), extra=["--label", "prot"])
    check("several databases with --label: labels are prefixed", os.path.isdir(os.path.join(out, "prot-alpha")) and os.path.isdir(os.path.join(out, "prot-beta")), p.stderr)
    p, out, calls, seen = run_real("samelabel", seqs, dbs=("/x/a/prot", "/y/b/prot"))
    check("two databases that would get the same label are refused", p.returncode != 0 and "same label" in p.stderr, p.stderr)
    # --- iterative filtering across databases in ONE command (--sequential)
    seq5 = [("g1", "MKVLAA"), ("g2", "MKVLAAG"), ("g3", "MKVLAAGG"), ("g4", "MKVLAAGGG"), ("g5", "MKVL")]
    hits_env = {"HITS_alpha": "g1 g2", "HITS_beta": "g2 g3", "HITS_gamma": "g4"}
    dbs3 = ("/x/alpha", "/x/beta", "/x/gamma")
    p, out, calls, seen = run_real("par", seq5, dbs=dbs3, env_extra=hits_env)
    par = {lab: statuses(os.path.join(out, lab)) for lab in ("alpha", "beta", "gamma")}
    check("WITHOUT --sequential every database searches every gene (independent runs)",
          p.returncode == 0 and all(set(v) == {"g1", "g2", "g3", "g4", "g5"} for v in par.values())
          and par["beta"]["g2"] == "HIT" and len(open(calls).read().splitlines()) == 9, (p.stderr, par))
    p, out, calls, seen = run_real("seq", seq5, dbs=dbs3, extra=["--sequential"], env_extra=hits_env)
    seq = {lab: statuses(os.path.join(out, lab)) for lab in ("alpha", "beta", "gamma")}
    check("--sequential succeeds", p.returncode == 0, p.stderr)
    check("--sequential: the first database searches all 5 genes (g1 and g2 hit it)",
          seq["alpha"] == {"g1": "HIT", "g2": "HIT", "g3": "NO_HIT", "g4": "NO_HIT", "g5": "NO_HIT"}, seq["alpha"])
    check("--sequential: the second database searches only the 3 survivors; g2, which already hit, is never searched again",
          seq["beta"] == {"g3": "HIT", "g4": "NO_HIT", "g5": "NO_HIT"}, seq["beta"])
    check("--sequential: the third database searches only the 2 survivors of the second",
          seq["gamma"] == {"g4": "HIT", "g5": "NO_HIT"}, seq["gamma"])
    check("--sequential: the genes that passed every database are the last run's survivors (g5 only)",
          open(os.path.join(out, "gamma", "passed_ids.tsv")).read().split() == ["species", "gene_id", "Sp", "g5"])
    check("--sequential: fewer BLAST calls than searching everything everywhere (6 instead of 9)",
          len(open(calls).read().splitlines()) == 6, open(calls).read())
    rj = json.load(open(os.path.join(out, "beta", "run.json")))
    check("--sequential: run.json records the position in the chain and where the input came from",
          rj["chain"] == {"sequential": True, "position": 2, "of": 3, "input_from": "alpha"}
          and json.load(open(os.path.join(out, "alpha", "run.json")))["chain"]["input_from"] is None, rj["chain"])
    inputs_n = {lab: [ln.split("\t") for ln in open(os.path.join(out, lab, "summary.tsv")).read().splitlines()][1][1] for lab in ("alpha", "beta", "gamma")}
    check("--sequential: summary.tsv shows how many genes each database received (5, 3, 2)", inputs_n == {"alpha": "5", "beta": "3", "gamma": "2"}, inputs_n)
    check("--sequential: the end says where the final survivors are", "genes that passed every database are in" in p.stderr, p.stderr)

    p, out, calls, seen = run_real("seqfail", seq5, dbs=("/x/alpha", "/x/broken_db", "/x/gamma"), extra=["--sequential"], env_extra=hits_env)
    check("--sequential with a broken database: the chain STOPS, and says which databases were not searched",
          p.returncode != 0 and "NOT searched: gamma" in p.stderr and "cannot pass a database that was never searched" in p.stderr, p.stderr)
    check("--sequential with a broken database: later databases produce no results, so nothing passes unchecked",
          os.path.isfile(os.path.join(out, "alpha", "gene_status.tsv")) and not os.path.exists(os.path.join(out, "gamma")), os.listdir(out))
    p, out, calls, seen = run_real("seqempty", seq5, dbs=dbs3, extra=["--sequential"],
                                   env_extra={"HITS_alpha": "g1 g2 g3 g4 g5", "HITS_beta": "g1", "HITS_gamma": "g1"})
    empty_in = {lab: [ln.split("\t") for ln in open(os.path.join(out, lab, "summary.tsv")).read().splitlines()][-1][1] for lab in ("beta", "gamma")}
    check("--sequential when the first database removes every gene: the others still produce (empty) results and run no BLAST",
          p.returncode == 0 and empty_in == {"beta": "0", "gamma": "0"} and len(open(calls).read().splitlines()) == 3, (p.stderr, empty_in))
    p, out, calls, seen = run_real("seqone", seq5, dbs=("/x/alpha",), extra=["--sequential"], env_extra=hits_env)
    check("--sequential with one database has no effect and says so", p.returncode == 0 and "needs at least two --db" in p.stderr, p.stderr)
    p = subprocess.run([sys.executable, SCRIPT, "-i", INPUT, "-o", os.path.join(tmp, "seqparse"), "--parse-only", MOCK, "--sequential", "--quiet"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    check("--sequential with --parse-only is refused", p.returncode != 0 and "cannot be used with --parse-only" in p.stderr, p.stderr)

    # --- a FASTA file as the database: the BLAST index is built from it (makeblastdb, local), cached, and reused
    fasta_db = os.path.join(tmp, "proteins.faa")
    open(fasta_db, "w").write(">a1__Ceratodon_purpureus\nMKVLAA\n>a2__Blasia_pusilla\nMKV\n")
    p, out, calls, seen = run_real("idx1", seqs, dbs=(fasta_db,))
    mk = [ln for ln in open(calls).read().splitlines() if ln.startswith("makeblastdb")]
    idxdir = os.path.join(out, "db_index")
    check("a FASTA file as --db: the index is built once, with -parse_seqids and the right type, inside the run's output folder",
          p.returncode == 0 and len(mk) == 1 and "-dbtype prot" in mk[0] and "-parse_seqids" in mk[0]
          and ("-in " + fasta_db) in mk[0] and idxdir in mk[0] and os.path.isdir(idxdir), (p.stderr, mk))
    check("a FASTA file as --db: run.json records the FASTA and where the index is",
          json.load(open(os.path.join(out, "run.json")))["database"] == os.path.abspath(fasta_db)
          and "db_index" in json.load(open(os.path.join(out, "run.json")))["database_index"])
    shared = os.path.join(tmp, "shared_index")
    p, out, calls1, seen = run_real("idx2", seqs, dbs=(fasta_db,), extra=["--index-dir", shared])
    p2, out2, calls2, seen2 = run_real("idx3", seqs, dbs=(fasta_db,), extra=["--index-dir", shared])
    n_built = sum(1 for c in (calls1, calls2) for ln in open(c).read().splitlines() if ln.startswith("makeblastdb"))
    check("a shared --index-dir: the second run reuses the index instead of building it again",
          p.returncode == 0 and p2.returncode == 0 and n_built == 1 and "reusing the BLAST index" in p2.stderr, (p2.stderr, n_built))
    open(fasta_db, "a").write(">a3__Blasia_pusilla\nMKVLAAGG\n")
    p3, out3, calls3, seen3 = run_real("idx4", seqs, dbs=(fasta_db,), extra=["--index-dir", shared])
    check("...and rebuilds it when the FASTA has changed",
          p3.returncode == 0 and any(ln.startswith("makeblastdb") for ln in open(calls3).read().splitlines()), p3.stderr)
    p, out, calls, seen = run_real("idxn", seqs, dbs=(fasta_db,), extra=["--program", "tblastn"], env_extra={"FAKE_DBTYPE": "nucl"})
    check("a FASTA as the database for tblastn is indexed as nucleotide",
          p.returncode == 0 and any("-dbtype nucl" in ln for ln in open(calls).read().splitlines()), p.stderr)
    p, out, calls, seen = run_real("idxfail", seqs, dbs=(fasta_db,), env_extra={"FAKE_MAKEDB_FAIL": "1"})
    check("makeblastdb failing is reported with its own message, and nothing is left behind",
          p.returncode != 0 and "makeblastdb failed" in p.stderr and "BLAST Database creation error" in p.stderr
          and not os.path.exists(out), p.stderr)
    p, out, calls, seen = run_real("idxnobin", seqs, dbs=(fasta_db,), env_extra={"PATH": bindir2 + os.pathsep + "/usr/bin:/bin"})
    check("indexing a FASTA without makeblastdb on PATH is reported", p.returncode != 0 and "makeblastdb not found" in p.stderr, p.stderr)

    p, out, calls, seen = run_real("nobin", seqs, env_extra={"PATH": "/usr/bin:/bin"})
    check("BLAST+ missing from PATH is reported", p.returncode != 0 and "not found in PATH" in p.stderr, p.stderr)

print()
print("FAILED: %d" % len(fails) if fails else "ALL PASSED")
sys.exit(1 if fails else 0)
