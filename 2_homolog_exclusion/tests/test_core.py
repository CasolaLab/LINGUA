"""
test_core.py - checks for the shared code in bin/hep_common.py, using small made-up FASTA and config files.

Run from anywhere:   python3 tests/test_core.py
It needs no external tools and no data, writes only to a temporary folder, and exits 0 if every check passes.
Each line printed is one check (PASS or FAIL). Run it first on a new machine to confirm the environment works.
"""
import json, os, sys, tempfile, io, contextlib
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "bin"))
import hep_common as h

fails = []
def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + (("  " + str(extra)) if (extra and not cond) else ""))
    if not cond:
        fails.append(name)

def raises(exc, fn, *a, **k):
    try:
        fn(*a, **k)
    except exc:
        return True
    except Exception as e:
        return "wrong exception: %r" % e
    return False

tmp = tempfile.mkdtemp()
def w(name, text):
    p = os.path.join(tmp, name)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w").write(text)
    return p

# --- coerce / as_list
check("coerce float", h.coerce("1e-3") == 0.001)
check("coerce int", h.coerce("50") == 50 and isinstance(h.coerce("50"), int))
check("coerce blank", h.coerce("  ") is None)
check("coerce bool", h.coerce("True") is True)
check("as_list", h.as_list("a, b,,c") == ["a", "b", "c"] and h.as_list(None) == [])

# --- defaults and config resolution
d = h.load_defaults("homology-jackhmmer")
check("defaults jackhmmer", d["iterations"] == 3 and d["incE"] == 0.001 and d["dbsize"] is None, d)
check("common merged", d["species_strip_suffixes"] == "_final")
preset = w("preset.conf", "[domain-cdd]\nevalue=1e-5\n[homology-jackhmmer]\nevalue=1e-4\nincE=1e-5\niterations=1\n")
cfg, nd = h.resolve_config("homology-jackhmmer", config_path=preset)
check("config jackhmmer", cfg["iterations"] == 1 and cfg["evalue"] == 1e-4 and cfg["incE"] == 1e-5, cfg)
check("non_default recorded", set(nd) == {"iterations", "evalue", "incE"}, nd)
cfg, nd = h.resolve_config("domain-cdd", config_path=preset)
check("config cdd", cfg["evalue"] == 1e-5 and cfg["min_domain_cov"] == 50)
cfg, nd = h.resolve_config("domain-cdd", config_path=preset, overrides={"evalue": "1e-8"})
check("override beats config", cfg["evalue"] == 1e-8)
cfg, nd = h.resolve_config("homology-blast", config_path=preset)
check("section for another module ignored", nd == {}, nd)
check("unknown override rejected", raises(h.HepError, h.resolve_config, "domain-cdd", None, {"evalu": "1"}) is True)
bad = w("bad.conf", "[domain-cdd]\nevalu=1\n")
check("unknown section key rejected", raises(h.HepError, h.resolve_config, "domain-cdd", bad) is True)
bad2 = w("bad2.conf", "not a pair\n")
check("bad line rejected", raises(h.HepError, h.resolve_config, "domain-cdd", bad2) is True)
glob_conf = w("glob.conf", "evalue=0.5\nfoo=1\n")
cfg, nd = h.resolve_config("domain-cdd", glob_conf)
check("global key applies, unknown global ignored", cfg["evalue"] == 0.5 and "foo" not in cfg)


# --- new defaults and the shipped organelle example
b = h.load_defaults("homology-blast")
check("blast defaults", b["program"] == "blastp" and b["evalue"] == 0.001 and b["min_qcov"] == 50 and b["min_tcov"] == 0
      and b["seg"] is True and b["ignore_in_phylogeny"] is True and b["ungapped"] is False and b["max_hsps"] is None, b)
org = os.path.join(h.HERE, os.pardir, "presets", "organelle_tblastn.conf")
cfg, nd = h.resolve_config("homology-blast", config_path=org)
check("organelle example", cfg["program"] == "tblastn" and cfg["ignore_in_phylogeny"] is False and cfg["ungapped"] is True
      and cfg["comp_based_stats"] == 0 and cfg["max_hsps"] == 5 and cfg["min_qcov"] == 0, cfg)
check("organelle example does not touch other modules", h.resolve_config("domain-cdd", config_path=org)[1] == {})
i = h.load_defaults("domain-interpro")
check("interpro applications: the lab's 12 databases plus MobiDBLite, Coils, ProSitePatterns and AntiFam",
      len(h.as_list(i["applications"])) == 16 and {"Pfam", "CDD", "MobiDBLite", "Coils", "ProSitePatterns", "AntiFam"} <= set(h.as_list(i["applications"])),
      i["applications"])
check("interpro ignore_domains empty", i["ignore_domains"] is None and h.as_list(i["ignore_domains"]) == [])
check("old module name 'homology-blastp' no longer exists", raises(h.HepError, h.load_defaults, "homology-blastp") is True)
j = h.load_defaults("homology-jackhmmer")
check("jackhmmer ignore_in_phylogeny", j["ignore_in_phylogeny"] is True)
check("match_level defaults to species", b["match_level"] == "species" and j["match_level"] == "species")
check("cdd ignore_sources empty by default", h.load_defaults("domain-cdd")["ignore_sources"] is None)

# --- ID handling
check("split_id suffix", h.split_id("Aa_G1.t1__Arabis_alpina") == ("Aa_G1.t1", "Arabis_alpina"))
check("split_id none", h.split_id("gene_a1") == ("gene_a1", None))
check("split_id pipes", h.split_id("sp|P12345|PROT_HUMAN__Homo_sapiens") == ("sp|P12345|PROT_HUMAN", "Homo_sapiens"))
check("split_id edge", h.split_id("__x") == ("__x", None) and h.split_id("x__") == ("x__", None))

# --- FASTA
fa = w("in/SpeciesA_final.faa",
       ">gene_a1 desc one\nMKV\nLLA\n>gene_a2__SpeciesA\nMMMM\n\n>gene_a3\nAAAAAAAAAA\n")
recs = list(h.read_fasta(fa))
check("read_fasta", [r[0] for r in recs] == ["gene_a1", "gene_a2__SpeciesA", "gene_a3"] and recs[0][2] == "MKVLLA")
check("load ids canonical", h.load_fasta_ids(fa) == ["gene_a1", "gene_a2", "gene_a3"])
dup = w("dup.faa", ">g1\nAA\n>g1__X\nAA\n")
check("duplicate ids rejected", raises(h.HepError, h.load_fasta_ids, dup) is True)
nohdr = w("nohdr.faa", "AAAA\n>g1\nAA\n")
check("sequence before header rejected", raises(h.HepError, list, h.read_fasta(nohdr)) is True)
chunks = h.chunk_fasta(fa, os.path.join(tmp, "chunks"), 2)
check("chunk_fasta", len(chunks) == 2 and len(list(h.read_fasta(chunks[0]))) == 2 and len(list(h.read_fasta(chunks[1]))) == 1)

lens = w("lens.faa", ">s1\nAA\n>s2\nAAAAAAAA\n>s3\nAAAA\n>s4\nAAAAAA\n>s5\nAAAAAAAAAA\n")
ch = h.chunk_fasta(lens, os.path.join(tmp, "chunks_sorted"), 2)
first = [r[0] for r in h.read_fasta(ch[0])]
check("chunks are longest first: chunk 0 holds the two longest sequences", first == ["s5", "s2"] and len(ch) == 3, first)
check("last chunk holds the shortest", [r[0] for r in h.read_fasta(ch[-1])] == ["s1"])
ch = h.chunk_fasta(lens, os.path.join(tmp, "chunks_plain"), 2, sort_by_length=False)
check("sorting can be switched off (input order kept)", [r[0] for r in h.read_fasta(ch[0])] == ["s1", "s2"])
tie = w("tie.faa", ">t1\nAAAA\n>t2\nAAAA\n>t3\nAAAA\n")
ch = h.chunk_fasta(tie, os.path.join(tmp, "chunks_tie"), 3)
check("equal lengths keep the input order", [r[0] for r in h.read_fasta(ch[0])] == ["t1", "t2", "t3"])

# --- species levels and the taxon matcher (subspecies, strains, genus)
check("taxon_key species level: subspecies and strain suffixes are dropped",
      h.taxon_key("Marchantia polymorpha subsp. ruderalis") == h.taxon_key("Marchantia_polymorpha") == "marchantia_polymorpha"
      and h.taxon_key("Ceratodon_purpureus_R40") == "ceratodon_purpureus")
check("taxon_key exact level keeps the whole name",
      h.taxon_key("Marchantia polymorpha subsp. ruderalis", "exact") == "marchantia_polymorpha_subsp._ruderalis")
check("taxon_key genus level", h.taxon_key("Marchantia polymorpha", "genus") == "marchantia" and h.taxon_key("Sphagnum", "genus") == "sphagnum")
check("an unknown match level is refused", raises(h.HepError, h.taxon_key, "A b", "family") is True)
m = h.TaxonMatcher(["Ceratodon purpureus R40", "Marchantia polymorpha subsp. ruderalis"], ["Sphagnum", "Physcomitrium patens"])
check("matcher (species level): a database tag without the strain matches the analysis species with it",
      m.match("Ceratodon_purpureus") == "Ceratodon purpureus R40" and m.match("Marchantia_polymorpha") is not None)
check("matcher: a different species of the same genus does not match at species level",
      m.match("Ceratodon_other") is None and m.match("Marchantia_paleacea") is None)
check("matcher: a single word in ignore_taxa is a genus and matches every species of it",
      m.match("Sphagnum_fallax") == "Sphagnum" and m.match("Sphagnum_palustre") == "Sphagnum")
check("matcher: a two-word ignore_taxa entry is a species", m.match("Physcomitrium_patens") == "Physcomitrium patens")
check("matcher: an unknown species never matches, so the hit counts", m.match(None) is None and m.match("") is None)
mx = h.TaxonMatcher(["Ceratodon purpureus R40"], [], "exact")
check("matcher (exact level): the strain must match too", mx.match("Ceratodon_purpureus") is None and mx.match("Ceratodon_purpureus_R40") is not None)
mg = h.TaxonMatcher(["Ceratodon purpureus"], [], "genus")
check("matcher (genus level): any species of the genus matches", mg.match("Ceratodon_other") is not None and mg.match("Marchantia_x") is None)

# --- species
check("species strip _final", h.species_from_path("x/Arabis_alpina_final.faa") == "Arabis_alpina")
check("species plain", h.species_from_path("x/SpeciesB.fasta") == "SpeciesB")
check("species custom suffix none", h.species_from_path("x/A_final.faa", strip_suffixes=()) == "A_final")
smap = w("map.tsv", "Species\tBasename\nArabis alpina\tAa_proteins\nBrassica rapa\tBr.faa\n")
m = h.load_species_map(smap)
check("species map", m == {"Aa_proteins": "Arabis alpina", "Br": "Brassica rapa"}, m)
check("species from map", h.species_from_path("z/Aa_proteins.faa", species_map=m) == "Arabis alpina")
check("species_key", h.species_key("Arabis  alpina") == "arabis_alpina" == h.species_key("Arabis_Alpina"))
w("in/SpeciesB.faa", ">gene_b1\nMK\n")
ins = h.discover_inputs(os.path.join(tmp, "in"))
check("discover_inputs", [s for s, _ in ins] == ["SpeciesA", "SpeciesB"], ins)
w("clash/A_final.faa", ">g\nA\n"); w("clash/A.faa", ">g\nA\n")
check("species clash rejected", raises(h.HepError, h.discover_inputs, os.path.join(tmp, "clash")) is True)
check("empty dir rejected", raises(h.HepError, h.discover_inputs, os.path.join(tmp, "nothing")) is True)

# --- classify / passes
st = h.classify_genes(["a", "b", "c", "d", "e", "f"], counted=["a", "f"], excluded=["b"], spurious=["c", "a"],
                      in_phylogeny=["d"], not_run=["f"])
check("classify", st == {"a": "HIT", "b": "EXCLUDED_ONLY", "c": "SPURIOUS", "d": "IN_PHYLOGENY_ONLY",
                         "e": "NO_HIT", "f": "NOT_RUN"}, st)
check("passes", h.passes("NO_HIT") and h.passes("EXCLUDED_ONLY") and h.passes("IN_PHYLOGENY_ONLY")
      and not h.passes("HIT") and not h.passes("NOT_RUN")
      and h.passes("SPURIOUS", "flag") and not h.passes("SPURIOUS", "remove"))

# --- ModuleOutput
out = os.path.join(tmp, "out", "homology-blast")
with h.ModuleOutput(out, "homology-blast", label="organelle") as mo:
    tally = mo.add_species("SpeciesA", fa,
                           [{"gene_id": "gene_a1", "target": "T1", "evalue": 1e-9, "counts_as_hit": True},
                            {"gene_id": "gene_a2", "target": "T2", "evalue": 0.001, "counts_as_hit": False}],
                           {"gene_a1": "HIT", "gene_a2": "IN_PHYLOGENY_ONLY"})   # gene_a3 missing -> NOT_RUN
    mo.close({"tool": {"name": "blastp", "version": "2.x"}}, keep_work=False, quiet=True)
check("tally", tally["HIT"] == 1 and tally["IN_PHYLOGENY_ONLY"] == 1 and tally["NOT_RUN"] == 1 and tally["passed"] == 1, tally)
rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "gene_status.tsv"))]
check("gene_status.tsv", rows[0] == h.STATUS_COLUMNS and rows[1] == ["SpeciesA", "gene_a1", "HIT", "1"]
      and rows[3] == ["SpeciesA", "gene_a3", "NOT_RUN", "0"], rows)
hits = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "hits.tsv"))]
check("hits.tsv", hits[0] == h.HITS_COLUMNS and hits[1][:2] == ["SpeciesA", "gene_a1"] and hits[1][10] == "yes" and hits[2][10] == "no", hits)
passed = open(os.path.join(out, "passed", "SpeciesA.faa")).read()
check("passed fasta keeps original header, only passers", passed.startswith(">gene_a2__SpeciesA\n") and "gene_a1" not in passed and "gene_a3" not in passed, passed)
check("passed_ids.tsv", open(os.path.join(out, "passed_ids.tsv")).read() == "species\tgene_id\nSpeciesA\tgene_a2\n")
rj = json.load(open(os.path.join(out, "run.json")))
check("run.json", rj["label"] == "organelle" and rj["module"] == "homology-blast" and rj["tool"]["name"] == "blastp" and rj["counts"]["SpeciesA"]["input"] == 3, rj)
check("work removed", not os.path.exists(os.path.join(out, "work")))

# --- end-of-run citation note
def run_close(name, quiet, tools):
    d = os.path.join(tmp, "out", name)
    err = io.StringIO()
    with h.ModuleOutput(d, "homology-blast") as m:
        m.add_species("SpeciesA", fa, [], {"gene_a1": "NO_HIT", "gene_a2": "NO_HIT", "gene_a3": "HIT"})
        with contextlib.redirect_stderr(err):
            m.close({"tools": tools} if tools is not None else {}, quiet=quiet)
    return err.getvalue(), open(os.path.join(d, "citations.txt")).read()

printed, saved = run_close("cite1", False, [{"name": "BLAST+ blastp", "version": "2.17.0"}])
check("note says how many passed", "2 of 3 genes passed" in printed, printed)
check("note names the tool and version", "BLAST+ blastp 2.17.0" in printed and "homolog-exclusion-pipeline" in printed, printed)
check("note points to docs/citing.md", "docs/citing.md" in printed)
check("citations.txt equals the printed note", printed.strip() == saved.strip(), (printed, saved))
printed, saved = run_close("cite2", True, [{"name": "BLAST+ blastp", "version": "2.17.0"}])
check("quiet prints nothing but still writes citations.txt", printed == "" and "BLAST+ blastp 2.17.0" in saved, (printed, saved))
printed, saved = run_close("cite3", False, None)
check("note works when no tools were recorded", "This run used: homolog-exclusion-pipeline" in printed, printed)
note_nr = h.citation_note("m", "o", [], {"A": {"input": 5, "passed": 2, "NOT_RUN": 2}, "B": {"input": 3, "passed": 1, "NOT_RUN": 1}})
note_ok = h.citation_note("m", "o", [], {"A": {"input": 5, "passed": 2, "NOT_RUN": 0}})
check("the end-of-run note warns, in plain words, when genes are NOT_RUN (3 here) and says they do not pass",
      "3 gene(s) are NOT_RUN" in note_nr and "do NOT pass" in note_nr, note_nr)
check("the note has no warning line when nothing is NOT_RUN", "WARNING" not in note_ok, note_ok)
printed, saved = run_close("cite4", False, [{"name": "HMMER jackhmmer"}])
check("missing tool version is stated, not invented", "HMMER jackhmmer (version unknown)" in printed, printed)

# --- summary.tsv and output-folder safety
srows = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(out, "summary.tsv"))]
check("summary.tsv: header and per-species numbers",
      srows[0] == ["species", "input", "passed", "HIT", "NO_HIT", "EXCLUDED_ONLY", "SPURIOUS", "IN_PHYLOGENY_ONLY", "NOT_RUN"]
      and srows[1] == ["SpeciesA", "3", "1", "1", "0", "0", "0", "1", "1"], srows)
check("summary.tsv: an ALL row totals the species", srows[-1] == ["ALL"] + srows[1][1:], srows)
blanks = [(name, n) for name in ("hits.tsv", "gene_status.tsv", "summary.tsv", "passed_ids.tsv")
          for n, line in enumerate(open(os.path.join(out, name)).read().splitlines(), 1) if "" in line.split("\t")]
check("no output table has an empty cell: missing values are written as NA (a rule for scripts, R and pandas)", blanks == [], blanks)

d = os.path.join(tmp, "folder_rules")
with h.ModuleOutput(d, "homology-blast", "runA") as m:
    m.close({}, quiet=True)
check("same folder, same run, again: refused", raises(h.HepError, h.ModuleOutput, d, "homology-blast", "runA") is True)
with h.ModuleOutput(d, "homology-blast", "runA", True) as m:   # resume
    pass
check("--resume is allowed into a folder of the same run", True)
with h.ModuleOutput(d, "homology-blast", "runA", False, True) as m:   # force
    pass
check("--force is allowed into a folder of the same run", True)
check("a folder of another module is refused, even with force",
      raises(h.HepError, h.ModuleOutput, d, "domain-cdd", None, False, True) is True)
check("the same module with another label is refused, even with force",
      raises(h.HepError, h.ModuleOutput, d, "homology-blast", "runB", False, True) is True)
stray = w("stray/notes.txt", "not ours")
check("a folder that already holds other files is refused, even with force",
      raises(h.HepError, h.ModuleOutput, os.path.dirname(stray), "domain-cdd", None, False, True) is True)
check("...and its files were not touched", open(stray).read() == "not ours" and os.listdir(os.path.dirname(stray)) == ["notes.txt"])
ab = os.path.join(tmp, "abort_new")
mo_a = h.ModuleOutput(ab, "homology-blast", "x")
mo_a.abort()
check("abort() removes a folder that this run created, so no header-only tables are left behind", not os.path.exists(ab))
ab2 = os.path.join(tmp, "abort_old")
with h.ModuleOutput(ab2, "homology-blast", "x") as m:
    m.close({}, quiet=True)
mo_b = h.ModuleOutput(ab2, "homology-blast", "x", resume=True)
mo_b.abort()
check("abort() leaves alone a folder that already held an earlier run", os.path.isfile(os.path.join(ab2, "run.json")))
os.makedirs(os.path.join(tmp, "empty_ok"))
with h.ModuleOutput(os.path.join(tmp, "empty_ok"), "domain-cdd") as m:
    pass
check("an existing empty folder is fine", True)

# --- runner
rc = h.run_command(["definitely_not_a_program_xyz"], os.path.join(tmp, "l.log"))
check("missing program -> 127", rc == 127)
check("run_command ok", h.run_command([sys.executable, "-c", "print('hi')"], os.path.join(tmp, "l2.log")) == 0
      and "hi" in open(os.path.join(tmp, "l2.log")).read())
check("run_command nonzero", h.run_command([sys.executable, "-c", "import sys; sys.exit(3)"]) == 3)
calls = []
def job(n):
    calls.append(n)
    if n == "bad":
        raise h.CommandError("boom")
markers = os.path.join(tmp, "markers")
res = h.run_jobs(["a", "b", "bad"], job, workers=3, marker_dir=markers)
check("run_jobs failure reported", res["done"] == ["a", "b"] and list(res["failed"]) == ["bad"], res)
calls[:] = []
res = h.run_jobs(["a", "b", "bad"], job, workers=2, marker_dir=markers, resume=True)
check("resume skips done, retries failed", res["skipped"] == ["a", "b"] and calls == ["bad"], (res, calls))

# --- start_module / CLI
ctx = h.start_module("domain-cdd", "test", ["-i", os.path.join(tmp, "in"), "-o", os.path.join(tmp, "o"),
                                             "--evalue", "1e-4", "--set", "min_qcov=25", "--label", "cdd2"])
check("start_module overrides", ctx["cfg"]["evalue"] == 1e-4 and ctx["cfg"]["min_qcov"] == 25 and ctx["args"].label == "cdd2")
check("start_module non_default", set(ctx["non_default"]) == {"evalue", "min_qcov"}, ctx["non_default"])
check("start_module inputs", [s for s, _ in ctx["inputs"]] == ["SpeciesA", "SpeciesB"])
info = h.run_info_base(ctx)
check("run_info_base", info["settings"]["evalue"] == 1e-4 and len(info["inputs"]) == 2)
buf = io.StringIO()
try:
    with contextlib.redirect_stdout(buf):
        h.start_module("domain-cdd", "t", ["-i", os.path.join(tmp, "in"), "-o", "x", "--print-config"])
except SystemExit as e:
    check("print-config exits 0", e.code == 0)
check("print-config output", "evalue=0.01" in buf.getvalue() and "min_domain_cov=50" in buf.getvalue(), buf.getvalue())
try:
    with contextlib.redirect_stderr(io.StringIO()):
        h.start_module("domain-cdd", "t", ["-i", os.path.join(tmp, "in"), "-o", "x", "--set", "nope=1"])
    check("bad --set exits", False)
except SystemExit as e:
    check("bad --set exits nonzero", e.code == 1)

# --- what counts as a species tag after '__'
check("looks_like_species: real names pass (Homo_sapiens, Ceratodon_purpureus_R40, O_sativa, Athaliana)",
      all(h.looks_like_species(t) for t in ("Homo_sapiens", "Ceratodon_purpureus_R40", "O_sativa", "Athaliana",
                                             "Marchantia_polymorpha_subsp._ruderalis")))
check("looks_like_species: random tails and empty do not (FwSs-9Ag, W3i5P, xh_7, u, empty, None)",
      not any(h.looks_like_species(t) for t in ("FwSs-9Ag", "W3i5P", "xh_7", "u", "", None, "iKiFYXF0334")))

# --- every command answers --version with the pipeline version
import subprocess
for script in ("domain_interpro", "domain_cdd", "homology_blast", "homology_jackhmmer", "make_homology_db", "combine", "db_report"):
    pv = subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "bin", script + ".py"),
                         "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    check("%s --version prints the pipeline version" % script, pv.returncode == 0 and pv.stdout.strip() == "%s %s" % (h.PIPELINE_NAME, h.__version__), pv.stdout + pv.stderr)

# --- repository hygiene: nothing personal or machine-specific may ship
import re
REPO = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
LEAK = re.compile(r"/home/(?!user\b)[A-Za-z0-9_.-]+|C:\\Users\\|/Users/(?!user\b)[A-Za-z0-9_.-]+|/mnt/c/Users/(?!user\b)[A-Za-z0-9_.-]+|"
                  r"[A-Za-z0-9._-]+@[A-Za-z0-9-]+\.(?:com|edu|org|net)")
leaks, pycache = [], []
for root, dirs, files in os.walk(REPO):
    dirs[:] = [d for d in dirs if d not in (".git",)]
    if "__pycache__" in dirs:
        dirs.remove("__pycache__")            # created by running the tests; .gitignore keeps it out of the repository
    for f in files:
        if f.endswith((".pyc", ".faa", ".fa")) and os.path.getsize(os.path.join(root, f)) > 2_000_000:
            continue
        path = os.path.join(root, f)
        try:
            text = open(path, encoding="utf-8").read()
        except (UnicodeDecodeError, OSError):
            continue
        for m in LEAK.finditer(text):
            if os.path.basename(path) == "test_core.py":
                continue                        # this file holds the pattern itself
            leaks.append("%s: %s" % (os.path.relpath(path, REPO), m.group(0)))
check("no file in the repository holds a local user path or an e-mail address (example data included)", not leaks, leaks[:5])
check(".gitignore keeps __pycache__ and compiled files out of the repository",
      os.path.isfile(os.path.join(REPO, ".gitignore"))
      and "__pycache__" in open(os.path.join(REPO, ".gitignore")).read())

print()
print("FAILED: %d" % len(fails) if fails else "ALL PASSED")
sys.exit(1 if fails else 0)
