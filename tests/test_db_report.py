"""
test_db_report.py - checks for bin/db_report.py. Run from anywhere:   python3 tests/test_db_report.py
No download: the test writes a small taxonomy in the exact NCBI names.dmp / nodes.dmp format (fields separated by TAB|TAB,
lines ending TAB|), with the same rank names NCBI uses today (domain, kingdom, class, order, family, genus, species, ...).
Expected answers were worked out by hand from this mini tree:

  root(1) - cellular organisms(131567) - Bacteria(2, domain) - Escherichia(561) - Escherichia coli(562)
                                       - Eukaryota(2759, domain) - Fungi(4751, kingdom) - Saccharomyces(4930) - S. cerevisiae(4932)
                                                              - Metazoa(33208, kingdom) - Homo(9605) - Homo sapiens(9606)
                                                              - Viridiplantae(33090, kingdom) - Magnoliopsida(3398, class)
                                                                   - Brassicales(3699, order) - Brassicaceae(3700, family)
                                                                          - Arabidopsis(3701) - A. thaliana(3702)
                                                                          - Brassica(9000) - B. rapa(9001)
                                                                        - Cleomaceae(28531) - Tarenaya(28532) - T. hassleriana(28533)
                                                                   - Poales(4527, order) - Oryza(4528) - Oryza sativa(4530)
                                                              - Ambiguus(7001, genus) under Fungi AND Ambiguus(7002, genus) under Metazoa
Clade = Brassicaceae (family). Its ladder at the default ranks: Brassicales (order), Magnoliopsida (class),
Viridiplantae (kingdom), Eukaryota (domain), then "outside all of these".
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
SCRIPT = os.path.join(ROOT, "bin", "db_report.py")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import db_report as d  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + ((" | " + str(extra)) if (extra and not cond) else ""))
    if not cond:
        fails.append(name)


tmp = tempfile.mkdtemp()
tax = os.path.join(tmp, "tax")
os.makedirs(tax)
NODES = [(1, 1, "no rank"), (131567, 1, "no rank"), (2, 131567, "domain"), (561, 2, "genus"), (562, 561, "species"),
         (2759, 131567, "domain"), (4751, 2759, "kingdom"), (4930, 4751, "genus"), (4932, 4930, "species"),
         (33208, 2759, "kingdom"), (9605, 33208, "genus"), (9606, 9605, "species"),
         (33090, 2759, "kingdom"), (3398, 33090, "class"), (3699, 3398, "order"), (3700, 3699, "family"),
         (3701, 3700, "genus"), (3702, 3701, "species"), (9000, 3700, "genus"), (9001, 9000, "species"),
         (28531, 3699, "family"), (28532, 28531, "genus"), (28533, 28532, "species"),
         (4527, 3398, "order"), (4528, 4527, "genus"), (4530, 4528, "species"),
         (7001, 4751, "genus"), (7002, 33208, "genus")]
NAMES = [(1, "root", "scientific name"), (131567, "cellular organisms", "scientific name"), (2, "Bacteria", "scientific name"),
         (561, "Escherichia", "scientific name"), (562, "Escherichia coli", "scientific name"),
         (2759, "Eukaryota", "scientific name"), (4751, "Fungi", "scientific name"), (4930, "Saccharomyces", "scientific name"),
         (4932, "Saccharomyces cerevisiae", "scientific name"), (33208, "Metazoa", "scientific name"),
         (9605, "Homo", "scientific name"), (9606, "Homo sapiens", "scientific name"), (9606, "human", "genbank common name"),
         (9606, "Homo sapiens Linnaeus, 1758", "authority"),
         (33090, "Viridiplantae", "scientific name"), (3398, "Magnoliopsida", "scientific name"),
         (3699, "Brassicales", "scientific name"), (3700, "Brassicaceae", "scientific name"), (3701, "Arabidopsis", "scientific name"),
         (3702, "Arabidopsis thaliana", "scientific name"), (3702, "thale cress", "genbank common name"),
         (9000, "Brassica", "scientific name"), (9001, "Brassica rapa", "scientific name"),
         (28531, "Cleomaceae", "scientific name"), (28532, "Tarenaya", "scientific name"),
         (28533, "Tarenaya hassleriana", "scientific name"), (28533, "Cleome hassleriana", "synonym"),
         (4527, "Poales", "scientific name"), (4528, "Oryza", "scientific name"), (4530, "Oryza sativa", "scientific name"),
         (7001, "Ambiguus", "scientific name"), (7002, "Ambiguus", "scientific name")]
with open(os.path.join(tax, "nodes.dmp"), "w") as fh:
    for t, par, rank in NODES:
        fh.write("\t|\t".join([str(t), str(par), rank, "", "0", "0", "1", "0", "0", "0", "0", "0", ""]) + "\t|\n")
with open(os.path.join(tax, "names.dmp"), "w") as fh:
    for t, name, cls in NAMES:
        fh.write("\t|\t".join([str(t), name, "", cls]) + "\t|\n")
open(os.path.join(tax, "TAXONOMY_SOURCE.txt"), "w").write("source=test\narchive_last_modified=2026-09-26\narchive_md5=abc\n")

# database: species, number of sequences
DB = [("Arabidopsis_thaliana", 3), ("Brassica_rapa_FPsc", 2), ("Tarenaya_hassleriana", 4), ("Oryza_sativa", 5),
      ("Saccharomyces_cerevisiae", 30), ("Homo_sapiens", 6), ("Escherichia_coli", 5), ("O_sativa", 1), ("Ambiguus_foo", 1)]
fa = os.path.join(tmp, "db.faa")
with open(fa, "w") as fh:
    for sp, n in DB:
        for i in range(n):
            fh.write(">g%d__%s\nMKV\n" % (i, sp))
    fh.write(">plain_untagged\nMKV\n>Acid_ENSB:xy__FwSs-9Ag\nMKV\n")     # no tag; '__' inside a random ID
counter = [0]


def run(args, ok=True):
    counter[0] += 1
    out = os.path.join(tmp, "rep%d" % counter[0])
    p = subprocess.run([sys.executable, SCRIPT, "-o", out, "--quiet"] + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       universal_newlines=True)
    return p, out


def table(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return [dict(zip(rows[0], r)) for r in rows[1:]]


p, out = run(["--db", fa, "--clade", "Brassicaceae", "--taxonomy", tax, "--min-seqs", "2", "--thin", "2"])
check("db-report succeeds", p.returncode == 0, p.stderr)
st = {r["species"]: r for r in table(os.path.join(out, "species_tiers.tsv"))}
tier_of = {sp: r["tier"] for sp, r in st.items()}
check("Arabidopsis thaliana and Brassica rapa (strain FPsc) are INSIDE the analysis clade",
      tier_of["Arabidopsis_thaliana"] == "INSIDE_CLADE" and tier_of["Brassica_rapa_FPsc"] == "INSIDE_CLADE", tier_of)
check("Tarenaya hassleriana: nearest outside tier = Brassicales (order)", tier_of["Tarenaya_hassleriana"] == "Brassicales (order)", tier_of)
check("Oryza sativa: Magnoliopsida (class)", tier_of["Oryza_sativa"] == "Magnoliopsida (class)")
check("yeast and human: Eukaryota (domain), because neither is a plant", tier_of["Saccharomyces_cerevisiae"] == "Eukaryota (domain)"
      and tier_of["Homo_sapiens"] == "Eukaryota (domain)")
check("E. coli: outside all of these", tier_of["Escherichia_coli"] == "outside all of these (rest of the tree)")
check("an abbreviated name the taxonomy cannot match is UNMATCHED with NA (never blank); an ambiguous genus is UNMATCHED too",
      tier_of["O_sativa"] == "UNMATCHED" and tier_of["Ambiguus_foo"] == "UNMATCHED" and st["O_sativa"]["taxid"] == "NA", tier_of)
check("a strain is matched through its species: Brassica_rapa_FPsc found as 'brassica rapa' at species level",
      st["Brassica_rapa_FPsc"]["match_level"] == "species" and st["Brassica_rapa_FPsc"]["taxid"] == "9001", st["Brassica_rapa_FPsc"])
check("common names and authority strings are not used for matching: 'thale cress' is not needed, exact scientific name wins",
      st["Arabidopsis_thaliana"]["match_level"] == "exact" and st["Arabidopsis_thaliana"]["taxid"] == "3702")
tiers = {r["tier"]: r for r in table(os.path.join(out, "tier_summary.tsv"))}
check("tier summary: Brassicales 1 species (4 seqs), Magnoliopsida 1 (5), Viridiplantae 0 = EMPTY, Eukaryota 2 (36), outside 1 (5)",
      (tiers["Brassicales (order)"]["n_species"], tiers["Brassicales (order)"]["n_sequences"]) == ("1", "4")
      and (tiers["Magnoliopsida (class)"]["n_species"], tiers["Magnoliopsida (class)"]["n_sequences"]) == ("1", "5")
      and (tiers["Viridiplantae (kingdom)"]["n_species"], tiers["Viridiplantae (kingdom)"]["status"]) == ("0", "EMPTY")
      and (tiers["Eukaryota (domain)"]["n_species"], tiers["Eukaryota (domain)"]["n_sequences"]) == ("2", "36")
      and tiers["outside all of these (rest of the tree)"]["n_species"] == "1", tiers)
check("a tier below --thin is THIN (Brassicales has 1 species, --thin 2); the inside tier is flagged INSIDE_CLADE",
      tiers["Brassicales (order)"]["status"] == "THIN" and tiers["INSIDE_CLADE"]["status"] == "INSIDE_CLADE", tiers)
check("tier order: inside first, then nearest outside group outward, then the rest, then unmatched",
      list(tiers) == ["INSIDE_CLADE", "Brassicales (order)", "Magnoliopsida (class)", "Viridiplantae (kingdom)",
                      "Eukaryota (domain)", "outside all of these (rest of the tree)", "UNMATCHED"], list(tiers))
check("yeast holds 30 of 57 tagged sequences = 52.6%: flagged DOMINANT", "DOMINANT" in st["Saccharomyces_cerevisiae"]["flags"]
      and st["Saccharomyces_cerevisiae"]["pct_of_database"] == "52.6", st["Saccharomyces_cerevisiae"])
check("a species below --min-seqs is FEW_SEQUENCES (Brassica rapa has 2, min 2 is fine; O_sativa 1 is below)",
      "FEW_SEQUENCES" in st["O_sativa"]["flags"] and "FEW_SEQUENCES" not in st["Brassica_rapa_FPsc"]["flags"], st["O_sativa"])
rep = json.load(open(os.path.join(out, "report.json")))
check("report.json records the clade, the taxonomy copy used and the warnings",
      "Brassicaceae" in rep["clade"] and rep["taxonomy"]["archive_last_modified"] == "2026-09-26"
      and any("INSIDE_CLADE" in w for w in rep["warnings"]) and any("EMPTY" in w for w in rep["warnings"]), rep["warnings"])
check("untagged IDs are counted and the accidental '__' inside a random ID is NOT read as a species",
      any("2 sequence IDs carry no usable species tag" in w and "random-looking" in w for w in rep["warnings"])
      and "FwSs-9Ag" not in st, rep["warnings"])
check("no table has an empty cell", all(v != "" for f in ("species_tiers.tsv", "tier_summary.tsv") for r in table(os.path.join(out, f)) for v in r.values()))

# --- the clade can also be given other ways
inp = os.path.join(tmp, "inputs")
os.makedirs(inp)
open(os.path.join(inp, "Arabidopsis_thaliana_final.faa"), "w").write(">a\nMK\n")
open(os.path.join(inp, "Brassica_rapa.faa"), "w").write(">b\nMK\n")
p, out2 = run(["--db", fa, "--clade-inputs", inp, "--taxonomy", tax])
check("--clade-inputs: the clade is the common ancestor of the analysis species (Brassicaceae), same tiers",
      p.returncode == 0 and "Brassicaceae" in json.load(open(os.path.join(out2, "report.json")))["clade"], p.stderr)
lst = os.path.join(tmp, "an.txt")
open(lst, "w").write("Arabidopsis thaliana\nTarenaya hassleriana\n")
p, out3 = run(["--db", fa, "--clade-species", lst, "--taxonomy", tax])
check("--clade-species: the analysis species include Tarenaya, so the clade is Brassicales and Tarenaya is INSIDE it",
      p.returncode == 0 and {r["species"]: r["tier"] for r in table(os.path.join(out3, "species_tiers.tsv"))}["Tarenaya_hassleriana"] == "INSIDE_CLADE", p.stderr)
p, out4 = run(["--db", fa, "--clade-taxid", "3700", "--taxonomy", tax])
check("--clade-taxid works", p.returncode == 0, p.stderr)
p, out5 = run(["--db", fa, "--clade", "Ambiguus", "--taxonomy", tax])
check("an ambiguous clade name is refused, not guessed", p.returncode != 0 and "matches 2 taxa" in p.stderr, p.stderr)
p, out5 = run(["--db", fa, "--clade", "Nothingus", "--taxonomy", tax])
check("an unknown clade name is refused", p.returncode != 0 and "not found" in p.stderr, p.stderr)
tm = os.path.join(tmp, "tm.tsv")
open(tm, "w").write("O_sativa\t4530\n")
p, out6 = run(["--db", fa, "--clade", "Brassicaceae", "--taxonomy", tax, "--taxid-map", tm])
check("--taxid-map fixes a name the taxonomy cannot match (O_sativa -> Oryza sativa, Magnoliopsida)",
      {r["species"]: r["tier"] for r in table(os.path.join(out6, "species_tiers.tsv"))}["O_sativa"] == "Magnoliopsida (class)", p.stderr)
p, out7 = run(["--db", fa, "--clade", "Brassicaceae", "--taxonomy", tax, "--ranks", "all"])
tiers7 = [r["tier"] for r in table(os.path.join(out7, "tier_summary.tsv"))]
check("--ranks all uses every ancestor (adds Eukaryota's 'cellular organisms' and the root)", p.returncode == 0 and len(tiers7) > 6, tiers7)

# --- tree
tr = os.path.join(tmp, "t.nwk")
p, out8 = run(["--db", fa, "--clade", "Brassicaceae", "--taxonomy", tax, "--tree", tr])
nw = open(tr).read().strip()
check("tree: valid Newick (balanced, ends with ;), with each matched species once and the analysis clade marked",
      p.returncode == 0 and nw.endswith(";") and nw.count("(") == nw.count(")") and nw.count("ANALYSIS_CLADE_Brassicaceae") == 1
      and all(sp in nw for sp in ("Arabidopsis_thaliana", "Tarenaya_hassleriana", "Escherichia_coli", "Homo_sapiens"))
      and "O_sativa" not in nw, nw)
check("tree: tips = 6 matched species + 2 more (Brassica_rapa_FPsc, Oryza_sativa, Saccharomyces) = 7 species + the clade marker",
      sum(nw.count(s) for s in ("Arabidopsis_thaliana", "Brassica_rapa_FPsc", "Tarenaya_hassleriana", "Oryza_sativa",
                                "Saccharomyces_cerevisiae", "Homo_sapiens", "Escherichia_coli")) == 7, nw)

# --- your own tiers table, no taxonomy
tt = os.path.join(tmp, "tiers.tsv")
open(tt, "w").write("Tarenaya_hassleriana\tclose\nOryza_sativa\tplants\nSaccharomyces_cerevisiae\tfungi\nHomo_sapiens\tanimals\n")
lst2 = os.path.join(tmp, "names.txt")
open(lst2, "w").write("Tarenaya_hassleriana\t4000\nOryza_sativa\t5000\nSaccharomyces_cerevisiae\t6000\nHomo_sapiens\t7000\nMystery_sp\n")
p, out9 = run(["--species-list", lst2, "--tiers", tt, "--expected-tiers", "close,plants,protists,fungi,animals"])
t9 = {r["tier"]: r for r in table(os.path.join(out9, "tier_summary.tsv"))}
check("--tiers with --species-list: no taxonomy needed; an expected tier with no species is EMPTY; an unlisted species is UNMATCHED",
      p.returncode == 0 and t9["protists"]["status"] == "EMPTY" and t9["UNMATCHED"]["n_species"] == "1" and t9["fungi"]["n_sequences"] == "6000", (p.stderr, t9))
check("a species with no count gives NA (never 0) in the sequence columns",
      {r["species"]: r for r in table(os.path.join(out9, "species_tiers.tsv"))}["Mystery_sp"]["n_sequences"] == "NA")
p, _ = run(["--species-list", lst2, "--tiers", tt, "--tree", tr + "2"])
check("--tree needs the taxonomy", p.returncode != 0 and "taxonomy" in p.stderr, p.stderr)

# --- common names and authority strings are not species names
cn = os.path.join(tmp, "common.txt")
open(cn, "w").write("\n".join(["Human", "Thale_cress", "Cleome_hassleriana"]) + "\n")
p, outc = run(["--species-list", cn, "--clade", "Brassicaceae", "--taxonomy", tax])
tc = {r["species"]: r["tier"] for r in table(os.path.join(outc, "species_tiers.tsv"))}
check("a common name (human, thale cress) is NOT accepted as a species name: UNMATCHED, never guessed; a real synonym is (Cleome hassleriana)",
      p.returncode == 0 and tc["Human"] == "UNMATCHED" and tc["Thale_cress"] == "UNMATCHED"
      and tc["Cleome_hassleriana"] == "Brassicales (order)", (p.stderr, tc))

# --- refusals
p, _ = run(["--db", fa, "--clade", "Brassicaceae"])
check("no taxonomy and no tiers table is refused with a pointer to the download script", p.returncode != 0 and "get_taxonomy" in p.stderr, p.stderr)
p, _ = run(["--db", fa, "--clade", "Brassicaceae", "--taxonomy", os.path.join(tmp, "empty")])
check("a taxonomy folder without the dump files is refused with the download command", p.returncode != 0 and "get_taxonomy.sh" in p.stderr, p.stderr)
p, _ = run(["--db", fa, "--taxonomy", tax])
check("a taxonomy without a clade is refused", p.returncode != 0 and "exactly one" in p.stderr, p.stderr)
p, _ = run(["--taxonomy", tax, "--clade", "Brassicaceae"])
check("no database is refused", p.returncode != 0 and "--db" in p.stderr, p.stderr)
untagged = os.path.join(tmp, "untagged.faa")
open(untagged, "w").write(">a\nMK\n>b\nMK\n")
p, _ = run(["--db", untagged, "--clade", "Brassicaceae", "--taxonomy", tax])
check("a database with no species tags is refused with the fix", p.returncode != 0 and "make-db" in p.stderr, p.stderr)
p2 = subprocess.run([sys.executable, SCRIPT, "--db", fa, "--clade", "Brassicaceae", "--taxonomy", tax, "-o", out, "--quiet"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
check("an earlier report folder is not overwritten without --force", p2.returncode != 0 and "already holds results" in p2.stderr, p2.stderr)

# --- pieces
check("species tags: a real name passes, a random tail does not", d.looks_like_species("Homo_sapiens") and d.looks_like_species("Ceratodon_purpureus_R40")
      and d.looks_like_species("Marchantia_polymorpha_subsp._ruderalis") and not d.looks_like_species("FwSs-9Ag") and not d.looks_like_species("u"))
check("the dump is read in NCBI's real layout (TAB|TAB): lineage of Homo sapiens goes up to the root",
      d.Taxonomy(tax).lineage(9606) == [9606, 9605, 33208, 2759, 131567, 1])

print()
print("FAILED: %d" % len(fails) if fails else "ALL PASSED")
sys.exit(1 if fails else 0)
