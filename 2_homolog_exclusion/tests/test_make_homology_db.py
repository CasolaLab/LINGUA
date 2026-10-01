"""
test_make_homology_db.py - checks for bin/make_homology_db.py (the `make-db` command). Run: python3 tests/test_make_homology_db.py
Small hand-made FASTA files; the expected answers were worked out by hand.
  Physcomitrium_patens.faa   2 sequences (one ends in '*'): the analysis species
  Ceratodon_purpureus.faa    2 sequences
  Marchantia_polymorpha_subsp._ruderalis.faa  1 sequence (a subspecies of the analysis species below in test 3)
  Blasia_pusilla.faa         1 sequence whose ID already carries __Blasia_pusilla (Stage 1 tag)
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
SCRIPT = os.path.join(ROOT, "bin", "make_homology_db.py")
fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + ((" | " + str(extra)) if (extra and not cond) else ""))
    if not cond:
        fails.append(name)


tmp = tempfile.mkdtemp()
src = os.path.join(tmp, "src")
os.makedirs(src)


def w(name, text, d=src):
    with open(os.path.join(d, name), "w") as fh:
        fh.write(text)


w("Physcomitrium_patens.faa", ">p1 first\nMKV*\n>p2\nMAAA\n")
w("Ceratodon_purpureus.faa", ">c1\nMKKK\n>c2 desc here\nMLLL\n")
w("Marchantia_polymorpha_subsp._ruderalis.faa", ">m1\nMGGG\n")
w("Blasia_pusilla.faa", ">b1__Blasia_pusilla\nMSSS\n")
n = [0]


def run(args):
    n[0] += 1
    out = os.path.join(tmp, "db%d.faa" % n[0])
    p = subprocess.run([sys.executable, SCRIPT, "-o", out, "--quiet"] + args, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
    return p, out


def ids(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


def rows(path):
    ls = [l.rstrip("\n").split("\t") for l in open(path + ".species.tsv")]
    return [dict(zip(ls[0], r)) for r in ls[1:]]


p, out = run(["-i", src])
check("builds", p.returncode == 0, p.stderr)
check("every ID is tagged, existing tags kept",
      ids(out) == ["b1__Blasia_pusilla", "c1__Ceratodon_purpureus", "c2__Ceratodon_purpureus",
                   "m1__Marchantia_polymorpha_subsp._ruderalis", "p1__Physcomitrium_patens", "p2__Physcomitrium_patens"]
      or sorted(ids(out)) == sorted(["b1__Blasia_pusilla", "c1__Ceratodon_purpureus", "c2__Ceratodon_purpureus",
                                     "m1__Marchantia_polymorpha_subsp._ruderalis", "p1__Physcomitrium_patens",
                                     "p2__Physcomitrium_patens"]), ids(out))
check("the description after the ID is kept", "c2__Ceratodon_purpureus desc here" in open(out).read())
check("'*' is stripped", "*" not in "".join(l for l in open(out) if not l.startswith(">")))
t = rows(out)
check("species table: 4 species, all included, no empty cells",
      len(t) == 4 and all(r["status"] == "INCLUDED" for r in t)
      and all(v != "" for r in t for v in r.values()), t)
check("table counts: patens 2, Blasia 1 already tagged",
      {r["species"]: (r["n_sequences"], r["n_ids_already_tagged"]) for r in t}["Blasia_pusilla"] == ("1", "1"))

p, out = run(["-i", src])
p2 = subprocess.run([sys.executable, SCRIPT, "-i", src, "-o", out, "--quiet"], stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, universal_newlines=True)
check("refuses to overwrite without --force", p2.returncode != 0 and "--force" in p2.stderr, p2.stderr)
p2 = subprocess.run([sys.executable, SCRIPT, "-i", src, "-o", out, "--quiet", "--force"], stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, universal_newlines=True)
check("--force overwrites", p2.returncode == 0, p2.stderr)

p, out = run(["-i", src, "--exclude-species", "Physcomitrium_patens"])
check("--exclude-species leaves the species out and lists it as EXCLUDED",
      p.returncode == 0 and not any("Physcomitrium" in i for i in ids(out))
      and {r["species"]: r["status"] for r in rows(out)}["Physcomitrium_patens"].startswith("EXCLUDED"), p.stderr)

ana = os.path.join(tmp, "ana")
os.makedirs(ana)
w("Marchantia_polymorpha_final.faa", ">q1\nMAAA\n", ana)
p, out = run(["-i", src, "--exclude-inputs", ana])
check("--exclude-inputs: a subspecies file is the same species as the analysis species",
      p.returncode == 0 and not any("Marchantia" in i for i in ids(out)), p.stderr)
ana2 = os.path.join(tmp, "ana2")
os.makedirs(ana2)
w("Mp1.faa", ">q1\nMAAA\n", ana2)          # an analysis file whose name is a code, not a species name
smap = os.path.join(tmp, "smap.tsv")
open(smap, "w").write("Species\tBasename\nMarchantia polymorpha\tMp1\n")
p, out = run(["-i", src, "--exclude-inputs", ana2])
check("--exclude-inputs with a coded file name (Mp1) and no map: the species is 'Mp1', so nothing is excluded",
      p.returncode == 0 and any("Marchantia" in i for i in ids(out)), p.stderr)
p, out = run(["-i", src, "--exclude-inputs", ana2, "--species-map", smap])
check("--species-map also names the analysis species for --exclude-inputs (Mp1 = Marchantia polymorpha, so its subspecies is left out)",
      p.returncode == 0 and not any("Marchantia" in i for i in ids(out)), p.stderr)
p, out = run(["-i", src, "--exclude-inputs", ana, "--match-level", "exact"])
check("--match-level exact keeps the subspecies", any("Marchantia" in i for i in ids(out)), p.stderr)

d = os.path.join(tmp, "dup")
os.makedirs(d)
w("A_one.faa", ">x\nMAA\n>x\nMCC\n", d)
p, out = run(["-i", d])
check("duplicate IDs are refused and no output is left",
      p.returncode != 0 and "duplicate" in p.stderr and not os.path.exists(out) and not os.path.exists(out + ".partial"), p.stderr)
d2 = os.path.join(tmp, "empty")
os.makedirs(d2)
w("B_two.faa", ">y\n", d2)
p, out = run(["-i", d2])
check("a record with no sequence is refused", p.returncode != 0 and not os.path.exists(out), p.stderr)
p, out = run(["-i", src, "--exclude-species", "Physcomitrium_patens", "Ceratodon_purpureus",
              "Marchantia_polymorpha", "Blasia_pusilla"])
check("excluding everything is refused", p.returncode != 0 and "excluded" in p.stderr, p.stderr)
p, out = run(["-i", os.path.join(tmp, "nowhere")])
check("missing input is reported", p.returncode != 0 and "not found" in p.stderr, p.stderr)

# an '__' inside a random-looking ID is not a species tag: make-db adds its own tag and keeps a real one
oddd = os.path.join(tmp, "oddtags")
os.makedirs(oddd)
w("Bacillus_subtilis.faa", ">Acid_ENSB:xy__FwSs-9Ag first\nMKV\n>c1__Bacillus_subtilis\nMKV\n>plain\nMKV\n", oddd)
p, out = run(["-i", oddd])
oids = ids(out)
check("a '__' inside a random-looking ID is NOT taken for a species tag: our own tag is added after it",
      p.returncode == 0 and "Acid_ENSB:xy__FwSs-9Ag__Bacillus_subtilis" in oids and "FwSs-9Ag" not in [i.rpartition("__")[2] for i in oids], oids)
check("a real tag is kept, and an untagged ID gets one", "c1__Bacillus_subtilis" in oids and "plain__Bacillus_subtilis" in oids, oids)
row = rows(out)[0]
check("the species table counts them apart: 1 already tagged, 1 retagged (random '__'), 3 in all",
      (row["n_sequences"], row["n_ids_already_tagged"], row["n_ids_retagged"]) == ("3", "1", "1"), row)
check("species read back from the new IDs is the real species for every ID",
      {i.rpartition("__")[2] for i in oids} == {"Bacillus_subtilis"}, oids)

print()
print("FAILED: %d" % len(fails) if fails else "ALL PASSED")
sys.exit(1 if fails else 0)
