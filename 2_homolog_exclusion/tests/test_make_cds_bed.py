"""
test_make_cds_bed.py - checks for bin/make_cds_bed.py. Run from anywhere:   python3 tests/test_make_cds_bed.py

Two hand-built mini GFF3 files, in the exact layout of the two real annotation sources this was checked against
(Phytozome-style, clean Name= on the mRNA; and MAKER-style, no isoform suffix, no Name= that matches an ID list
that uses accession-style IDs, forcing the CDS-level fallback match). Expected answers were worked out by hand.

phytozome.gff3 (Phytozome-style, mirrors real Athaliana.gff): one gene AT1G01010 with two mRNA isoforms.
  g1: AT1G01010.1, target list asks for "AT1G01010.1" -> matched via Name= on the mRNA line, two CDS segments
      (60 and 40 bp), one exon-only (non-coding) exon in between (3 exon rows, 2 CDS rows)
  g2: AT1G01010.2, NOT in the target list -> ignored entirely
  a gene AT1G02020 with one mRNA AT1G02020.1, one CDS line duplicated verbatim (a known real-world GFF quirk):
      must count once in the BED and once in the tally, not twice
  a target ID "AT9G99999.1" that matches nothing in the file -> reported in not_found.tsv

maker.gff3 (MAKER-style, mirrors real Aarabicum.gff): mRNA Name= equals the bare gene name, shared by every
  isoform of the gene (no per-isoform suffix) and the target list uses that gene name directly, so this file
  matches through the primary mRNA path, same as the real data.

fallback.gff3: an mRNA line with NO Name= and an ID unrelated to the target list (so the primary path finds
  nothing), whose CDS lines carry their own protein_id= that IS in the target list: must be found only through
  the CDS-level fallback, and (this is the known, documented gap) its exon rows are not counted because exons
  have no fallback path of their own.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
SCRIPT = os.path.join(ROOT, "bin", "make_cds_bed.py")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import make_cds_bed as m  # noqa: E402
import hep_common as h  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + ((" | " + str(extra)) if (extra and not cond) else ""))
    if not cond:
        fails.append(name)


tmp = tempfile.mkdtemp()


def w(path, text):
    with open(path, "w", newline="\n") as fh:
        fh.write(text)


PHYTOZOME = """##gff-version 3
Chr1\tphytozome\tgene\t1\t500\t.\t+\t.\tID=AT1G01010;Name=AT1G01010
Chr1\tphytozome\tmRNA\t1\t500\t.\t+\t.\tID=AT1G01010.1.v1;Name=AT1G01010.1;Parent=AT1G01010
Chr1\tphytozome\texon\t1\t60\t.\t+\t.\tID=e1;Parent=AT1G01010.1.v1
Chr1\tphytozome\tCDS\t1\t60\t.\t+\t0\tID=c1;Parent=AT1G01010.1.v1
Chr1\tphytozome\texon\t100\t120\t.\t+\t.\tID=e2;Parent=AT1G01010.1.v1
Chr1\tphytozome\texon\t200\t239\t.\t+\t.\tID=e3;Parent=AT1G01010.1.v1
Chr1\tphytozome\tCDS\t200\t239\t.\t+\t2\tID=c2;Parent=AT1G01010.1.v1
Chr1\tphytozome\tmRNA\t1\t500\t.\t+\t.\tID=AT1G01010.2.v1;Name=AT1G01010.2;Parent=AT1G01010
Chr1\tphytozome\texon\t1\t80\t.\t+\t.\tID=e4;Parent=AT1G01010.2.v1
Chr1\tphytozome\tCDS\t1\t80\t.\t+\t0\tID=c3;Parent=AT1G01010.2.v1
Chr2\tphytozome\tgene\t1\t300\t.\t+\t.\tID=AT1G02020;Name=AT1G02020
Chr2\tphytozome\tmRNA\t1\t300\t.\t+\t.\tID=AT1G02020.1.v1;Name=AT1G02020.1;Parent=AT1G02020
Chr2\tphytozome\texon\t1\t100\t.\t+\t.\tID=e5;Parent=AT1G02020.1.v1
Chr2\tphytozome\tCDS\t1\t100\t.\t+\t0\tID=c4;Parent=AT1G02020.1.v1
Chr2\tphytozome\tCDS\t1\t100\t.\t+\t0\tID=c4;Parent=AT1G02020.1.v1
"""
w(os.path.join(tmp, "phytozome.gff3"), PHYTOZOME)
w(os.path.join(tmp, "phytozome_ids.tsv"), "Species\tGene_ID\n"
  "Arabidopsis thaliana\tAT1G01010.1\nArabidopsis thaliana\tAT1G02020.1\nArabidopsis thaliana\tAT9G99999.1\n")

MAKER = """##gff-version 3
LG1\tmaker\tgene\t1\t400\t.\t+\t.\tID=Aa1G10
LG1\tmaker\tmRNA\t1\t400\t.\t+\t.\tID=Aa1G10-mRNA-1;Name=Aa1G10;Parent=Aa1G10
LG1\tmaker\texon\t1\t50\t.\t+\t.\tID=e1;Parent=Aa1G10-mRNA-1
LG1\tmaker\tCDS\t1\t50\t.\t+\t0\tID=c1;Parent=Aa1G10-mRNA-1
LG1\tmaker\texon\t120\t170\t.\t+\t.\tID=e2;Parent=Aa1G10-mRNA-1
LG1\tmaker\tCDS\t120\t170\t.\t+\t0\tID=c2;Parent=Aa1G10-mRNA-1
LG2\tmaker\tmRNA\t1\t90\t.\t+\t.\tID=orphanRNA;Name=orphanRNA
LG2\tmaker\tCDS\t1\t90\t.\t+\t0\tID=c3;Parent=orphanRNA
"""
w(os.path.join(tmp, "maker.gff3"), MAKER)
w(os.path.join(tmp, "maker_ids.tsv"), "Species\tGene_ID\nAethionema arabicum\tAa1G10\nAethionema arabicum\torphanRNA\n")

FALLBACK = """##gff-version 3
scaffold1\tsrc\tgene\t1\t200\t.\t+\t.\tID=geneX
scaffold1\tsrc\tmRNA\t1\t200\t.\t+\t.\tID=mrnaX;Parent=geneX
scaffold1\tsrc\tCDS\t1\t90\t.\t+\t0\tID=cdsX1;Parent=mrnaX;protein_id=XP_00001.1
scaffold1\tsrc\texon\t1\t90\t.\t+\t.\tID=exX1;Parent=mrnaX
"""
w(os.path.join(tmp, "fallback.gff3"), FALLBACK)
w(os.path.join(tmp, "fallback_ids.tsv"), "Species\tGene_ID\nSp\tXP_00001.1\n")

n = [0]


def run(args, ok=True):
    n[0] += 1
    out = os.path.join(tmp, "out%d" % n[0])
    p = subprocess.run([sys.executable, SCRIPT, "-o", out, "--quiet"] + args, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
    return p, out


def table(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return [dict(zip(rows[0], r)) for r in rows[1:]]


# --- Phytozome-style: primary mRNA-Name match, duplicate CDS line, an unmatched target
p, out = run(["--ids", os.path.join(tmp, "phytozome_ids.tsv"), "--gff", os.path.join(tmp, "phytozome.gff3"), "--species-name", "Arabidopsis_thaliana"])
check("Phytozome-style file: succeeds (exit 2: one target not found)", p.returncode == 2, p.stderr)
mt = {r["matched_id"]: r for r in table(os.path.join(out, "master.tsv"))}
check("AT1G01010.1: 2 CDS segments (60+40=100bp), 2 exons counted (the exon-only middle exon is not a CDS exon "
      "but IS counted as an exon), duplicate CDS line skipped in AT1G02020.1 not this one",
      mt["AT1G01010.1"]["n_cds_segments"] == "2" and mt["AT1G01010.1"]["cds_length"] == "100"
      and mt["AT1G01010.1"]["n_exons"] == "3", mt["AT1G01010.1"])
check("AT1G01010.2 (not a target) contributes nothing", "AT1G01010.2" not in mt)
check("AT1G02020.1: the identical duplicate CDS line is NOT double-counted (1 segment, 100bp, not 2/200)",
      mt["AT1G02020.1"]["n_cds_segments"] == "1" and mt["AT1G02020.1"]["cds_length"] == "100", mt["AT1G02020.1"])
bed_lines = open(os.path.join(out, "Arabidopsis_thaliana.bed")).read().strip().splitlines()
check("the BED itself also has exactly one row for the duplicated CDS, not two",
      sum(1 for l in bed_lines if "AT1G02020.1" in l) == 1, bed_lines)
check("AT9G99999.1 is reported in not_found.tsv, not silently dropped",
      [r["gene_id"] for r in table(os.path.join(out, "not_found.tsv"))] == ["AT9G99999.1"])
check("gene_length is looked up from the real 'gene' feature span (500bp for AT1G01010's gene row)",
      mt["AT1G01010.1"]["gene_length"] == "500")
check("no table has an empty cell", all(v != "" for f in ("master.tsv", "not_found.tsv", "stats.tsv")
                                        for r in table(os.path.join(out, f)) for v in r.values()))

# --- MAKER-style: matches through the same primary path (Name= on mRNA), no isoform suffix
p, out = run(["--ids", os.path.join(tmp, "maker_ids.tsv"), "--gff", os.path.join(tmp, "maker.gff3"), "--species-name", "Aethionema_arabicum"])
check("MAKER-style file: succeeds, exit 0 (found)", p.returncode == 0, p.stderr)
mt = {r["matched_id"]: r for r in table(os.path.join(out, "master.tsv"))}
check("Aa1G10: 2 CDS segments, 50+51=101bp, 2 exons", mt["Aa1G10"]["n_cds_segments"] == "2"
      and mt["Aa1G10"]["cds_length"] == "101" and mt["Aa1G10"]["n_exons"] == "2", mt["Aa1G10"])
check("orphanRNA has no Parent= at all: gene_id and gene_length are genuinely NA, never a fake 0",
      mt["orphanRNA"]["gene_id"] == "NA" and mt["orphanRNA"]["gene_length"] == "NA", mt["orphanRNA"])

# --- the fallback path itself (no real data ever exercised this)
p, out = run(["--ids", os.path.join(tmp, "fallback_ids.tsv"), "--gff", os.path.join(tmp, "fallback.gff3")])
check("fallback match (CDS protein_id=, no usable mRNA Name=): still found, exit 0", p.returncode == 0, p.stderr)
mt = {r["matched_id"]: r for r in table(os.path.join(out, "master.tsv"))}
check("matched under its protein_id, with the right CDS coordinates", "XP_00001.1" in mt
      and mt["XP_00001.1"]["n_cds_segments"] == "1" and mt["XP_00001.1"]["cds_length"] == "90", mt)
check("documented gap: its exon row is NOT counted, because exon lines have no fallback match of their own",
      mt["XP_00001.1"]["n_exons"] == "0", mt["XP_00001.1"])
check("gene_length is still resolved for a fallback-matched CDS: the parent-of-parent lookup (mRNA -> gene) runs "
      "for every mRNA, matched or not, so geneX's real 200bp span is found even though the match itself came "
      "through the CDS-level fallback, not the mRNA",
      mt["XP_00001.1"]["gene_length"] == "200", mt["XP_00001.1"])

# --- batch mode: two species, matched by the heuristic (no map)
batch_ids = os.path.join(tmp, "batch_ids")
batch_gff = os.path.join(tmp, "batch_gff")
os.makedirs(batch_ids); os.makedirs(batch_gff)
w(os.path.join(batch_ids, "Arabidopsis_thaliana.tsv"), "Species\tGene_ID\nArabidopsis thaliana\tAT1G01010.1\n")
w(os.path.join(batch_ids, "Aethionema_arabicum.tsv"), "Species\tGene_ID\nAethionema arabicum\tAa1G10\n")
w(os.path.join(batch_gff, "Athaliana.gff3"), PHYTOZOME)
w(os.path.join(batch_gff, "Aarabicum.gff3"), MAKER)
p, out = run(["--ids-dir", batch_ids, "--gff-dir", batch_gff])
check("batch mode: both species found via the guessed-prefix heuristic (Athaliana, Aarabicum), exit 0",
      p.returncode == 0, p.stderr)
st = {r["species"]: r for r in table(os.path.join(out, "stats.tsv"))}
check("batch mode: species names come from the ID file names (underscores as spaces)",
      "Arabidopsis thaliana" in st and "Aethionema arabicum" in st, st)
check("batch mode: BED filenames are filesystem-safe (no literal space)",
      set(os.listdir(out)) >= {"Arabidopsis_thaliana.bed", "Aethionema_arabicum.bed"}, os.listdir(out))

# --- batch mode with an explicit --species-map (exact pairing, not guessed)
smap = os.path.join(tmp, "smap.tsv")
w(smap, "Species\tBasename\nArabidopsis thaliana\tAthaliana\nAethionema arabicum\tAarabicum\n")
w(os.path.join(batch_gff, "Athaliana.gff3"), PHYTOZOME)  # re-affirm same content
p, out = run(["--ids-dir", batch_ids, "--gff-dir", batch_gff, "--species-map", smap])
check("--species-map gives the same pairing as the heuristic here (sanity)", p.returncode == 0, p.stderr)

# --- match_files collision handling (real code, not just mocked)
amb_ids = os.path.join(tmp, "amb_ids"); amb_gff = os.path.join(tmp, "amb_gff")
os.makedirs(amb_ids); os.makedirs(amb_gff)
w(os.path.join(amb_ids, "Foo_bar.tsv"), "Species\tGene_ID\nFoo bar\tx.1\n")
w(os.path.join(amb_gff, "Fbar_one.gff"), "")
w(os.path.join(amb_gff, "Fbar_two.gff"), "")
try:
    m.match_files(amb_ids, amb_gff)
    check("two GFFs matching the same guessed prefix is refused, not silently resolved", False)
except h.HepError as e:
    check("two GFFs matching the same guessed prefix is refused, not silently resolved", "more than one GFF" in str(e), str(e))

os.makedirs(os.path.join(tmp, "amb2_ids")); os.makedirs(os.path.join(tmp, "amb2_gff"))
w(os.path.join(tmp, "amb2_ids", "Foo_bar.tsv"), "Species\tGene_ID\nFoo bar\tx.1\n")
w(os.path.join(tmp, "amb2_ids", "Foo_baz.tsv"), "Species\tGene_ID\nFoo baz\ty.1\n")
w(os.path.join(tmp, "amb2_gff", "Fbar.gff"), "")
try:
    m.match_files(os.path.join(tmp, "amb2_ids"), os.path.join(tmp, "amb2_gff"))
except h.HepError:
    check("two different ID files racing for one GFF is a non-issue when only one legitimately matches", False)
else:
    check("two different ID files racing for one GFF is a non-issue when only one legitimately matches", True)

# a genuine race: two DIFFERENT id-file stems whose guessed candidates both single-match the SAME one gff
os.makedirs(os.path.join(tmp, "amb3_ids")); os.makedirs(os.path.join(tmp, "amb3_gff"))
w(os.path.join(tmp, "amb3_ids", "Foo_bar.tsv"), "Species\tGene_ID\nFoo bar\tx.1\n")     # candidates: Foo_bar, Foo, Fbar
w(os.path.join(tmp, "amb3_ids", "Foo_baz.tsv"), "Species\tGene_ID\nFoo baz\ty.1\n")     # candidates: Foo_baz, Foo, Fbaz
w(os.path.join(tmp, "amb3_gff", "Foo.gff"), "")                                          # starts with "Foo": matches BOTH
try:
    m.match_files(os.path.join(tmp, "amb3_ids"), os.path.join(tmp, "amb3_gff"))
    check("two id files that both single-match the same GFF ('Foo') is refused, not silently overwritten", False)
except h.HepError as e:
    check("two id files that both single-match the same GFF ('Foo') is refused, not silently overwritten",
          "matched the same GFF file" in str(e), str(e))

# --- refusals
p, out = run(["--gff", os.path.join(tmp, "maker.gff3")])
check("--ids without --gff is refused", p.returncode != 0 and "needs both" in p.stderr, p.stderr)
p, out = run(["--ids", os.path.join(tmp, "maker_ids.tsv")])
check("--gff without --ids is refused", p.returncode != 0 and "needs both" in p.stderr, p.stderr)
p, out = run(["--ids", os.path.join(tmp, "maker_ids.tsv"), "--gff", os.path.join(tmp, "maker.gff3"),
             "--ids-dir", batch_ids])
check("mixing single and batch flags is refused", p.returncode != 0 and "not a mix" in p.stderr, p.stderr)
p, out = run(["--ids", os.path.join(tmp, "nowhere.tsv"), "--gff", os.path.join(tmp, "maker.gff3")])
check("a missing ID file does not crash the run: reported per-species, not a Python traceback",
      "Traceback" not in p.stderr, p.stderr)
p1, out1 = run(["--ids", os.path.join(tmp, "maker_ids.tsv"), "--gff", os.path.join(tmp, "maker.gff3")])
p2 = subprocess.run([sys.executable, SCRIPT, "--ids", os.path.join(tmp, "maker_ids.tsv"), "--gff",
                    os.path.join(tmp, "maker.gff3"), "-o", out1, "--quiet"], stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, universal_newlines=True)
check("an existing output folder is not overwritten without --force", p2.returncode != 0 and "already holds results" in p2.stderr, p2.stderr)
p3 = subprocess.run([sys.executable, SCRIPT, "--ids", os.path.join(tmp, "maker_ids.tsv"), "--gff",
                    os.path.join(tmp, "maker.gff3"), "-o", out1, "--quiet", "--force"], stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, universal_newlines=True)
check("--force overwrites it", p3.returncode == 0, p3.stderr)

# --- --id-column and a comma-separated ID list
w(os.path.join(tmp, "csv_ids.csv"), "Protein,Note\nAT1G01010.1,foo\n")
p, out = run(["--ids", os.path.join(tmp, "csv_ids.csv"), "--gff", os.path.join(tmp, "phytozome.gff3"), "--id-column", "Protein"])
check("--id-column reads a comma-separated file by an explicit column name", p.returncode == 0
      and table(os.path.join(out, "master.tsv"))[0]["matched_id"] == "AT1G01010.1", p.stderr)
p, out = run(["--ids", os.path.join(tmp, "csv_ids.csv"), "--gff", os.path.join(tmp, "phytozome.gff3"), "--id-column", "Nope"])
check("an unknown --id-column is loud (warned) and counted as incomplete (exit 2), never silent",
      p.returncode == 2 and "no column" in p.stderr, p.stderr)
check("...and the species is recorded in stats.tsv as an error, not silently dropped from the run",
      table(os.path.join(out, "stats.tsv"))[0]["status"].startswith("error:"))

# --- --version
pv = subprocess.run([sys.executable, SCRIPT, "--version"], stdout=subprocess.PIPE, universal_newlines=True)
check("--version prints the pipeline version", pv.stdout.strip() == "%s %s" % (h.PIPELINE_NAME, h.__version__))

# --- pieces
check("normalize_id strips a trailing '.p' and nothing else", m.normalize_id("AT1G01010.1.p") == "AT1G01010.1"
      and m.normalize_id("AT1G01010.1") == "AT1G01010.1" and m.normalize_id("") == "")

print()
print("FAILED: %d" % len(fails) if fails else "ALL PASSED")
sys.exit(1 if fails else 0)
