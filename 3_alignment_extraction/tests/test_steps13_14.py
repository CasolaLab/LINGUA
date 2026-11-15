"""
test_steps13_14.py

Toy checks of stage 3 steps 13 and 14 as changed on 2026-11-15, each answer worked out by hand. One 33-nt reference CDS
(one exon, chr1:100-132, + strand) and four other rows:
  N4  lacks the first two bases:   step 13 reads the reference frame, so its stops are TAA at 7, TAG at 16 and the
                                   shared stop TAA at 31 (reading its own frame from base 3 would give only TGA at 12);
                                   step 14 does not report its 2-nt gap at 1-2 (an unaligned end)
  N6  1-nt deletion at 10:         internal, outside the mask, a frameshift -> the only reported indel
  N9  5-nt deletion at 20-24:      overlaps the missing-data mask (CDS 22-30) -> not reported, and not trimmed to 20-21
  N3  4-nt deletion at 30-33:      reaches the last base (an unaligned end) -> not reported
Usage: python3 3_alignment_extraction/tests/test_steps13_14.py   (or pytest); exits non-zero on failure
"""
import csv
import os
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(os.path.dirname(HERE), "scripts")
R = "Arabidopsis_thaliana"
REF = "ATGGCTCAAAGTGGAAAGCTTCATGGTGCATAA"
ROWS = {
    "N4": "--GGCTTAAAGTGAATAGCTTCATGGTGCATAA",
    "N6": REF[:9] + "-" + REF[10:],
    "N9": REF[:19] + "-----" + REF[24:],
    "N3": REF[:29] + "----",
}
FN = f"{R}_g1_chr1_100_132"


def w(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def setup(d):
    w(f"{d}/step9/{FN}.aln", f">{R}.chr1\n{REF}\n" + "".join(f">{k}.x\n{v}\n" for k, v in ROWS.items()))
    w(f"{d}/ref.bed", "chr1\t99\t132\tg1_exon1\t0\t+\tg1\n")
    w(f"{d}/summary.tsv", f"file_name\tgff_transcriptID\tchromosome\tstrand\tactual_nt\n{FN}\tg1\tchr1\t+\t33\n")
    w(f"{d}/mask.tsv", f"alignment\tspecies\texon\tref_start\tref_end\treason\n{FN}\tN9\t1\t120\t129\tabsent\n")


def rows(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def test_step13_reference_frame():
    assert all(len(v) == len(REF) for v in ROWS.values())
    with tempfile.TemporaryDirectory() as d:
        setup(d)
        subprocess.run(["bash", f"{S}/step13_identify_stop_codons.sh", f"{d}/step9", f"{d}/step13", f"{d}/ref.bed", R,
                        f"{d}/summary.tsv"], check=True, capture_output=True)
        got = sorted((r["Species"], int(r["StopCodonStart_RefPos"]), int(r["StopCodonEnd_RefPos"]), r["StopCodon"])
                     for r in rows(f"{d}/step13/master_stop.tsv"))
        exp = sorted([("N4", 7, 9, "TAA"), ("N4", 16, 18, "TAG"), ("N4", 31, 33, "TAA"),
                      ("N6", 31, 33, "TAA"), ("N9", 31, 33, "TAA")])  # N3 lacks the last codon
        assert got == exp, got
        r0 = rows(f"{d}/step13/master_stop.tsv")[0]
        assert (r0["gff_transcriptID"], r0["Chromosome"], r0["gff_strand_ref"], r0["reference_length"]) == ("g1", "chr1", "+", "33")
        assert os.path.exists(f"{d}/step13/{FN}_stops.tsv")


def test_step14_ends_and_mask():
    with tempfile.TemporaryDirectory() as d:
        setup(d)
        p = subprocess.run(["bash", f"{S}/step14_identify_indels.sh", f"{d}/step9", f"{d}/step14", f"{d}/ref.bed", R,
                            f"{d}/summary.tsv", f"{d}/mask.tsv", f"{d}/not_reported.tsv"],
                           check=True, capture_output=True, text=True)
        got = [(r["Species"], r["IndelStart_RefPos"], r["IndelEnd_RefPos"], r["IndelLength"]) for r in rows(f"{d}/step14/master_indels.tsv")]
        assert got == [("N6", "10", "10", "1")], got
        aud = {(r["Species"], r["action"]) for r in rows(f"{d}/not_reported.tsv")}
        assert aud == {("N4", "unaligned_end"), ("N9", "masked_out"), ("N3", "unaligned_end")}, aud
        assert "unaligned ends) not reported: 2" in p.stdout and "1 deletion runs overlapping them not reported" in p.stdout, p.stdout


if __name__ == "__main__":
    for t in (test_step13_reference_frame, test_step14_ends_and_mask):
        t()
        print(f"{t.__name__} PASSED")
