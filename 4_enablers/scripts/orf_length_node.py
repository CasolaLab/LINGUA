#!/usr/bin/env python3
"""
orf_length_node.py  (stage 4, step 4b; method option orf_test = "length")

Replace the tBLASTn-based ancestral ORF presence of step 4 with a length-based test (his principle:
ORF length matters far more than protein conservation). For every locus and every ancestor
(N label) in the stage 3 alignment, the longest ORF in the three forward frames of the ancestor's
gap-free row is found, starting at ATG or at the focal gene's own first codon (so annotated non-ATG
starts such as CUG are honoured) and ending at the first stop. The ancestor has the ORF when that
ORF is >= MIN_LEN % of the focal CDS length and spans >= MIN_COV % of the focal CDS positions. No
identity cutoff, any frame.

Rewrites ANCSEQ_TSV (same columns: geneID, Node_AncSeq_ORF = oldest focal-lineage ancestor with the
ORF or <focal>-specific, then one 0/1 column per ancestor), keeping the step 4 table as
<name>.tblastn.tsv (read from there on reruns), and writes <name>.orf_length_details.tsv with, per
locus and ancestor, the ORF length %, coverage %, in-frame % of the row (longest unbroken focal-frame stretch) and
in-frame % of the ORF (share of the focal CDS that the ORF reads in the focal frame). With orf_min_in_frame the
call also needs the ORF itself to be mostly in the focal frame (orf_frame_measure "orf", the default since
2026-11-14), or, with orf_frame_measure "row", the row stretch to reach the cutoff (the 2026-11-13 runs).
orf_frame_measure "contiguous" (T2T, 2026-11-15, his decision: disablers are the main signal of de novo birth, focal-frame
conservation decides whether an ancestor still has the ORF): the ancestor has the ORF when, from a start codon (ATG or the
focal first codon) read in the focal frame (the row bases aligned to one focal codon), the focal frame runs unbroken by a
stop codon or a frameshifting gap over >= orf_min_in_frame % of the focal CDS positions, and that stretch holds row bases
(not '-' or N) at >= MIN_LEN % of the focal CDS positions. Pieces of the focal frame separated by a frameshift are never
added together. Breaks are those of in_frame_pct (the focal stop codon itself is not a break).

Usage: python3 orf_length_node.py ANCSEQ_TSV ALN_DIR --focal-species Homo_sapiens
                                  [--min-len 80] [--min-cov 80]
"""
import argparse
import csv
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lingua_config as cfg  # noqa: E402

STOP = {"TAA", "TAG", "TGA"}


def die(m):
    sys.exit(f"ERROR: {m}")


def read_aln(path):
    rows, k = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split(".")[0]
            rows.setdefault(k, "")
        elif k:
            rows[k] += line.upper()
    return rows


def longest_orf(row, starts, count_stop=False):
    """Longest ORF of the 3 forward frames of the gap-free row: (length in nt, first column, last column).
    count_stop (method option orf_count_stop): an ORF ended by a stop codon includes it, as the focal CDS length
    does, so an ORF identical to the focal one scores 100% (without it, (L - 3) / L)."""
    cols = [i for i, c in enumerate(row) if c != "-"]
    s = "".join(row[i] for i in cols)
    best = (0, None, None)
    for f in range(3):
        i = f
        while i + 3 <= len(s):
            if s[i:i + 3] in starts:
                j = i
                while j + 3 <= len(s) and s[j:j + 3] not in STOP:
                    j += 3
                k = j + 3 if count_stop and j + 3 <= len(s) else j
                if k - i > best[0]:
                    best = (k - i, cols[i], cols[k - 1])
                i = j + 3
            else:
                i += 3
    return best


def all_orfs(row, starts, count_stop=False):
    """Every ORF of the 3 forward frames of the gap-free row, as (length in nt, first column, last column, columns of
    the first base of each codon); an ORF starts at the first start codon after the previous stop of its frame."""
    cols = [i for i, c in enumerate(row) if c != "-"]
    s = "".join(row[i] for i in cols)
    out = []
    for f in range(3):
        i = f
        while i + 3 <= len(s):
            if s[i:i + 3] in starts:
                j = i
                while j + 3 <= len(s) and s[j:j + 3] not in STOP:
                    j += 3
                k = j + 3 if count_stop and j + 3 <= len(s) else j
                if k > i:
                    out.append((k - i, cols[i], cols[k - 1], [cols[x] for x in range(i, k - 2, 3)]))
                i = j + 3
            else:
                i += 3
    return out


def orf_call(ref, row, min_len, min_cov, min_in_frame=0.0, count_stop=False, measure="orf"):
    """Step 4b ORF call of one ancestral row. Returns (has_orf, length %, coverage %, in-frame % of the row,
    in-frame % of the ORF), the last four for the ORF that decided the call (the longest ORF when none passes).
    measure "orf" (default; fix 2026-11-14): an ORF counts when it passes the length and coverage tests AND reads at
    least min_in_frame % of the focal CDS in the focal frame (its codons start on focal codon boundaries); every
    ORF is tried, so a long ORF in another frame no longer hides a shorter one in the focal frame. measure "row" (the
    2026-11-13 runs): the longest ORF passes length and coverage and, separately, the longest unbroken focal-frame
    stretch of the row (in_frame_pct, not tied to that ORF) is >= min_in_frame %. measure "contiguous" (T2T 2026-11-15):
    see focal_frame_orf; min_cov is not used (the stretch span is the coverage)."""
    hcols = [i for i, c in enumerate(ref) if c != "-"]
    n = len(hcols)
    pos = {c: p for p, c in enumerate(hcols)}
    inf = in_frame_pct(ref, row)

    def score(o):
        L, c0, c1, codons = o
        cov = 100.0 * sum(1 for i in hcols if c0 <= i <= c1) / n
        fr = 100.0 * 3 * sum(1 for c in codons if c in pos and pos[c] % 3 == 0) / n
        return 100.0 * L / n, cov, fr

    starts = {"ATG", "".join(ref[hcols[k]] for k in range(3))}
    orfs = sorted(all_orfs(row, starts, count_stop), key=lambda o: -o[0])
    if measure == "contiguous":
        # (has, length % = row bases in the stretch, coverage % = stretch span, row in-frame %, stretch span %)
        has, span, bases = focal_frame_orf(ref, row, starts, min_in_frame, min_len)
        return has, bases, span, inf, span
    if not orfs:
        return False, 0.0, 0.0, inf, 0.0
    if measure == "row":
        lp, cov, fr = score(orfs[0])
        return lp >= min_len and cov >= min_cov and inf >= min_in_frame, lp, cov, inf, fr
    if measure != "orf":
        die(f"orf_frame_measure must be orf, row or contiguous, not {measure}")
    for o in orfs:
        lp, cov, fr = score(o)
        if lp >= min_len and cov >= min_cov and fr >= min_in_frame:
            return True, lp, cov, inf, fr
    lp, cov, fr = score(orfs[0])
    return False, lp, cov, inf, fr


def in_frame_pct(ref, row):
    """Longest stretch of the focal CDS with no stop in focal codon columns and no frameshifting
    gap run, as % of the CDS (same measure as frame_check_orf_nodes.py)."""
    cols = [i for i, c in enumerate(ref) if c != "-"]
    n = len(cols)
    cuts = set()
    for k in range(0, n - 3 - (n % 3), 3):
        if "".join(row[cols[k + j]] for j in range(3)) in STOP:
            cuts.add(k + 1)
    i, refpos = 0, 0
    while i < len(ref):
        if ref[i] != "-":
            refpos += 1
        if ref[i] != "-" and row[i] == "-":
            j = i
            while j < len(ref) and ref[j] != "-" and row[j] == "-":
                j += 1
            if (j - i) % 3:
                cuts.add(refpos)
            refpos += j - i - 1
            i = j
            continue
        if ref[i] == "-" and row[i] != "-":
            j = i
            while j < len(ref) and ref[j] == "-" and row[j] != "-":
                j += 1
            if (j - i) % 3:
                cuts.add(refpos)
            i = j
            continue
        i += 1
    best, prev = 0, 0
    for c in sorted(cuts):
        best = max(best, c - 1 - prev)
        prev = c
    return 100.0 * max(best, n - prev) / n


def frame_breaks(ref, row):
    """Focal CDS positions where the focal frame of `row` breaks (contiguous measure): the first base of a stop codon read in
    the focal frame (the focal terminal codon excepted), the first focal position of a row gap run whose length is not a
    multiple of 3, and the focal position before a row insertion (bases in focal-gap columns) not a multiple of 3. The same
    rules as in_frame_pct, written separately so that function stays exactly as it was."""
    cols = [i for i, c in enumerate(ref) if c != "-"]
    n = len(cols)
    cuts = set()
    for k in range(0, n - 3 - (n % 3), 3):
        if "".join(row[cols[k + j]] for j in range(3)) in STOP:
            cuts.add(k + 1)
    refpos, i = 0, 0
    while i < len(ref):
        if ref[i] != "-":
            refpos += 1
            if row[i] == "-":
                j = i
                while j < len(ref) and ref[j] != "-" and row[j] == "-":
                    j += 1
                if (j - i) % 3:
                    cuts.add(refpos)
                refpos += j - i - 1
                i = j
                continue
        elif row[i] != "-":
            j = i
            while j < len(ref) and ref[j] == "-" and row[j] != "-":
                j += 1
            if (j - i) % 3:
                cuts.add(max(refpos, 1))
            i = j
            continue
        i += 1
    return sorted(cuts)


def focal_frame_orf(ref, row, starts, min_span, min_bases):
    """Contiguous measure. Returns (has_orf, span %, bases %) of the best focal-frame ORF of `row`: for every focal codon whose
    three aligned row bases form a start codon, the stretch runs from its first base to the base before the next break
    (frame_breaks), or to the CDS end; span = stretch length, bases = row bases (not '-' or N) inside it, both as % of the
    focal CDS. has_orf when some stretch reaches min_span and min_bases; else the longest stretch is reported."""
    cols = [i for i, c in enumerate(ref) if c != "-"]
    n = len(cols)
    cuts = frame_breaks(ref, row)
    based = [0]
    for c in cols:
        based.append(based[-1] + (row[c] not in "-N"))
    best = (False, 0.0, 0.0)
    for k in range(0, n - 2, 3):
        if "".join(row[cols[k + j]] for j in range(3)) not in starts:
            continue
        s = k + 1
        nxt = [c for c in cuts if c > s]
        end = nxt[0] - 1 if nxt else n
        span = 100.0 * (end - s + 1) / n
        bases = 100.0 * (based[end] - based[s - 1]) / n
        ok = span >= min_span and bases >= min_bases
        if ok and not best[0]:
            best = (True, span, bases)
        elif ok == best[0] and span > best[1]:
            best = (ok, span, bases)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("ancseq")
    ap.add_argument("aln_dir")
    ap.add_argument("--focal-species", required=True)
    ap.add_argument("--min-len", type=float, default=80.0)
    ap.add_argument("--min-cov", type=float, default=80.0)
    a = ap.parse_args()
    F = a.focal_species
    backup = a.ancseq[:-4] + ".tblastn.tsv"
    details = a.ancseq[:-4] + ".orf_length_details.tsv"
    if not os.path.exists(backup):
        if not os.path.exists(a.ancseq):
            die(f"not found: {a.ancseq}")
        shutil.copy2(a.ancseq, backup)
    src = list(csv.DictReader(open(backup), delimiter="\t"))
    if not src:
        die(f"empty: {backup}")
    fields = list(src[0].keys())
    nodes = [c for c in fields if c not in ("geneID", "Node_AncSeq_ORF")]
    lineage = sorted([n for n, c in cfg.NODE_CLADES.items() if F in c], key=lambda n: -len(cfg.NODE_CLADES[n]))
    spec = f"{F}-specific"
    min_inf = float(cfg.opt("orf_min_in_frame", 0) or 0)
    measure = cfg.OPTIONS.get("orf_frame_measure", "orf")  # absent = "orf", the fixed measure

    out, det = [], []
    changed, no_aln = 0, []
    for r in src:
        g = r["geneID"]
        gk = g[:-2] if g.endswith(".p") else g  # Phytozome IDs carry ".p" in the protein tables, not in alignments
        p = glob.glob(os.path.join(a.aln_dir, f"{F}_{gk}_*.aln"))
        if len(p) > 1:
            die(f"{len(p)} alignments for {g} in {a.aln_dir}")
        if not p:
            # stage 3 wrote no alignment for this locus, so step 4 searched nothing either: its row is kept as it is
            # (Brassicaceae 2026-11-13: 8 loci, all with an ORF in no ancestor); one with an ORF somewhere is an error
            if any(r[c] not in ("0", "0.0") for c in nodes):
                die(f"no alignment for {g} in {a.aln_dir}, but its step 4 row has an ORF in an ancestor")
            no_aln.append(g)
            out.append(dict(r))
            det.append({"geneID": g, "node": "NA", "row": "no_alignment", "orf_len_pct": "NA",
                        "orf_cov_pct": "NA", "in_frame_pct": "NA", "orf_in_frame_pct": "NA", "has_orf": "0"})
            continue
        aln = read_aln(p[0])
        ref = aln[F]
        hcols = [i for i, c in enumerate(ref) if c != "-"]
        n = len(hcols)
        starts = {"ATG", "".join(ref[hcols[k]] for k in range(3))}
        new = {"geneID": g}
        for node in nodes:
            row = aln.get(node)
            if row is None:
                new[node] = "0"
                det.append({"geneID": g, "node": node, "row": "absent", "orf_len_pct": "NA",
                            "orf_cov_pct": "NA", "in_frame_pct": "NA", "orf_in_frame_pct": "NA", "has_orf": "0"})
                continue
            # orf_min_in_frame (T2T 2026-11-13): the ancestral ORF counts only when it encodes mostly the same protein;
            # one definition of "ORF present" for dating, step 6c and step 7c. orf_frame_measure (2026-11-14) says how
            # the in-frame share is measured (see orf_call)
            has, lp, cov, inf, fr = orf_call(ref, row, a.min_len, a.min_cov, min_inf, cfg.opt("orf_count_stop"), measure)
            new[node] = "1" if has else "0"
            det.append({"geneID": g, "node": node, "row": "present", "orf_len_pct": f"{lp:.1f}",
                        "orf_cov_pct": f"{cov:.1f}", "in_frame_pct": f"{inf:.1f}", "orf_in_frame_pct": f"{fr:.1f}",
                        "has_orf": new[node]})
        new["Node_AncSeq_ORF"] = cfg.orf_node_from_cells(lambda x: new.get(x) == "1", lineage, spec)
        changed += new["Node_AncSeq_ORF"] != r["Node_AncSeq_ORF"]
        out.append(new)
    with open(a.ancseq, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out)
    with open(details, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(det[0]), delimiter="\t")
        w.writeheader()
        w.writerows(det)
    print(f"{os.path.basename(a.ancseq)}: length-based ORF test (>= {a.min_len:g}% length, >= {a.min_cov:g}% "
          f"coverage, any frame); ORF node changed for {changed} of {len(out)} loci"
          + (f"; {len(no_aln)} loci without an alignment kept as in step 4 ({', '.join(no_aln[:5])})" if no_aln else ""),
          file=sys.stderr)


if __name__ == "__main__":
    main()
