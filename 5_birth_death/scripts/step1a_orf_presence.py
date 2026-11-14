#!/usr/bin/env python3
"""
step1a_orf_presence.py - presence of each candidate's ORF in every extant genome by the stage 4 ORF rule, in place of the
stage 4 tBLASTn presence (from the T2T mammal analysis, LINGUA_t2t stage5_t2t/build_length_presence_t2t.py, 2026-11-13).

For every row of each stage 4 tBLASTn summary (<prefix>_cdngs_tblastn_<n>_summary.tsv) and every species column, from the
stage 3 alignment of the locus (<Focal>_<gene>_*.aln, in the alignment folder whose files start with "<Focal>_"):
  present (1)  the genome's gap-free row holds an ORF that passes the stage 4 step 4b test with the thresholds of the stage 4
               method options (orf_min_len, orf_min_cov, orf_min_in_frame, orf_count_stop), computed with the functions of
               4_enablers/scripts/orf_length_node.py
  unknown      not present, and the row is missing or has bases (not '-' or N) at fewer than --unknown-min-aligned % of
               the focal CDS columns
  absent (0)   otherwise
The focal species is always present. A locus with no stage 3 alignment keeps the focal species only and lists every other
genome as unknown. Phytozome IDs ending in ".p" are looked up without the suffix.

Outputs in OUTDIR:
  <prefix>_cdngs_tblastn_<n>_summary.tsv  same name and columns as the tBLASTn summary, 1/0 (unknown written as 0)
  unknown_presence.tsv                    species, geneID, unknown_in (';'-joined genomes), for step 5 --unknown
  presence_comparison.tsv                 per focal species, genome cells by tBLASTn call and new call

Usage: python3 step1a_orf_presence.py --tree-config TREE_CONFIG --tblastn TBLASTN_DIR --alignments ALIGNMENTS_DIR
           --stage4-options METHOD_OPTIONS_JSON --outdir OUTDIR [--stage4-scripts DIR] [--unknown-min-aligned 50]
"""
import argparse
import collections
import csv
import glob
import json
import os
import sys

csv.field_size_limit(10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__))


def die(msg):
    sys.exit(f"ERROR: {msg}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--tree-config", required=True)
    ap.add_argument("--tblastn", required=True, help="folder of stage 4 *_cdngs_tblastn_*_summary.tsv")
    ap.add_argument("--alignments", required=True, help="folder holding one stage 3 alignment folder per focal species")
    ap.add_argument("--stage4-options", required=True, help="the stage 4 method options file of the run")
    ap.add_argument("--stage4-scripts", default=os.path.join(HERE, "..", "..", "4_enablers", "scripts"))
    ap.add_argument("--unknown-min-aligned", type=float, default=50.0)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    for p in (a.tree_config, a.tblastn, a.alignments, a.stage4_options, a.stage4_scripts):
        if not os.path.exists(p):
            die(f"not found: {p}")
    if os.path.exists(a.outdir):
        die(f"{a.outdir} exists; reruns go to a new folder")
    opts = json.load(open(a.stage4_options))
    if opts.get("orf_test") != "length":
        die(f"{a.stage4_options} does not use the length ORF test (orf_test = {opts.get('orf_test')!r}); presence by the "
            "ORF rule would not match stage 4, so keep the tBLASTn presence")
    min_len, min_cov = float(opts.get("orf_min_len", 80)), float(opts.get("orf_min_cov", 80))
    min_inf, count_stop = float(opts.get("orf_min_in_frame", 0) or 0), bool(opts.get("orf_count_stop", False))
    os.environ["LINGUA_CONFIG"] = a.tree_config
    os.environ["LINGUA_OPTIONS"] = a.stage4_options
    sys.path.insert(0, a.stage4_scripts)
    from orf_length_node import longest_orf, in_frame_pct, read_aln  # noqa: E402
    al = json.load(open(a.tree_config))["aliases"]
    print(f"ORF rule: >= {min_len:g}% length, >= {min_cov:g}% coverage, >= {min_inf:g}% in frame, stop counted: {count_stop}; "
          f"unknown below {a.unknown_min_aligned:g}% aligned", file=sys.stderr)

    folders = [d for d in sorted(glob.glob(os.path.join(a.alignments, "*"))) if os.path.isdir(d)]
    summaries = sorted(glob.glob(os.path.join(a.tblastn, "*_cdngs_tblastn_*_summary.tsv")))
    if not summaries:
        die(f"no *_cdngs_tblastn_*_summary.tsv in {a.tblastn}")
    os.makedirs(a.outdir)
    unk_rows, comp = [], collections.Counter()
    for summ in summaries:
        pre = os.path.basename(summ).split("_cdngs_tblastn")[0]
        if pre not in al:
            die(f"{pre}: not a species or alias of {a.tree_config}")
        F = al[pre]
        alnd = [d for d in folders if any(os.path.basename(f).startswith(F + "_") for f in glob.glob(os.path.join(d, "*.aln"))[:50])]
        if len(alnd) != 1:
            die(f"{F}: {len(alnd)} alignment folders in {a.alignments}")
        rows = list(csv.DictReader(open(summ), delimiter="\t"))
        if not rows:
            die(f"empty: {summ}")
        fields = list(rows[0].keys())
        spcols = [c for c in fields if c not in ("geneID", "Node_tblastn")]
        new, no_aln = [], 0
        for r in rows:
            g = r["geneID"]
            gk = g[:-2] if g.endswith(".p") else g
            p = glob.glob(os.path.join(alnd[0], f"{F}_{gk}_*.aln"))
            if len(p) > 1:
                die(f"{F} {g}: {len(p)} alignments")
            aln = read_aln(p[0]) if p else None
            no_aln += aln is None
            if aln is not None:
                ref = aln[F]
                hc = [i for i, ch in enumerate(ref) if ch != "-"]
                n = len(hc)
                starts = {"ATG", "".join(ref[hc[k]] for k in range(3))}
            o = {"geneID": g, "Node_tblastn": r["Node_tblastn"]}
            unknown = []
            for sp in spcols:
                spc = al.get(sp, sp)
                if spc == F:
                    call = "1"
                elif aln is None or spc not in aln:
                    call = "unknown"
                else:
                    row = aln[spc]
                    L, c0, c1 = longest_orf(row, starts, count_stop)
                    cov = 100.0 * sum(1 for i in hc if c0 is not None and c0 <= i <= c1) / n
                    if 100.0 * L / n >= min_len and cov >= min_cov and in_frame_pct(ref, row) >= min_inf:
                        call = "1"
                    elif 100.0 * sum(1 for i in hc if row[i] not in "-N") / n < a.unknown_min_aligned:
                        call = "unknown"
                    else:
                        call = "0"
                o[sp] = "1" if call == "1" else "0"
                if call == "unknown":
                    unknown.append(spc)
                if spc != F:
                    comp[(F, "1" if r[sp].strip() == "1" else "0", call)] += 1
            new.append(o)
            if unknown:
                unk_rows.append([F, g, ";".join(sorted(unknown))])
        with open(os.path.join(a.outdir, os.path.basename(summ)), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(new)
        print(f"{F}: {len(new)} candidates" + (f", {no_aln} without an alignment" if no_aln else ""), file=sys.stderr)
    with open(os.path.join(a.outdir, "unknown_presence.tsv"), "w") as fh:
        fh.write("species\tgeneID\tunknown_in\n")
        for x in unk_rows:
            fh.write("\t".join(x) + "\n")
    with open(os.path.join(a.outdir, "presence_comparison.tsv"), "w") as fh:
        fh.write("focal_species\ttblastn\tnew\tcells\n")
        for (F, t, nw), c in sorted(comp.items()):
            fh.write(f"{F}\t{t}\t{nw}\t{c}\n")


if __name__ == "__main__":
    main()
