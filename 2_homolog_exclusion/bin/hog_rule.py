#!/usr/bin/env python3
"""
hog_rule.py - homology through a hierarchical orthogroup (HOG).

Question answered: a candidate passed every search, but OrthoFinder placed it in the same HOG as a gene that a search removed.
Is the candidate homologous to the part of that gene that produced its hits?

For every gene that combine marked PASS and that shares a HOG with a gene removed by an evidence run, the candidate is aligned
with BLASTp to those removed members. Each member's evidence region is the union of the query intervals of its counted hits
(column q_intervals of the runs' hits.tsv, written since version 0.2.0). The candidate is removed (REMOVED_HOG) when one
alignment to a member overlaps that member's evidence region by at least min_overlap_aa amino acids. Genes not in such a HOG
keep their PASS. Nothing else is changed: genes combine did not pass are not touched.

  run:         python3 hog_rule.py -i proteins/ -o out/hog --combine out/combine --hogs OrthoFinder/.../N0.tsv \
                   --set evidence_runs=blast,jackhmmer,cdd,interpro
  parse only:  python3 hog_rule.py ... --parse-only earlier_out/work/hog_blastp.tsv
Python 3.11+, standard library only; needs BLAST+ (blastp, makeblastdb) unless --parse-only is given.
"""
import csv
import fnmatch
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402

MODULE = "hog-rule"
DESCRIPTION = __doc__.strip().split("\n\n")[0]
OUTFMT = "qseqid sseqid evalue bitscore qstart qend sstart send"
PASS, REMOVED_HOG, REGION_UNKNOWN, NOT_RUN = "PASS", "REMOVED_HOG", "REGION_UNKNOWN", "NOT_RUN"


def read_hogs(path):
    """{gene_id: hog} from an OrthoFinder Phylogenetic_Hierarchical_Orthogroups N0.tsv (HOG, OG, parent clade, species...)."""
    if not os.path.isfile(path):
        raise h.HepError("HOG table not found: %s" % path)
    out = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd, None)
        if not head or head[0] != "HOG":
            raise h.HepError("%s: expected an OrthoFinder N0.tsv with a 'HOG' first column" % path)
        first = 3 if len(head) > 2 and head[2].startswith("Gene Tree") else 2
        for row in rd:
            for cell in row[first:]:
                for g in cell.split(","):
                    g = h.canonical_id(g.strip())
                    if g:
                        if g in out and out[g] != row[0]:
                            raise h.HepError("%s: gene %s is in two HOGs (%s, %s)" % (path, g, out[g], row[0]))
                        out[g] = row[0]
    return out


def read_matrix(combine_dir):
    path = os.path.join(combine_dir, "gene_matrix.tsv")
    if not os.path.isfile(path):
        raise h.HepError("%s not found: give the output folder of combine" % path)
    rows = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[(r["species"], r["gene_id"])] = r
    return rows


def evidence_regions(run_paths):
    """{(species, gene): [intervals]} from the counted hits of the given run folders, plus the set of genes with a counted hit
    but no coordinates."""
    regions, unknown = {}, set()
    for p in run_paths:
        with open(os.path.join(p, "hits.tsv"), "r", encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("counts_as_hit") != "yes":
                    continue
                key = (r["species"], r["gene_id"])
                iv = h.parse_intervals(r.get("q_intervals"))
                if iv:
                    regions.setdefault(key, []).extend(iv)
                else:
                    unknown.add(key)
    return {k: h.merge_intervals(v) for k, v in regions.items()}, unknown - set(regions)


def decide(alignments, regions, unknown, min_overlap, unknown_region="fail"):
    """
    alignments: {candidate: {member: [(sstart, send), ...]}} (subject = member coordinates); regions: {member: intervals};
    unknown: members with evidence but no coordinates. Returns {candidate: (verdict, best_member, best_overlap_aa)}.
    The overlap with a member is the overlap of the union of the candidate's alignments on that member with its region.
    """
    out = {}
    for cand, by_member in alignments.items():
        best, best_m, saw_unknown = 0, None, False
        for m, iv in sorted(by_member.items()):
            if m in regions:
                ov = h.overlap_length(iv, regions[m])
                if ov > best:
                    best, best_m = ov, m
            elif m in unknown:
                saw_unknown = True
        if best >= min_overlap:
            out[cand] = (REMOVED_HOG, best_m, best)
        elif saw_unknown and unknown_region == "fail":
            out[cand] = (REGION_UNKNOWN, best_m, best)
        else:
            out[cand] = (PASS, best_m, best)
    return out


def parse_hog_blast(path, pairs_ok, evalue):
    """{query: {subject: [(sstart, send)]}} for alignments with E <= evalue between a candidate and a member of its HOG."""
    confirmed, unconfirmed, rows = h.read_blast_tabular(path, len(OUTFMT.split()), OUTFMT)
    if unconfirmed:
        raise h.HepError("%s: the within-HOG BLASTp did not finish for %d candidate(s), e.g. %s"
                         % (path, len(unconfirmed), sorted(unconfirmed)[0]))
    al = {}
    for c in rows:
        q, s = h.canonical_id(c[0]), h.canonical_id(c[1])
        if float(c[2]) <= evalue and (q, s) in pairs_ok:
            al.setdefault(q, {}).setdefault(s, []).append((int(c[6]), int(c[7])))
    return al, confirmed


def add_args(p):
    p.add_argument("--combine", required=True, help="output folder of combine (gene_matrix.tsv, run.json)")
    p.add_argument("--hogs", required=True, help="OrthoFinder Phylogenetic_Hierarchical_Orthogroups/N0.tsv")
    p.add_argument("--parse-only", default=None, help="read this earlier within-HOG BLASTp output instead of running BLAST+")


def main(argv=None):
    ctx = h.start_module(MODULE, DESCRIPTION, argv, extra=add_args)
    a, cfg = ctx["args"], ctx["cfg"]
    try:
        if cfg.get("unknown_region") not in ("fail", "pass"):
            raise h.HepError("unknown_region must be 'fail' or 'pass'")
        if int(cfg["min_overlap_aa"]) < 1:
            raise h.HepError("min_overlap_aa must be at least 1")
        matrix = read_matrix(a.combine)
        with open(os.path.join(a.combine, "run.json"), "r", encoding="utf-8") as fh:
            cinfo = json.load(fh)
        runs = {r["label"]: r["path"] for r in cinfo.get("runs", [])}
        pats = h.as_list(cfg.get("evidence_runs")) or ["*"]
        missing = [p for p in pats if not fnmatch.filter(runs, p)]
        want = [l for l in runs if any(fnmatch.fnmatchcase(l, p) for p in pats)]
        if missing:
            raise h.HepError("evidence_runs not found in combine's run.json: %s (known: %s)" % (", ".join(missing), ", ".join(runs)))
        if not cfg.get("evidence_runs"):
            h.warn("evidence_runs not set: every run counts as evidence (%s)" % ", ".join(runs))
        regions, unknown = evidence_regions([runs[l] for l in want])
        hogs = read_hogs(a.hogs)
        seqs = {}
        for sp, path in ctx["inputs"]:
            for rid, header, seq in h.read_fasta(path):
                seqs[(sp, h.canonical_id(rid))] = seq
        gene_sp = {}
        for (sp, g) in matrix:
            if g in gene_sp and gene_sp[g] != sp:
                raise h.HepError("gene ID %s occurs in two species (%s, %s); HOG membership would be ambiguous" % (g, gene_sp[g], sp))
            gene_sp[g] = sp
        evid = {g for (sp, g), r in matrix.items() if r["verdict"] == "REMOVED" and set(r["removed_by"].split(";")) & set(want)}
        hog_members = {}
        for g in evid:
            if g in hogs:
                hog_members.setdefault(hogs[g], set()).add(g)
        cands = {g for (sp, g), r in matrix.items() if r["verdict"] == "PASS" and hogs.get(g) in hog_members}
        pairs = {(c, m) for c in cands for m in hog_members[hogs[c]]}
        h.log("%d PASS genes share a HOG with %d removed evidence members" % (len(cands), len({m for _, m in pairs})))
        h.check_output_folder(a.output, MODULE, a.label or MODULE, resume=False, force=a.force)
    except h.HepError as e:
        h.die(str(e))
    os.makedirs(os.path.join(a.output, "work"), exist_ok=True)
    with open(os.path.join(a.output, h.RUN_MARKER), "w", encoding="utf-8") as fh:
        fh.write("module=%s\nlabel=%s\n" % (MODULE, a.label or MODULE))
    workdir = os.path.join(a.output, "work")
    tools = []
    if cands:
        if a.parse_only:
            bl = a.parse_only
        else:
            for t in ("blastp", "makeblastdb"):
                if not h.find_tool(t):
                    h.die("%s not found on PATH (or use --parse-only)" % t)
            mfa, cfa, bl = (os.path.join(workdir, n) for n in ("members.faa", "candidates.faa", "hog_blastp.tsv"))
            for path, genes in ((mfa, {m for _, m in pairs}), (cfa, cands)):
                with open(path, "w", encoding="utf-8") as fh:
                    for g in sorted(genes):
                        if (gene_sp[g], g) not in seqs:
                            h.die("no sequence for %s in the input FASTA files" % g)
                        h.write_fasta_record(fh, g, seqs[(gene_sp[g], g)])
            if h.run_command(["makeblastdb", "-in", mfa, "-dbtype", "prot", "-out", mfa + ".db"], os.path.join(workdir, "makeblastdb.log")):
                h.die("makeblastdb failed (see work/makeblastdb.log)")
            cmd = ["blastp", "-query", cfa, "-db", mfa + ".db", "-evalue", "%g" % float(cfg["evalue"]), "-seg",
                   "yes" if cfg.get("seg") else "no", "-dbsize", str(int(cfg["dbsize"])), "-max_target_seqs",
                   str(int(cfg["max_target_seqs"])), "-num_threads", str(a.threads), "-outfmt", "7 " + OUTFMT, "-out", bl]
            if h.run_command(cmd, os.path.join(workdir, "blastp.log")):
                h.die("blastp failed (see work/blastp.log)")
            tools.append({"name": "BLAST+", "version": h.tool_version(["blastp", "-version"])})
        try:
            al, searched = parse_hog_blast(bl, pairs, float(cfg["evalue"]))
        except h.HepError as e:
            h.die(str(e))
        dec = decide(al, {g: iv for (_, g), iv in regions.items()}, {g for _, g in unknown}, int(cfg["min_overlap_aa"]),
                     cfg["unknown_region"])
        done = {h.canonical_id(q) for q in searched}
        lost = sorted(cands - done)
        if lost:
            h.warn("%d candidate(s) are not in the within-HOG BLASTp output (e.g. %s): NOT_RUN, they do not pass"
                   % (len(lost), lost[0]))
            for g in lost:
                dec[g] = (NOT_RUN, None, 0)
    else:
        dec = {}
    regions_g = {g: iv for (_, g), iv in regions.items()}
    rows = [["species", "gene_id", "hog", "combine_verdict", "verdict", "best_member", "overlap_aa", "member_region"]]
    for (sp, g), r in sorted(matrix.items()):
        if r["verdict"] != "PASS":
            continue
        v, m, ov = dec.get(g, (PASS, None, 0))
        rows.append([sp, g, hogs.get(g, "NA"), r["verdict"], v, m or "NA", ov, h.format_intervals(regions_g.get(m, [])) or "NA"])
    with open(os.path.join(a.output, "hog_rule.tsv"), "w", encoding="utf-8", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    passed = [(x[0], x[1]) for x in rows[1:] if x[4] == PASS]
    with open(os.path.join(a.output, "passed_ids.tsv"), "w", encoding="utf-8") as fh:
        fh.write("species\tgene_id\n")
        for sp, g in passed:
            fh.write("%s\t%s\n" % (sp, g))
    pdir = os.path.join(a.output, "passed"); os.makedirs(pdir, exist_ok=True)
    keep = set(passed)
    for sp, path in ctx["inputs"]:
        with open(os.path.join(pdir, sp + ".faa"), "w", encoding="utf-8") as fh:
            for rid, header, seq in h.read_fasta(path):
                if (sp, h.canonical_id(rid)) in keep:
                    h.write_fasta_record(fh, header, seq)
    tally = {}
    for x in rows[1:]:
        tally.setdefault(x[0], {PASS: 0, REMOVED_HOG: 0, REGION_UNKNOWN: 0, NOT_RUN: 0})[x[4]] += 1
    with open(os.path.join(a.output, "summary.tsv"), "w", encoding="utf-8") as fh:
        fh.write("species\tcombine_pass\t%s\t%s\t%s\t%s\n" % (PASS, REMOVED_HOG, REGION_UNKNOWN, NOT_RUN))
        for sp in sorted(tally):
            t = tally[sp]
            fh.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (sp, sum(t.values()), t[PASS], t[REMOVED_HOG], t[REGION_UNKNOWN], t[NOT_RUN]))
    info = h.run_info_base(ctx)
    info.update(module=MODULE, pipeline=h.PIPELINE_NAME, version=h.__version__, tools=tools,
                combine=os.path.abspath(a.combine), hogs=os.path.abspath(a.hogs), evidence_runs=want,
                n_candidates=len(cands), n_removed_hog=sum(v[0] == REMOVED_HOG for v in dec.values()),
                n_region_unknown=sum(v[0] == REGION_UNKNOWN for v in dec.values()), n_pass=len(passed),
                n_not_run=sum(v[0] == NOT_RUN for v in dec.values()),
                members_without_region=len(unknown))
    with open(os.path.join(a.output, "run.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
    if not a.keep_work:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)
    h.log("hog-rule: %d candidates, %d removed (overlap >= %s aa), %d region unknown, %d pass"
          % (len(cands), info["n_removed_hog"], cfg["min_overlap_aa"], info["n_region_unknown"], len(passed)))
    if info["n_region_unknown"] or info["n_not_run"]:
        h.note_incomplete(info["n_region_unknown"] + info["n_not_run"])


if __name__ == "__main__":
    main()
    sys.exit(h.exit_status())
