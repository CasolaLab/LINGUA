#!/usr/bin/env python3
"""
paralog_rule.py - a paralog among the established genes of the same genome.

Question answered: a candidate passed every search (and the HOG rule), but does its own genome hold a similar copy among the
genes that are not candidates? Such a copy means the candidate is a diverged duplicate of an established gene, not a gene
without homologs. OrthoFinder can miss these pairs (short or low-complexity proteins left unclustered), so they are searched
for directly.

Each passing gene is aligned with BLASTp (SEG off) to the full proteome of its species. Hits to itself and to other candidates
(any gene in the -i FASTA files) are ignored, because a copy that is itself a candidate makes a lineage-specific family, not a
duplicate of an established gene. The gene is removed (REMOVED_PARALOG) when one alignment to an established gene has at least
min_pident % identity over at least min_qcov % of the gene's length at E <= evalue.

  run:         python3 paralog_rule.py -i candidates/ -o out/paralog --passed out/hog --proteomes proteomes/
  parse only:  python3 paralog_rule.py ... --parse-only earlier_out/work/
Python 3.11+, standard library only; needs BLAST+ (blastp, makeblastdb) unless --parse-only is given.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402

MODULE = "paralog-rule"
DESCRIPTION = __doc__.strip().split("\n\n")[0]
OUTFMT = "qseqid sseqid pident evalue bitscore qstart qend qlen"
PASS, REMOVED_PARALOG, NOT_RUN = "PASS", "REMOVED_PARALOG", "NOT_RUN"


def read_passed(folder):
    path = os.path.join(folder, "passed_ids.tsv")
    if not os.path.isfile(path):
        raise h.HepError("%s not found: give the output folder of hog-rule (or of combine)" % path)
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return {(r["species"], h.canonical_id(r["gene_id"])) for r in csv.DictReader(fh, delimiter="\t")}


def best_paralogs(rows, candidates, evalue):
    """
    rows: BLAST lines (OUTFMT columns). candidates: IDs of every candidate of the species (the query itself is one of them).
    Returns {query: [(subject, pident, qcov_pct, evalue), ...]}, every alignment to a non-candidate gene with E <= evalue.
    qcov_pct is one alignment's span on the query over the query length.
    """
    out = {}
    for c in rows:
        q, s = h.canonical_id(c[0]), h.canonical_id(c[1])
        if s == q or s in candidates:
            continue
        ev = float(c[3])
        if ev > evalue:
            continue
        qcov = 100.0 * (int(c[6]) - int(c[5]) + 1) / int(c[7])
        out.setdefault(q, []).append((s, float(c[2]), qcov, ev))
    return out


def decide(alignments, min_pident, min_qcov):
    """{query: [(subject, pident, qcov_pct, evalue)]} -> {query: (verdict, subject, pident, qcov_pct, evalue)}.
    One alignment must pass both cuts; the reported paralog is the qualifying one with the highest identity, or, when none
    qualifies, the alignment with the highest identity."""
    out = {}
    for q, al in alignments.items():
        ok = [a for a in al if a[1] >= min_pident and a[2] >= min_qcov]
        pick = max(ok or al, key=lambda a: (a[1], a[2], -a[3]))
        out[q] = ((REMOVED_PARALOG if ok else PASS),) + pick
    return out


def add_args(p):
    p.add_argument("--passed", required=True, help="output folder of the previous step (hog-rule or combine), with passed_ids.tsv")
    p.add_argument("--proteomes", required=True,
                   help="folder of full proteomes, one FASTA per species, named like the -i files (longest isoform per gene)")
    p.add_argument("--parse-only", default=None,
                   help="folder of earlier BLASTp outputs (<species>.tsv, outfmt 7 %s) instead of running BLAST+" % OUTFMT)


def main(argv=None):
    ctx = h.start_module(MODULE, DESCRIPTION, argv, extra=add_args)
    a, cfg = ctx["args"], ctx["cfg"]
    try:
        min_pident, min_qcov, evalue = float(cfg["min_pident"]), float(cfg["min_qcov"]), float(cfg["evalue"])
        if not (0 < min_pident <= 100 and 0 < min_qcov <= 100):
            raise h.HepError("min_pident and min_qcov must be percentages above 0")
        passed = read_passed(a.passed)
        cand = {}
        seqs = {}
        for sp, path in ctx["inputs"]:
            for rid, header, seq in h.read_fasta(path):
                g = h.canonical_id(rid)
                cand.setdefault(sp, set()).add(g)
                seqs[(sp, g)] = seq
        lost = sorted(k for k in passed if k not in seqs)
        if lost:
            raise h.HepError("%d passed gene(s) are not in the -i FASTA files, e.g. %s %s" % (len(lost), lost[0][0], lost[0][1]))
        prot = dict(h.discover_inputs(a.proteomes, h.as_list(cfg.get("species_strip_suffixes")), ctx["species_map"]))
        species = sorted({sp for sp, _ in passed})
        nop = [sp for sp in species if sp not in prot]
        if nop:
            raise h.HepError("no proteome in %s for species: %s" % (a.proteomes, ", ".join(nop)))
        for sp in species:
            ids = {h.canonical_id(rid) for rid, _, _ in h.read_fasta(prot[sp])}
            absent = [g for (s, g) in passed if s == sp and g not in ids]
            if absent:
                h.warn("%s: %d passed gene(s) are not in the proteome (e.g. %s); their self hit cannot occur, the search still runs"
                       % (sp, len(absent), absent[0]))
        h.log("%d passing genes in %d species" % (len(passed), len(species)))
        h.check_output_folder(a.output, MODULE, a.label or MODULE, resume=False, force=a.force)
    except h.HepError as e:
        h.die(str(e))
    workdir = os.path.join(a.output, "work")
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(a.output, h.RUN_MARKER), "w", encoding="utf-8") as fh:
        fh.write("module=%s\nlabel=%s\n" % (MODULE, a.label or MODULE))
    tools = []
    if not a.parse_only:
        for t in ("blastp", "makeblastdb"):
            if not h.find_tool(t):
                h.die("%s not found on PATH (or use --parse-only)" % t)
        tools.append({"name": "BLAST+", "version": h.tool_version(["blastp", "-version"])})
    dec = {}
    for sp in species:
        genes = sorted(g for (s, g) in passed if s == sp)
        if a.parse_only:
            bl = os.path.join(a.parse_only, sp + ".tsv")
            if not os.path.isfile(bl):
                h.die("--parse-only: %s not found" % bl)
        else:
            qfa, bl, db = (os.path.join(workdir, sp + x) for x in (".query.faa", ".tsv", ".db"))
            with open(qfa, "w", encoding="utf-8") as fh:
                for g in genes:
                    h.write_fasta_record(fh, g, seqs[(sp, g)])
            if h.run_command(["makeblastdb", "-in", prot[sp], "-dbtype", "prot", "-out", db],
                             os.path.join(workdir, sp + ".makeblastdb.log")):
                h.die("makeblastdb failed for %s (see work/%s.makeblastdb.log)" % (sp, sp))
            cmd = ["blastp", "-query", qfa, "-db", db, "-evalue", "%g" % evalue, "-seg", "yes" if cfg.get("seg") else "no",
                   "-max_target_seqs", str(int(cfg["max_target_seqs"])), "-num_threads", str(a.threads),
                   "-outfmt", "7 " + OUTFMT, "-out", bl]
            if h.run_command(cmd, os.path.join(workdir, sp + ".blastp.log")):
                h.die("blastp failed for %s (see work/%s.blastp.log)" % (sp, sp))
        try:
            confirmed, unconfirmed, rows = h.read_blast_tabular(bl, len(OUTFMT.split()), OUTFMT)
        except h.HepError as e:
            h.die(str(e))
        d = decide(best_paralogs(rows, cand[sp], evalue), min_pident, min_qcov)
        done = {h.canonical_id(q) for q in confirmed}
        for g in genes:
            if g not in done:
                dec[(sp, g)] = (NOT_RUN, None, None, None, None)
            else:
                dec[(sp, g)] = d.get(g, (PASS, None, None, None, None))
        h.log("%s: %d genes, %d removed" % (sp, len(genes), sum(dec[(sp, g)][0] == REMOVED_PARALOG for g in genes)))
    nr = sum(v[0] == NOT_RUN for v in dec.values())
    if nr:
        h.warn("%d gene(s) have no finished record in the BLASTp output: NOT_RUN, they do not pass" % nr)
    rows_out = [["species", "gene_id", "verdict", "paralog", "pident", "qcov_pct", "evalue"]]
    for (sp, g), v in sorted(dec.items()):
        rows_out.append([sp, g, v[0], v[1] or "NA", "NA" if v[2] is None else "%.1f" % v[2],
                         "NA" if v[3] is None else "%.1f" % v[3], "NA" if v[4] is None else "%.2g" % v[4]])
    with open(os.path.join(a.output, "paralog_rule.tsv"), "w", encoding="utf-8", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows_out)
    keep = {k for k, v in dec.items() if v[0] == PASS}
    with open(os.path.join(a.output, "passed_ids.tsv"), "w", encoding="utf-8") as fh:
        fh.write("species\tgene_id\n")
        for sp, g in sorted(keep):
            fh.write("%s\t%s\n" % (sp, g))
    pdir = os.path.join(a.output, "passed")
    os.makedirs(pdir, exist_ok=True)
    for sp, path in ctx["inputs"]:
        with open(os.path.join(pdir, sp + ".faa"), "w", encoding="utf-8") as fh:
            for rid, header, seq in h.read_fasta(path):
                if (sp, h.canonical_id(rid)) in keep:
                    h.write_fasta_record(fh, header, seq)
    with open(os.path.join(a.output, "summary.tsv"), "w", encoding="utf-8") as fh:
        fh.write("species\tjudged\t%s\t%s\t%s\n" % (PASS, REMOVED_PARALOG, NOT_RUN))
        for sp in species:
            v = [dec[k][0] for k in dec if k[0] == sp]
            fh.write("%s\t%d\t%d\t%d\t%d\n" % (sp, len(v), v.count(PASS), v.count(REMOVED_PARALOG), v.count(NOT_RUN)))
    info = h.run_info_base(ctx)
    info.update(module=MODULE, pipeline=h.PIPELINE_NAME, version=h.__version__, tools=tools,
                passed_from=os.path.abspath(a.passed), proteomes={sp: os.path.abspath(prot[sp]) for sp in species},
                n_judged=len(dec), n_removed_paralog=sum(v[0] == REMOVED_PARALOG for v in dec.values()), n_pass=len(keep),
                n_not_run=nr)
    with open(os.path.join(a.output, "run.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
    if not a.keep_work and not a.parse_only:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)
    h.log("paralog-rule: %d judged, %d removed (paralog >= %g%% identity over >= %g%% of length), %d pass"
          % (len(dec), info["n_removed_paralog"], min_pident, min_qcov, len(keep)))
    if nr:
        h.note_incomplete(nr)


if __name__ == "__main__":
    main()
    sys.exit(h.exit_status())
