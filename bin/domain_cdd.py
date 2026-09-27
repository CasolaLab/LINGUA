#!/usr/bin/env python3
"""
domain_cdd.py - domain search with RPS-BLAST against the NCBI Conserved Domain Database (CDD).

Question answered: does the candidate gene match a known conserved domain model?

A CDD hit counts against a candidate when it passes the E-value and covers enough of the domain model
(min_domain_cov, default 50% of the model) and of the query (min_qcov, default 0). Every hit's CDD accession and
source database (Pfam, SMART, COG, ...) is written to hits.tsv. Sources listed in ignore_sources do not count.

Two ways to use it:
  run:         python3 domain_cdd.py -i proteins/ -o out/ --db /path/to/Cdd
  parse only:  python3 domain_cdd.py -i proteins/ -o out/ --parse-only earlier_results/
               (re-filters existing RPS-BLAST tabular results, e.g. when the search was run elsewhere)

Python 3.11+, standard library only. Needs `rpsblast` (BLAST+) on PATH for a run, not for --parse-only.
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402

MODULE = "domain-cdd"
DESCRIPTION = __doc__.strip().split("\n\n")[0]

# The exact output format this module asks rpsblast for. Parsing requires this column order.
OUTFMT_FIELDS = ("qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
                 "evalue bitscore qlen slen stitle")
NCOLS = len(OUTFMT_FIELDS.split())

ACCESSION_RE = re.compile(r"^\s*([^,\s]+)\s*,")                   # 'pfam00069, Pkinase, ...' in the subject title
SOURCE_RE = re.compile(r"^([A-Za-z]+)")                           # letters at the start of an accession


# --------------------------------------------------------------------------
# Parsing RPS-BLAST tabular output (outfmt 7 with the fields above)
# --------------------------------------------------------------------------
def source_of(accession):
    """CDD source database from an accession: pfam00069 -> PFAM, cd00001 -> CD, COG0515 -> COG, smart00220 -> SMART."""
    m = SOURCE_RE.match(accession or "")
    return m.group(1).upper() if m else "NA"


def parse_rps_output(path):
    """
    Parse an RPS-BLAST outfmt-7 file. Returns (confirmed, unconfirmed, rows):
      confirmed   - raw query IDs proven complete
      unconfirmed - raw query IDs that are NOT proven complete (a truncated or killed run, or a count mismatch);
                    the caller must treat these genes as NOT_RUN, never as "no hit"
      rows        - one dict per alignment line

    Completeness is proven by BLAST+'s own footer, checked against the number of queries seen; see
    hep_common.read_blast_tabular for exactly how.
    """
    confirmed, unconfirmed, raw_rows = h.read_blast_tabular(path, NCOLS, OUTFMT_FIELDS)
    rows = []
    for lineno, cols in enumerate(raw_rows, 1):
        try:
            rows.append({
                "qid": cols[0], "sid": cols[1], "pident": float(cols[2]),
                "qstart": int(cols[6]), "qend": int(cols[7]),
                "sstart": int(cols[8]), "send": int(cols[9]),
                "evalue": float(cols[10]), "bitscore": float(cols[11]),
                "qlen": int(cols[12]), "slen": int(cols[13]), "stitle": cols[14],
            })
        except ValueError as e:
            raise h.HepError("%s: cannot read a number in an alignment line for %s (%s)" % (path, cols[0], e))
    return confirmed, unconfirmed, rows


_union_len = h.union_length   # the shared implementation; kept under this name for the tests


# --------------------------------------------------------------------------
# Deciding what counts
# --------------------------------------------------------------------------
def evaluate(rows, cfg):
    """
    Group alignment lines into one hit per (gene, CDD model), apply the rules, and return
    (hit_rows, counted_genes, ignored_only_genes).
      - the E-value is applied to each alignment line;
      - domain coverage = union of the model positions covered by the alignments that passed, as a percent of the
        model length; query coverage likewise for the query;
      - a hit counts if it passes the E-value, both coverage minimums, and its source is not in ignore_sources.
    """
    evalue_max = float(cfg["evalue"])
    min_dom = float(cfg["min_domain_cov"] or 0)
    min_q = float(cfg["min_qcov"] or 0)
    ignored = {s.upper() for s in h.as_list(cfg.get("ignore_sources"))}

    groups = {}
    for r in rows:
        groups.setdefault((h.canonical_id(r["qid"]), r["sid"]), []).append(r)

    hit_rows, counted, qualifying_ignored = [], set(), set()
    for (gene, _sid), rs in sorted(groups.items()):
        passing = [r for r in rs if r["evalue"] <= evalue_max]
        best = min(rs, key=lambda r: (r["evalue"], -r["bitscore"]))
        m = ACCESSION_RE.match(best["stitle"])
        accession = m.group(1) if m else best["sid"]
        source = source_of(accession) if m else "NA"
        row = {"gene_id": gene, "target": accession, "source": source, "evalue": best["evalue"],
               "bitscore": best["bitscore"], "pident": best["pident"], "counts_as_hit": False}
        if not passing:
            row["note"] = "E-value above %g" % evalue_max
        else:
            qcov = 100.0 * _union_len([(r["qstart"], r["qend"]) for r in passing]) / max(1, passing[0]["qlen"])
            tcov = 100.0 * _union_len([(r["sstart"], r["send"]) for r in passing]) / max(1, passing[0]["slen"])
            row["qcov"], row["tcov"] = round(qcov, 1), round(tcov, 1)
            best_passing = min(passing, key=lambda r: (r["evalue"], -r["bitscore"]))
            row.update(evalue=best_passing["evalue"], bitscore=best_passing["bitscore"],
                       pident=best_passing["pident"])
            if tcov < min_dom:
                row["note"] = "domain coverage %.1f%% below %g%%" % (tcov, min_dom)
            elif qcov < min_q:
                row["note"] = "query coverage %.1f%% below %g%%" % (qcov, min_q)
            elif source in ignored:
                row["note"] = "ignored source %s" % source
                qualifying_ignored.add(gene)
            else:
                row["counts_as_hit"] = True
                counted.add(gene)
        hit_rows.append(row)
    return hit_rows, counted, qualifying_ignored - counted


# --------------------------------------------------------------------------
# Running rpsblast
# --------------------------------------------------------------------------
def cdd_release(db):
    """CDD release text if a cdd.info file sits next to the database, else None."""
    info = os.path.join(os.path.dirname(os.path.abspath(db)), "cdd.info")
    try:
        with open(info, "r", encoding="utf-8") as fh:
            return fh.readline().strip() or None
    except OSError:
        return None


def rps_command(chunk, db, out_tsv, cfg, threads):
    cmd = ["rpsblast", "-query", chunk, "-db", db, "-evalue", "%g" % float(cfg["evalue"]),
           "-outfmt", "7 " + OUTFMT_FIELDS, "-num_threads", str(int(threads)), "-out", out_tsv]
    if cfg.get("comp_based_stats") is not None:
        cmd += ["-comp_based_stats", str(int(cfg["comp_based_stats"]))]
    return cmd


def run_search(species, fasta, cfg, args, out):
    """
    Run rpsblast on the chunks of one species. Returns (raw_path, not_run_ids): the merged raw result file, and the
    IDs of genes in chunks that failed or did not finish (they become NOT_RUN, never NO_HIT).
    """
    chunk_dir = os.path.join(out.workdir, species, "chunks")
    res_dir = os.path.join(out.workdir, species, "results")
    log_dir = os.path.join(out.workdir, species, "logs")
    for d in (res_dir, log_dir):
        os.makedirs(d, exist_ok=True)
    chunks = h.chunk_fasta(fasta, chunk_dir, int(cfg["chunk_size"]))
    by_name = {os.path.splitext(os.path.basename(c))[0]: c for c in chunks}

    def job(name):
        tsv = os.path.join(res_dir, name + ".tsv")
        rc = h.run_command(rps_command(by_name[name], args.db, tsv, cfg, args.threads),
                           os.path.join(log_dir, name + ".log"))
        h.check_rc(rc, "rpsblast on %s/%s" % (species, name))
        confirmed, unconfirmed, _rows = parse_rps_output(tsv)
        expected = {rid for rid, _h, _s in h.read_fasta(by_name[name])}
        if unconfirmed or not expected <= confirmed:
            raise h.CommandError("rpsblast output for %s/%s is incomplete (no completion footer)" % (species, name))

    res = h.run_jobs(sorted(by_name), job, workers=int(cfg["parallel_jobs"]),
                     marker_dir=os.path.join(out.workdir, species, "markers"), resume=args.resume)
    not_run = set()
    for name in res["failed"]:
        h.warn("%s/%s failed: %s" % (species, name, res["failed"][name]))
        not_run |= {h.canonical_id(rid) for rid, _h, _s in h.read_fasta(by_name[name])}
    raw_dir = os.path.join(out.outdir, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    raw_path = os.path.join(raw_dir, species + ".rpsblast.tsv")
    with open(raw_path, "w", encoding="utf-8") as merged:
        for name in sorted(by_name):
            tsv = os.path.join(res_dir, name + ".tsv")
            if name not in res["failed"] and os.path.exists(tsv):
                with open(tsv, "r", encoding="utf-8") as fh:
                    merged.write(fh.read())
    return raw_path, not_run


def find_existing(parse_only, species):
    """Path of an earlier result for this species (<species>.rpsblast.tsv in a folder, or the file itself)."""
    if os.path.isfile(parse_only):
        return parse_only
    p = os.path.join(parse_only, species + ".rpsblast.tsv")
    return p if os.path.isfile(p) else None


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def add_args(parser):
    g = parser.add_argument_group("domain-cdd options")
    g.add_argument("--db", default=None,
                   help="path prefix of the formatted CDD database for rpsblast (e.g. /db/cdd/Cdd); "
                        "required unless --parse-only")
    g.add_argument("--parse-only", default=None, metavar="PATH",
                   help="do not run rpsblast: read existing results, a folder holding <species>.rpsblast.tsv "
                        "(or one file for a single species), produced with outfmt '7 " + OUTFMT_FIELDS + "'")


def main(argv=None):
    ctx = h.start_module(MODULE, DESCRIPTION, argv, extra=add_args)
    args, cfg = ctx["args"], ctx["cfg"]
    if not args.parse_only:
        if not args.db:
            h.die("--db is required (or use --parse-only to read existing results)")
        if not h.find_tool("rpsblast"):
            h.die("rpsblast not found in PATH (install BLAST+; see docs/setup.md)")
        if not any(f.startswith(os.path.basename(args.db) + ".") for f in
                   os.listdir(os.path.dirname(os.path.abspath(args.db)) or ".")):
            h.die("no database files found for prefix %s" % args.db)

    ignored = {s.upper() for s in h.as_list(cfg.get("ignore_sources"))}
    try:
        out = h.ModuleOutput(args.output, MODULE, args.label, resume=args.resume, force=args.force)
    except h.HepError as e:
        h.die(str(e))
    with out:
        for species, fasta in ctx["inputs"]:
            h.log("%s: %s" % (species, "parsing existing results" if args.parse_only else "running rpsblast"))
            try:
                all_ids = h.load_fasta_ids(fasta)
            except h.HepError as e:
                h.die(str(e))
            not_run, hit_rows, counted, excluded = set(all_ids), [], set(), set()
            if args.parse_only:
                raw = find_existing(args.parse_only, species)
                if raw is None:
                    h.warn("%s: no earlier result file found; all its genes are NOT_RUN" % species)
                    seen = set()
                else:
                    confirmed, unconfirmed, rows = parse_rps_output(raw)
                    seen = {h.canonical_id(q) for q in confirmed}
                    if unconfirmed:
                        h.warn("%s: %d gene(s) are in an unfinished result block and are NOT_RUN"
                               % (species, len(unconfirmed)))
                    hit_rows, counted, excluded = evaluate(
                        [r for r in rows if r["qid"] in confirmed], cfg)
                not_run = set(all_ids) - seen
            else:
                raw, failed_ids = run_search(species, fasta, cfg, args, out)
                confirmed, unconfirmed, rows = parse_rps_output(raw)
                seen = {h.canonical_id(q) for q in confirmed}
                hit_rows, counted, excluded = evaluate(rows, cfg)
                not_run = (set(all_ids) - seen) | failed_ids
            for r in hit_rows:
                r["species"] = species
            statuses = h.classify_genes(all_ids, counted=counted, excluded=excluded, not_run=not_run)
            tally = out.add_species(species, fasta, hit_rows, statuses)
            h.log("%s: %d genes, %d passed, %d with a counted CDD hit, %d NOT_RUN"
                  % (species, tally["input"], tally["passed"], tally[h.HIT], tally[h.NOT_RUN]))

        tools = []
        if args.parse_only:
            tools.append({"name": "RPS-BLAST (BLAST+), results parsed from an earlier run", "version": None})
        else:
            tools.append({"name": "RPS-BLAST (BLAST+)",
                          "version": h.tool_version(["rpsblast", "-version"], r"\d+\.\d+\.\d+\+?")})
            tools.append({"name": "NCBI Conserved Domain Database (CDD)", "version": cdd_release(args.db)})
        info = h.run_info_base(ctx)
        info.update(tools=tools, database=(None if args.parse_only else os.path.abspath(args.db)),
                    parsed_from=(os.path.abspath(args.parse_only) if args.parse_only else None),
                    ignore_sources=sorted(ignored))
        out.close(info, keep_work=args.keep_work, quiet=args.quiet)


if __name__ == "__main__":
    main()
    sys.exit(h.exit_status())
