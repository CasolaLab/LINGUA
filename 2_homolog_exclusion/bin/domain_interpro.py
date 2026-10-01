#!/usr/bin/env python3
"""
domain_interpro.py - domain search with InterProScan.

Question answered: does the candidate gene carry a known protein domain or family?

Every match is sorted by the analysis that produced it:
  - domain and family databases (Pfam, SMART, CDD, PRINTS, ...) are domain evidence and count against a candidate;
  - MobiDB-Lite, Coils and PROSITE patterns find short patterns or low-complexity regions, not domains: they never count
    (a candidate with only these is reported EXCLUDED_ONLY and passes);
  - AntiFam flags a spurious gene model, not a domain (a candidate with only that is reported SPURIOUS, and passes
    unless antifam_action=remove);
  - domains listed in ignore_domains (known lineage-specific domains for your clade) do not count either.

InterProScan is always run with -dp: no pre-calculated lookup. Everything is computed on your machine, nothing is sent to
EBI, and no internet connection is needed.

Two ways to use it:
  run:         python3 domain_interpro.py -i proteins/ -o out/
  parse only:  python3 domain_interpro.py -i proteins/ -o out/ --parse-only earlier_results/ [--assume-complete]

Python 3.11+, standard library only. Needs `interproscan.sh` on PATH (or --iprscan-bin) for a run, not for --parse-only.
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402

MODULE = "domain-interpro"
DESCRIPTION = __doc__.strip().split("\n\n")[0]

NCOLS = 15   # InterProScan's TSV: 15 tab-separated columns, no header line
# Analyses whose column 9 is an E-value. For the others it is a score (Hamap, ProSiteProfiles) or '-' (Coils,
# MobiDB-Lite, ProSitePatterns), so max_evalue is never applied to those.
E_VALUE_ANALYSES = {"pfam", "prints", "smart", "cdd", "pirsf", "pirsr", "sfld", "antifam", "gene3d", "superfamily",
                    "ncbifam", "funfam", "panther"}
DONE_RE = re.compile(r"100% done")
ANALYSES_RE = re.compile(r"^\[([^\]]+)\]\s*$")          # the line that lists the analyses and their versions
BAD_APPL_RE = re.compile(r"(Invalid input specified for -appl[^\n]*\n[^\n]*)")


def norm(name):
    """Compare names ignoring case and punctuation: MobiDB-Lite == MobiDBLite == mobidblite."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


# --------------------------------------------------------------------------
# Reading InterProScan's TSV
# --------------------------------------------------------------------------
def parse_tsv(path):
    """Return the rows of an InterProScan TSV file (15 columns, no header). A missing or empty file has no rows."""
    rows = []
    if not os.path.isfile(path):
        return rows
    with open(path, "r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            c = line.split("\t")
            if len(c) != NCOLS:
                raise h.HepError("%s line %d: expected %d tab-separated columns (InterProScan TSV, which has no header "
                                 "line), got %d" % (path, lineno, NCOLS, len(c)))
            try:
                rows.append({"qid": c[0], "length": int(c[2]), "analysis": c[3], "sig": c[4], "desc": c[5],
                             "start": int(c[6]), "stop": int(c[7]), "score": c[8], "status": c[9],
                             "ipr": c[11], "ipr_desc": c[12]})
            except ValueError as e:
                raise h.HepError("%s line %d: cannot read a number (%s); is this an InterProScan TSV?" % (path, lineno, e))
    return rows


def _float(text):
    try:
        return float(text)
    except ValueError:
        return None


# --------------------------------------------------------------------------
# Deciding what counts
# --------------------------------------------------------------------------
def evaluate(rows, cfg):
    """
    Sort every match into: counts, excluded, spurious, or ignored. Returns
    (hit_rows, counted_genes, excluded_genes, spurious_genes).
    Decided by ANALYSIS, never by whether a match has an InterPro accession (real domain matches from CDD, SFLD, PIRSF
    and others often have none).
    """
    excluded_an = {norm(a) for a in h.as_list(cfg.get("exclude_analyses"))}
    spurious_an = {norm(a) for a in h.as_list(cfg.get("spurious_analyses"))}
    ignored_dom = {norm(d) for d in h.as_list(cfg.get("ignore_domains"))}
    max_e = cfg.get("max_evalue")
    max_e = float(max_e) if max_e not in (None, "") else None

    hit_rows, counted, excluded, spurious = [], set(), set(), set()
    for r in rows:
        gene = h.canonical_id(r["qid"])
        a = norm(r["analysis"])
        value = _float(r["score"])
        row = {"gene_id": gene, "target": r["sig"], "source": r["analysis"], "counts_as_hit": False,
               "qcov": round(100.0 * (r["stop"] - r["start"] + 1) / r["length"], 1) if r["length"] > 0 else None}
        if value is not None:
            row["evalue" if a in E_VALUE_ANALYSES else "bitscore"] = value
        if a in spurious_an:
            row["note"] = "spurious gene-model signal (%s), not a domain" % r["analysis"]
            spurious.add(gene)
        elif a in excluded_an:
            row["note"] = "not domain evidence (%s)" % r["analysis"]
            excluded.add(gene)
        elif r["status"] != "T":
            row["note"] = "match status '%s', not a confirmed match" % r["status"]
        elif ignored_dom & {norm(r["sig"]), norm(r["ipr"]), norm(r["desc"])} - {""}:
            row["note"] = "ignored domain (ignore_domains)"
            excluded.add(gene)
        elif max_e is not None and a in E_VALUE_ANALYSES and value is not None and value > max_e:
            row["note"] = "E-value above %g" % max_e
        else:
            row["counts_as_hit"] = True
            row["note"] = ("InterPro %s" % r["ipr"]) if r["ipr"] not in ("", "-") else "no InterPro entry"
            counted.add(gene)
        hit_rows.append(row)
    return hit_rows, counted, excluded, spurious


# --------------------------------------------------------------------------
# Evidence that a run finished
# --------------------------------------------------------------------------
def log_says_complete(path):
    """
    InterProScan's TSV lists only the proteins that have matches, so a protein with no rows cannot be told from one that
    was never searched. Proof that the run finished comes from elsewhere: InterProScan's own log ('100% done ...'),
    or the status file this module writes ('hep_status=complete').
    """
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return False
    return bool(DONE_RE.search(text)) or "hep_status=complete" in text


def analyses_from_log(text):
    """['Pfam-38.2', 'CDD-3.21', ...] from InterProScan's 'Running the following analyses' line, if present."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if "Running the following analyses" in line:
            for nxt in lines[i:i + 3]:
                m = ANALYSES_RE.match(nxt.strip())
                if m:
                    return [x.strip() for x in m.group(1).split(",") if x.strip()]
    return []


# --------------------------------------------------------------------------
# Running InterProScan
# --------------------------------------------------------------------------
def strip_stops(path):
    """Remove '*' (stop symbols) from a chunk file in place; InterProScan rejects them. Returns how many sequences had any."""
    records = list(h.read_fasta(path))
    n = 0
    with open(path, "w", encoding="utf-8") as fh:
        for _rid, header, seq in records:
            if "*" in seq:
                n += 1
                seq = seq.replace("*", "")
            h.write_fasta_record(fh, header, seq)
    return n


def ips_command(iprscan, chunk, out_tsv, tmp_dir, cfg, threads):
    """InterProScan command line. -dp (no pre-calculated lookup) is always present: it keeps the run local and offline."""
    apps = h.as_list(cfg.get("applications"))
    cmd = [iprscan, "-i", chunk, "-f", "tsv", "-o", out_tsv, "-dp", "-cpu", str(int(threads)), "-T", tmp_dir]
    if apps:
        cmd += ["-appl", ",".join(apps)]
    return cmd


def run_search(species, fasta, cfg, args, out, iprscan):
    """
    Run InterProScan on the chunks of one species. Returns (raw_tsv, complete, failed_ids, analyses):
      complete    - every chunk finished (exit code 0 and InterProScan's '100% done' line);
      failed_ids  - IDs of genes in chunks that did not finish (they become NOT_RUN, never NO_HIT).
    """
    base = os.path.join(out.workdir, species)
    chunk_dir, res_dir, log_dir, tmp_root = (os.path.join(base, d) for d in ("chunks", "results", "logs", "tmp"))
    for d in (res_dir, log_dir, tmp_root):
        os.makedirs(d, exist_ok=True)
    chunks = h.chunk_fasta(fasta, chunk_dir, int(cfg["chunk_size"]))
    stripped = sum(strip_stops(c) for c in chunks)
    if stripped:
        h.warn("%s: '*' stop symbols were removed from %d sequence(s) before searching (InterProScan rejects them)"
               % (species, stripped))
    by_name = {os.path.splitext(os.path.basename(c))[0]: c for c in chunks}

    def job(name):
        tsv = os.path.join(res_dir, name + ".tsv")
        log = os.path.join(log_dir, name + ".log")
        tmp = os.path.join(tmp_root, name)      # a separate temp folder per job: parallel jobs must not share one
        os.makedirs(tmp, exist_ok=True)
        if os.path.exists(tsv):
            os.remove(tsv)
        if os.path.exists(log):
            os.remove(log)
        rc = h.run_command(ips_command(iprscan, by_name[name], tsv, tmp, cfg, args.threads), log)
        h.check_rc(rc, "InterProScan on %s/%s" % (species, name))
        if not log_says_complete(log):
            raise h.CommandError("InterProScan on %s/%s ended without its completion line ('100%% done')"
                                 % (species, name))

    res = h.run_jobs(sorted(by_name), job, workers=int(cfg["parallel_jobs"]),
                     marker_dir=os.path.join(base, "markers"), resume=args.resume)
    failed_ids = set()
    for name, msg in res["failed"].items():
        h.warn("%s/%s failed: %s" % (species, name, msg))
        failed_ids |= {h.canonical_id(rid) for rid, _h, _s in h.read_fasta(by_name[name])}
    # An invalid analysis name is the user's mistake and fails every chunk the same way: say so once, clearly.
    for name in res["failed"]:
        try:
            with open(os.path.join(log_dir, name + ".log"), "r", encoding="utf-8", errors="replace") as fh:
                m = BAD_APPL_RE.search(fh.read())
        except OSError:
            m = None
        if m:
            h.die("InterProScan rejected the analyses list: %s\nUse names your installation accepts (see docs/domain-interpro.md)."
                  % " ".join(m.group(1).split()))

    raw_dir = os.path.join(out.outdir, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    raw_tsv = os.path.join(raw_dir, species + ".interproscan.tsv")
    analyses = []
    with open(raw_tsv, "w", encoding="utf-8") as merged:
        for name in sorted(by_name):
            if name in res["failed"]:
                continue
            tsv = os.path.join(res_dir, name + ".tsv")
            if os.path.exists(tsv):
                with open(tsv, "r", encoding="utf-8") as fh:
                    text = fh.read()
                if text and not text.endswith("\n"):
                    text += "\n"
                merged.write(text)
            try:
                with open(os.path.join(log_dir, name + ".log"), "r", encoding="utf-8", errors="replace") as fh:
                    analyses = analyses or analyses_from_log(fh.read())
            except OSError:
                pass
    complete = not res["failed"]
    with open(os.path.join(raw_dir, species + ".interproscan.status"), "w", encoding="utf-8") as fh:
        fh.write("hep_status=%s\nchunks_finished=%d of %d\n"
                 % ("complete" if complete else "incomplete", len(by_name) - len(res["failed"]), len(by_name)))
    return raw_tsv, complete, failed_ids, analyses


def find_existing(parse_only, species):
    """(tsv path, evidence-of-completion path or None) for a species from an earlier run."""
    if os.path.isfile(parse_only):
        tsv = parse_only
        stem = tsv[:-len(".tsv")] if tsv.endswith(".tsv") else tsv
    else:
        tsv = os.path.join(parse_only, species + ".interproscan.tsv")
        stem = os.path.join(parse_only, species + ".interproscan")
    if not os.path.isfile(tsv):
        return None, None
    for ext in (".status", ".log"):
        if os.path.isfile(stem + ext) and log_says_complete(stem + ext):
            return tsv, stem + ext
    return tsv, None


def tool_version(iprscan):
    line = h.tool_version([iprscan, "--version"], r"InterProScan version \S+")
    return line.replace("InterProScan version ", "") if line != "unknown" else None


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def add_args(parser):
    g = parser.add_argument_group("domain-interpro options")
    g.add_argument("--iprscan-bin", default=None,
                   help="path to interproscan.sh (default: interproscan.sh or iprscan found on PATH)")
    g.add_argument("--parse-only", default=None, metavar="PATH",
                   help="do not run InterProScan: read existing results, a folder holding <species>.interproscan.tsv "
                        "(or one file for a single species)")
    g.add_argument("--assume-complete", action="store_true",
                   help="with --parse-only: state that the earlier run finished (no completion evidence is needed). "
                        "Without it, a protein with no rows is NOT_RUN, because InterProScan's TSV lists only proteins "
                        "that have matches")


def main(argv=None):
    ctx = h.start_module(MODULE, DESCRIPTION, argv, extra=add_args)
    args, cfg = ctx["args"], ctx["cfg"]
    if cfg.get("antifam_action") not in ("flag", "remove"):
        h.die("antifam_action must be 'flag' or 'remove', not %r" % cfg.get("antifam_action"))
    iprscan = None
    if not args.parse_only:
        iprscan = args.iprscan_bin or h.find_tool("interproscan.sh", "iprscan")
        if not iprscan or not os.path.isfile(iprscan):
            h.die("interproscan.sh not found (put it on PATH or give --iprscan-bin; see docs/setup.md)")
    try:
        out = h.ModuleOutput(args.output, MODULE, args.label, resume=args.resume, force=args.force)
    except h.HepError as e:
        h.die(str(e))

    analyses_used = []
    with out:
        for species, fasta in ctx["inputs"]:
            h.log("%s: %s" % (species, "parsing existing results" if args.parse_only else "running InterProScan"))
            try:
                all_ids = h.load_fasta_ids(fasta)
            except h.HepError as e:
                h.die(str(e))
            failed_ids, evidence = set(), None
            if args.parse_only:
                tsv, evidence = find_existing(args.parse_only, species)
                complete = bool(evidence) or args.assume_complete
                if tsv is None:
                    h.warn("%s: no earlier result file found; all its genes are NOT_RUN" % species)
                    rows, complete = [], False
                else:
                    rows = parse_tsv(tsv)
                    if not complete:
                        h.warn("%s: nothing shows that the earlier run finished (no InterProScan log or status file, "
                               "no --assume-complete), so only genes with a counted domain are decided; the rest are "
                               "NOT_RUN" % species)
            else:
                tsv, _all_chunks_ok, failed_ids, analyses = run_search(species, fasta, cfg, args, out, iprscan)
                analyses_used = analyses_used or analyses
                rows = parse_tsv(tsv)
                # In a run, each chunk is known to have finished (exit code 0 and '100% done') or not. Genes of the
                # chunks that finished and have no rows really have no matches; genes of failed chunks are NOT_RUN.
                complete = True

            ids = set(all_ids)
            unknown = {h.canonical_id(r["qid"]) for r in rows} - ids
            if unknown:
                h.warn("%s: %d protein(s) in the results are not in the input FASTA and were ignored (e.g. %s)"
                       % (species, len(unknown), ", ".join(sorted(unknown)[:3])))
            rows = [r for r in rows if h.canonical_id(r["qid"]) in ids]
            hit_rows, counted, excluded, spurious = evaluate(rows, cfg)
            for r in hit_rows:
                r["species"] = species
            if complete:
                not_run = failed_ids
            else:
                # Without proof that the run finished, only a counted domain is certain: any other verdict could
                # change if more of the run's results were missing.
                not_run = (ids - counted) | failed_ids
                excluded, spurious = excluded - not_run, spurious - not_run
            statuses = h.classify_genes(all_ids, counted=counted, excluded=excluded, spurious=spurious, not_run=not_run)
            tally = out.add_species(species, fasta, hit_rows, statuses, cfg["antifam_action"])
            h.log("%s: %d genes, %d passed, %d with a counted domain, %d NOT_RUN"
                  % (species, tally["input"], tally["passed"], tally[h.HIT], tally[h.NOT_RUN]))

        tools = []
        if args.parse_only:
            tools.append({"name": "InterProScan, results parsed from an earlier run", "version": None})
        else:
            tools.append({"name": "InterProScan", "version": tool_version(iprscan)})
            tools.extend({"name": a.rsplit("-", 1)[0], "version": a.rsplit("-", 1)[1]} if "-" in a
                         else {"name": a, "version": None} for a in analyses_used)
        info = h.run_info_base(ctx)
        info.update(tools=tools, analyses_run=analyses_used, applications=h.as_list(cfg.get("applications")),
                    lookup="disabled (-dp): computed locally, nothing sent to EBI",
                    parsed_from=(os.path.abspath(args.parse_only) if args.parse_only else None),
                    assumed_complete=bool(args.parse_only and args.assume_complete))
        out.close(info, keep_work=args.keep_work, quiet=args.quiet)


if __name__ == "__main__":
    main()
    sys.exit(h.exit_status())
