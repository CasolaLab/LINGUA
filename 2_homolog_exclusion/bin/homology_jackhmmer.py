#!/usr/bin/env python3
"""
homology_jackhmmer.py - sensitive protein homology search with jackhmmer (HMMER 3).

Question answered: does the candidate gene have a remote homolog in the user's database? jackhmmer builds a profile from
each candidate and searches again (iterations), so it finds relatives BLAST+ misses.

The rule is the same as homology-blast: ANY hit counts against a candidate (subject to the E-value and the coverage
settings). The database is one species-tagged FASTA (IDs ending __Species_Name; make-db builds it) and must not contain
the analysis species or close relatives: hits from the analysis species (the input file names) and from ignore_taxa
species are reported but not counted (status IN_PHYLOGENY_ONLY).

Several databases can be searched in one command (--db given more than once), each as its own labelled run, or one after
another with --sequential.

Nothing here uses the network. Python 3.11+, standard library only. Needs jackhmmer on PATH for a run, not for --parse-only.

  run:         python3 homology_jackhmmer.py -i proteins/ -o out/ --db db/mosses.faa
  parse only:  python3 homology_jackhmmer.py -i proteins/ -o out/ --parse-only earlier_out/raw
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402
import homology_blast as b  # noqa: E402  (shared: database problems, species map, inventory warnings)

MODULE = "homology-jackhmmer"
DESCRIPTION = __doc__.strip().split("\n\n")[0]
DbProblem = b.DbProblem
OK_FOOTER = "# [ok]"      # HMMER writes this as the last line of a table file only when the search finished
RAW_SUFFIX = ".jackhmmer.domtbl"
SEARCHED_SUFFIX = ".jackhmmer.searched"


# --------------------------------------------------------------------------
# Reading jackhmmer output (--domtblout)
# --------------------------------------------------------------------------
def parse_domtbl(path):
    """
    (complete, rows) for a jackhmmer --domtblout file. `complete` is True only if the last line is HMMER's '# [ok]'.
    jackhmmer writes one row per domain and nothing at all for a query without hits, so a finished file is the only
    proof that a query with no rows was searched. Rows are dicts; 'hmm' coordinates are positions on the QUERY (the
    query is the profile), 'ali' coordinates are positions on the target.
    """
    rows, last = [], ""
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            last = line.strip()
            if line.startswith("#"):
                continue
            c = line.split(None, 22)
            if len(c) < 22:
                raise h.HepError("%s: a domain line has %d columns, expected 22 or more: %r" % (path, len(c), line[:80]))
            try:
                rows.append({"tid": c[0], "tlen": int(c[2]), "qid": c[3], "qlen": int(c[5]), "evalue": float(c[6]),
                             "bitscore": float(c[7]), "ievalue": float(c[12]), "qstart": int(c[15]), "qend": int(c[16]),
                             "sstart": int(c[17]), "send": int(c[18])})
            except ValueError as e:
                raise h.HepError("%s: cannot read a number in a domain line for %s (%s)" % (path, c[3], e))
    return last == OK_FOOTER, rows


def read_searched(path):
    """IDs of the queries a hep run recorded as searched (one per line), or None if there is no such file."""
    side = re.sub(re.escape(RAW_SUFFIX) + r"$", SEARCHED_SUFFIX, path)
    if side == path or not os.path.isfile(side):
        return None
    with open(side, "r", encoding="utf-8") as fh:
        return {line.strip() for line in fh if line.strip()}


# --------------------------------------------------------------------------
# Deciding what counts
# --------------------------------------------------------------------------
def evaluate(rows, cfg, matcher, tmap, label):
    """
    Group domain rows into one hit per (gene, database sequence) and apply the rules. Returns
    (hit_rows, counted_genes, in_phylogeny_genes).
      - the E-value is the whole-sequence E-value (the same on every row of a hit), against `evalue`;
      - coverage uses only the domains whose own E-value also passes `evalue`: query coverage = percent of the query
        covered by the union of their query intervals; target coverage likewise on the target;
      - a hit counts if it passes the E-value, min_qcov and min_tcov;
      - a hit that qualifies but comes from an analysis species or an ignore_taxa species does not count when
        ignore_in_phylogeny is on; its gene becomes IN_PHYLOGENY_ONLY unless another hit counts.
    """
    evalue_max = float(cfg["evalue"])
    min_q = float(cfg["min_qcov"] or 0)
    min_t = float(cfg["min_tcov"] or 0)
    exempt = bool(cfg["ignore_in_phylogeny"])
    groups = {}
    for r in rows:
        groups.setdefault((h.canonical_id(r["qid"]), r["tid"]), []).append(r)

    hit_rows, counted, in_phylo = [], set(), set()
    for (gene, tid), rs in sorted(groups.items()):
        sp = b.target_species(tid, cfg, tmap)
        row = {"gene_id": gene, "target": tid, "target_species": sp, "source": label, "counts_as_hit": False,
               "evalue": rs[0]["evalue"], "bitscore": rs[0]["bitscore"], "pident": None}
        if rs[0]["evalue"] > evalue_max:
            row["note"] = "E-value above %g" % evalue_max
            hit_rows.append(row)
            continue
        good = [r for r in rs if r["ievalue"] <= evalue_max]
        qcov = 100.0 * h.union_length([(r["qstart"], r["qend"]) for r in good]) / max(1, rs[0]["qlen"])
        tcov = 100.0 * h.union_length([(r["sstart"], r["send"]) for r in good]) / max(1, rs[0]["tlen"])
        row.update(qcov=round(qcov, 1), tcov=round(tcov, 1))
        if qcov < min_q:
            row["note"] = "query coverage %.1f%% below %g%%" % (qcov, min_q)
        elif tcov < min_t:
            row["note"] = "target coverage %.1f%% below %g%%" % (tcov, min_t)
        else:
            match = matcher.match(sp) if exempt else None
            if match:
                row["note"] = "in-phylogeny hit (%s), not counted" % match
                in_phylo.add(gene)
            else:
                row["counts_as_hit"] = True
                row["note"] = None if sp else "target species unknown, counted"
                counted.add(gene)
        hit_rows.append(row)
    return hit_rows, counted, in_phylo - counted


# --------------------------------------------------------------------------
# The database (a species-tagged FASTA)
# --------------------------------------------------------------------------
def inspect_db(db):
    """The database must be a protein FASTA file. Returns {'type', 'sequences'}; raises DbProblem otherwise."""
    if not os.path.isfile(db):
        raise DbProblem("%s is not a file. jackhmmer reads a species-tagged protein FASTA directly (make-db builds one); "
                        "a BLAST database prefix cannot be used" % db)
    with open(db, "r", encoding="utf-8", errors="replace") as fh:
        first = next((ln for ln in fh if ln.strip()), "")
    if not first.startswith(">"):
        raise DbProblem("%s does not look like a FASTA file (it does not start with '>'). If it is a BLAST database file, "
                        "give jackhmmer the FASTA it was made from" % db)
    return {"type": "fasta", "sequences": None}


def fasta_inventory(db):
    """Count how many sequences carry each '__Species_Name' tag (same result shape as homology-blast)."""
    inv = {"ids": 0, "tagged": 0, "species": {}, "first_ids": []}
    with open(db, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            sid = line[1:].split(None, 1)[0] if line[1:].strip() else ""
            if not sid:
                continue
            inv["ids"] += 1
            if len(inv["first_ids"]) < 5:
                inv["first_ids"].append(sid)
            sp = h.split_id(sid)[1]
            if sp:
                inv["tagged"] += 1
                inv["species"][sp] = inv["species"].get(sp, 0) + 1
    return inv


# --------------------------------------------------------------------------
# Running jackhmmer
# --------------------------------------------------------------------------
def jackhmmer_command(chunk, db, out_tbl, cfg, threads=1):
    cmd = ["jackhmmer", "-N", str(int(cfg["iterations"])), "-E", "%g" % float(cfg["evalue"]),
           "--incE", "%g" % float(cfg["incE"]), "--noali", "--cpu", str(max(1, int(threads))), "-o", os.devnull, "--domtblout", out_tbl]
    if cfg.get("dbsize") is not None:
        cmd += ["-Z", str(int(cfg["dbsize"]))]
    return cmd + [chunk, db]


def check_settings(cfg):
    if int(cfg["iterations"]) < 1:
        raise h.HepError("iterations must be at least 1")
    for key in ("evalue", "incE"):
        if float(cfg[key]) <= 0:
            raise h.HepError("%s must be above 0" % key)
    if cfg.get("species_source") not in ("header", "map"):
        raise h.HepError("species_source must be 'header' or 'map', not %r" % cfg.get("species_source"))
    if cfg.get("species_source") == "map" and not cfg.get("taxon_map"):
        raise h.HepError("species_source=map needs taxon_map (a file with accession and species)")
    if cfg.get("match_level") not in h.MATCH_LEVELS:
        raise h.HepError("match_level must be one of %s, not %r" % (", ".join(h.MATCH_LEVELS), cfg.get("match_level")))


def run_search(species, fasta, db, cfg, args, out):
    """
    Run jackhmmer on the chunks of one species against one database. Returns (raw_path, searched_ids, not_run_ids):
    the merged result file of the chunks that finished, the IDs of the queries in those chunks (proven searched), and the
    IDs in chunks that failed or did not finish (NOT_RUN, never NO_HIT).
    """
    base = os.path.join(out.workdir, species)
    chunk_dir, res_dir, log_dir = (os.path.join(base, d) for d in ("chunks", "results", "logs"))
    for d in (res_dir, log_dir):
        os.makedirs(d, exist_ok=True)
    chunks = h.chunk_fasta(fasta, chunk_dir, int(cfg["chunk_size"]))
    by_name = {os.path.splitext(os.path.basename(c))[0]: c for c in chunks}

    def job(name):
        tbl = os.path.join(res_dir, name + ".domtbl")
        rc = h.run_command(jackhmmer_command(by_name[name], db, tbl, cfg, args.threads), os.path.join(log_dir, name + ".log"))
        h.check_rc(rc, "jackhmmer on %s/%s" % (species, name))
        if not os.path.exists(tbl) or not parse_domtbl(tbl)[0]:
            raise h.CommandError("jackhmmer output for %s/%s is incomplete (no '# [ok]' at the end)" % (species, name))

    res = h.run_jobs(sorted(by_name), job, workers=int(cfg["parallel_jobs"]),
                     marker_dir=os.path.join(base, "markers"), resume=args.resume)
    not_run, searched = set(), set()
    for name in sorted(by_name):
        ids = {h.canonical_id(rid) for rid, _h, _s in h.read_fasta(by_name[name])}
        if name in res["failed"]:
            h.warn("%s/%s failed: %s" % (species, name, res["failed"][name]))
            not_run |= ids
        else:
            searched |= ids
    raw_dir = os.path.join(out.outdir, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    raw_path = os.path.join(raw_dir, species + RAW_SUFFIX)
    with open(raw_path, "w", encoding="utf-8") as merged:
        for name in sorted(by_name):
            tbl = os.path.join(res_dir, name + ".domtbl")
            if name not in res["failed"] and os.path.exists(tbl):
                with open(tbl, "r", encoding="utf-8") as fh:
                    merged.write(fh.read())
    with open(os.path.join(raw_dir, species + SEARCHED_SUFFIX), "w", encoding="utf-8") as fh:
        fh.write("".join(i + "\n" for i in sorted(searched)))
    return raw_path, searched, not_run


def find_existing(parse_only, species):
    if os.path.isfile(parse_only):
        return parse_only
    p = os.path.join(parse_only, species + RAW_SUFFIX)
    return p if os.path.isfile(p) else None


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def add_args(parser):
    g = parser.add_argument_group("homology-jackhmmer options")
    g.add_argument("--db", action="append", default=[], metavar="FASTA",
                   help="the database: a species-tagged protein FASTA file (make-db builds one). Give it more than once to "
                        "search several databases, each as its own labelled run. Required unless --parse-only")
    g.add_argument("--skip-db-inventory", action="store_true",
                   help="do not read all the database's IDs to count species and look for your analysis species in it")
    g.add_argument("--parse-only", default=None, metavar="PATH",
                   help="do not run jackhmmer: read existing results, a folder holding <species>%s (or one file for a "
                        "single species), made with --domtblout. A file must end with HMMER's '# [ok]' line" % RAW_SUFFIX)
    g.add_argument("--sequential", action="store_true",
                   help="with several --db: search them in the order given, each one only against the genes that passed "
                        "the previous one (iterative filtering). Without it every database searches every gene. If one "
                        "database fails, the chain stops, because genes cannot pass a database that was never searched")


def run_one_database(ctx, db, label, outdir, matcher, tmap, inputs=None, chain=None):
    """One database (or one parse-only read) as a complete labelled run; removes its output folder if the database is unusable."""
    args = ctx["args"]
    try:
        out = h.ModuleOutput(outdir, MODULE, label, resume=args.resume, force=args.force)
    except h.HepError as e:
        raise DbProblem(str(e))
    try:
        _run_one_database_body(ctx, db, label, out, matcher, tmap, ctx["inputs"] if inputs is None else inputs, chain)
    except (DbProblem, h.HepError):
        out.abort()
        raise


def _run_one_database_body(ctx, db, label, out, matcher, tmap, inputs, chain):
    args, cfg = ctx["args"], ctx["cfg"]
    dbinfo, inventory = None, None
    with out:
        if not args.parse_only:
            dbinfo = inspect_db(db)
            if not args.skip_db_inventory:
                inventory = fasta_inventory(db)
                dbinfo["sequences"] = inventory["ids"]
                b.check_inventory(inventory, db, label, cfg, matcher)
        for species, fasta in inputs:
            h.log("%s [%s]: %s" % (species, label, "parsing existing results" if args.parse_only else "running jackhmmer"))
            all_ids = h.load_fasta_ids(fasta)
            if args.parse_only:
                raw = find_existing(args.parse_only, species)
                if raw is None:
                    h.warn("%s: no earlier result file found; all its genes are NOT_RUN" % species)
                    searched, rows = set(), []
                else:
                    complete, rows = parse_domtbl(raw)
                    searched = read_searched(raw)
                    if searched is None:
                        # a file made by hand: one finished search of the whole input, or nothing can be trusted
                        searched = set(all_ids) if complete else set()
                        if not complete:
                            h.warn("%s: %s does not end with '# [ok]', so the search may have been cut short; all its "
                                   "genes are NOT_RUN" % (species, raw))
                    rows = [r for r in rows if h.canonical_id(r["qid"]) in searched]
                failed_ids = set()
            else:
                raw, searched, failed_ids = run_search(species, fasta, db, cfg, args, out)
                _complete, rows = parse_domtbl(raw)
            hit_rows, counted, in_phylo = evaluate(rows, cfg, matcher, tmap, label)
            for r in hit_rows:
                r["species"] = species
            not_run = (set(all_ids) - searched) | failed_ids
            statuses = h.classify_genes(all_ids, counted=counted, in_phylogeny=in_phylo, not_run=not_run)
            tally = out.add_species(species, fasta, hit_rows, statuses)
            h.log("%s [%s]: %d genes, %d passed, %d with a counted hit, %d in-phylogeny only, %d NOT_RUN"
                  % (species, label, tally["input"], tally["passed"], tally[h.HIT], tally[h.IN_PHYLOGENY_ONLY],
                     tally[h.NOT_RUN]))
        if args.parse_only:
            tools = [{"name": "jackhmmer (HMMER), results parsed from an earlier run", "version": None}]
        else:
            tools = [{"name": "jackhmmer (HMMER)",
                      "version": h.tool_version(["jackhmmer", "-h"], r"HMMER \d+\.\d+(\.\d+)?")}]
        info = h.run_info_base(ctx)
        info.update(tools=tools, database=(os.path.abspath(db) if db else None), database_info=dbinfo,
                    database_inventory=inventory,
                    parsed_from=(os.path.abspath(args.parse_only) if args.parse_only else None),
                    chain=chain, analysis_species=[sp for sp, _ in ctx["inputs"]], match_level=cfg["match_level"],
                    ignore_in_phylogeny=bool(cfg["ignore_in_phylogeny"]),
                    ignore_taxa=h.as_list(cfg.get("ignore_taxa")))
        out.close(info, keep_work=args.keep_work, quiet=args.quiet)


def main(argv=None):
    ctx = h.start_module(MODULE, DESCRIPTION, argv, extra=add_args)
    args, cfg = ctx["args"], ctx["cfg"]
    try:
        check_settings(cfg)
        tmap = b.read_taxon_map(cfg["taxon_map"]) if cfg.get("taxon_map") else None
    except h.HepError as e:
        h.die(str(e))
    matcher = h.TaxonMatcher([sp for sp, _ in ctx["inputs"]], h.as_list(cfg.get("ignore_taxa")), cfg["match_level"])

    if args.parse_only:
        if args.db:
            h.warn("--db is ignored with --parse-only")
        runs = [(None, args.label or MODULE, args.output)]
    else:
        if not args.db:
            h.die("--db is required (or use --parse-only to read existing results)")
        if not h.find_tool("jackhmmer"):
            h.die("jackhmmer not found in PATH (install HMMER; see docs/setup.md)")
        multi = len(args.db) > 1
        runs, used = [], set()
        for db in args.db:
            label = b.db_label(db)
            if multi:
                label = "%s-%s" % (args.label, label) if args.label else label
                if label in used:
                    h.die("two databases would get the same label '%s'; give them different names" % label)
                used.add(label)
                runs.append((db, label, os.path.join(args.output, label)))
            else:
                runs.append((db, args.label or MODULE, args.output))
    sequential = bool(args.sequential)
    if sequential and args.parse_only:
        h.die("--sequential cannot be used with --parse-only")
    if sequential and len(runs) < 2:
        h.warn("--sequential needs at least two --db; it has no effect here")
        sequential = False
    failed, not_searched = [], []
    inputs, prev = ctx["inputs"], None
    for pos, (db, label, outdir) in enumerate(runs, 1):
        chain = {"sequential": True, "position": pos, "of": len(runs), "input_from": prev} if sequential else None
        try:
            run_one_database(ctx, db, label, outdir, matcher, tmap, inputs, chain)
        except (DbProblem, h.HepError) as e:
            h.warn("%s: %s" % (label, e))
            failed.append(label)
            if sequential:
                not_searched = [r[1] for r in runs[pos:]]
                break
            continue
        if sequential:
            inputs = [(sp, os.path.join(outdir, "passed", sp + ".faa")) for sp, _ in inputs]
            prev = label
    if sequential and not failed:
        h.log("sequential run finished; the genes that passed every database are in %s"
              % os.path.join(runs[-1][2], "passed"))
    if failed:
        extra = ""
        if not_searched:
            extra = ". The sequential chain stopped there; NOT searched: %s. Genes cannot pass a database that was never searched" % ", ".join(not_searched)
        h.die("%d of %d database run(s) did not complete: %s%s" % (len(failed), len(runs), ", ".join(failed), extra))


if __name__ == "__main__":
    main()
    sys.exit(h.exit_status())
