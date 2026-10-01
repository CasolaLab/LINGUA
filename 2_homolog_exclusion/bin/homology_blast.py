#!/usr/bin/env python3
"""
homology_blast.py - protein homology search with BLAST+ (blastp, or tblastn against nucleotide databases).

Question answered: does the candidate gene have a homolog in the user's database?

  - Against a protein database of other species (blastp): a young gene should not match genes from other species.
  - Against organelle genomes (tblastn, a nucleotide database): a hit means a leaked organellar gene posing as a
    lineage-specific gene.

ANY hit counts against a candidate (subject to the E-value and the optional coverage and identity settings). The
database must not contain the analysis species or close relatives: hits from the analysis species (the input file
names) and from species in ignore_taxa are reported but not counted (status IN_PHYLOGENY_ONLY). For an organelle
database, set ignore_in_phylogeny=no: a hit to the analysis species' own organelle genome is exactly the signal.

Several databases can be searched in one command (--db given more than once); each is its own labelled run.

Nothing here uses the network. Python 3.11+, standard library only. Needs BLAST+ (blastp or tblastn, blastdbcmd) on
PATH for a run, not for --parse-only.

  run:         python3 homology_blast.py -i proteins/ -o out/ --db /data/db/proteins --species-map species.tsv
  organelle:   python3 homology_blast.py -i proteins/ -o out/ --db /data/db/organelle --config presets/organelle_tblastn.conf --label organelle
  parse only:  python3 homology_blast.py -i proteins/ -o out/ --parse-only earlier_out/raw
"""

import hashlib
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402

MODULE = "homology-blast"
DESCRIPTION = __doc__.strip().split("\n\n")[0]

# The exact output format this module asks BLAST+ for. Parsing requires this column order.
OUTFMT_FIELDS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
NCOLS = len(OUTFMT_FIELDS.split())
PROGRAMS = ("blastp", "tblastn")
DB_SEQ_RE = re.compile(r"([\d,]+) sequences")
DB_TITLE_RE = re.compile(r"^Database:\s*(.*)$", re.M)


class DbProblem(Exception):
    """A problem with one database; the others are still searched."""


# --------------------------------------------------------------------------
# Reading BLAST+ output
# --------------------------------------------------------------------------
def parse_blast_output(path):
    """(confirmed, unconfirmed, rows) for a BLAST+ outfmt-7 file made with OUTFMT_FIELDS; rows are dicts."""
    confirmed, unconfirmed, raw_rows = h.read_blast_tabular(path, NCOLS, OUTFMT_FIELDS)
    rows = []
    for cols in raw_rows:
        try:
            rows.append({"qid": cols[0], "sid": cols[1], "pident": float(cols[2]),
                         "qstart": int(cols[6]), "qend": int(cols[7]), "sstart": int(cols[8]), "send": int(cols[9]),
                         "evalue": float(cols[10]), "bitscore": float(cols[11]),
                         "qlen": int(cols[12]), "slen": int(cols[13])})
        except ValueError as e:
            raise h.HepError("%s: cannot read a number in an alignment line for %s (%s)" % (path, cols[0], e))
    return confirmed, unconfirmed, rows


def read_taxon_map(path):
    """Two columns, accession then species, tab- or comma-separated. Lines starting with '#' and a header are skipped."""
    if not os.path.isfile(path):
        raise h.HepError("taxon map not found: %s" % path)
    mapping = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split(",")
            if len(parts) < 2 or not parts[0].strip() or not parts[1].strip():
                raise h.HepError("%s: expected 'accession<TAB>species', got: %r" % (path, line))
            mapping[parts[0].strip()] = parts[1].strip()
    if not mapping:
        raise h.HepError("%s: the taxon map is empty" % path)
    return mapping


def target_species(sid, cfg, tmap):
    """Species of a database hit: the '__Species_Name' suffix of its ID (species_source=header) or the taxon map."""
    if cfg.get("species_source") == "map":
        return (tmap or {}).get(sid) or (tmap or {}).get(h.canonical_id(sid))
    return h.split_id(sid)[1]


# --------------------------------------------------------------------------
# Deciding what counts
# --------------------------------------------------------------------------
def evaluate(rows, cfg, program, matcher, tmap, label):
    """
    Group alignment lines into one hit per (gene, database sequence), apply the rules, and return
    (hit_rows, counted_genes, in_phylogeny_genes).
      - the E-value is applied to each alignment line;
      - query coverage = percent of the query covered by the union of that hit's alignments; target coverage likewise
        for the target (blastp only: for tblastn the target coordinates are nucleotides, so it is not reported);
      - a hit counts if it passes the E-value, min_qcov, min_tcov and min_pident;
      - a hit that qualifies but comes from an analysis species or an ignore_taxa species does not count when
        ignore_in_phylogeny is on; its gene becomes IN_PHYLOGENY_ONLY unless another hit counts.
    """
    evalue_max = float(cfg["evalue"])
    min_q = float(cfg["min_qcov"] or 0)
    min_t = float(cfg["min_tcov"] or 0)
    min_id = cfg.get("min_pident")
    min_id = float(min_id) if min_id is not None else None
    exempt = bool(cfg["ignore_in_phylogeny"])

    groups = {}
    for r in rows:
        groups.setdefault((h.canonical_id(r["qid"]), r["sid"]), []).append(r)

    hit_rows, counted, in_phylo = [], set(), set()
    for (gene, sid), rs in sorted(groups.items()):
        sp = target_species(sid, cfg, tmap)
        passing = [r for r in rs if r["evalue"] <= evalue_max]
        best = min(rs, key=lambda r: (r["evalue"], -r["bitscore"]))
        row = {"gene_id": gene, "target": sid, "target_species": sp, "source": label, "counts_as_hit": False,
               "evalue": best["evalue"], "bitscore": best["bitscore"], "pident": best["pident"]}
        if not passing:
            row["note"] = "E-value above %g" % evalue_max
            hit_rows.append(row)
            continue
        best_p = min(passing, key=lambda r: (r["evalue"], -r["bitscore"]))
        row.update(evalue=best_p["evalue"], bitscore=best_p["bitscore"], pident=best_p["pident"])
        qcov = 100.0 * h.union_length([(r["qstart"], r["qend"]) for r in passing]) / max(1, passing[0]["qlen"])
        row["qcov"] = round(qcov, 1)
        tcov = None
        if program == "blastp":
            tcov = 100.0 * h.union_length([(r["sstart"], r["send"]) for r in passing]) / max(1, passing[0]["slen"])
            row["tcov"] = round(tcov, 1)
        if qcov < min_q:
            row["note"] = "query coverage %.1f%% below %g%%" % (qcov, min_q)
        elif tcov is not None and tcov < min_t:
            row["note"] = "target coverage %.1f%% below %g%%" % (tcov, min_t)
        elif min_id is not None and best_p["pident"] < min_id:
            row["note"] = "identity %.1f%% below %g%%" % (best_p["pident"], min_id)
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
# Databases
# --------------------------------------------------------------------------
def db_label(db):
    stem = os.path.basename(db.rstrip("/\\"))
    for ext in (".fasta", ".faa", ".fa", ".fna"):
        if stem.lower().endswith(ext):
            stem = stem[:-len(ext)]
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", stem) or "database"


def inspect_db(db, program):
    """
    Look at a formatted BLAST database with blastdbcmd (local, no network). Returns
    {'type', 'title', 'sequences'} and checks that its type suits the program. Raises DbProblem otherwise.
    """
    want = "prot" if program == "blastp" else "nucl"
    info = None
    for dbtype in ("prot", "nucl"):
        try:
            p = subprocess.run(["blastdbcmd", "-db", db, "-dbtype", dbtype, "-info"], stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, universal_newlines=True, timeout=120)
        except (OSError, subprocess.SubprocessError) as e:
            raise DbProblem("could not run blastdbcmd (is BLAST+ installed?): %s" % e)
        if p.returncode == 0 and DB_SEQ_RE.search(p.stdout):
            info = (dbtype, p.stdout)
            break
    if info is None:
        hint = " It looks like a FASTA file: format it first with makeblastdb -parse_seqids." if os.path.isfile(db) else ""
        raise DbProblem("%s is not a formatted BLAST database.%s" % (db, hint))
    dbtype, text = info
    if dbtype != want:
        raise DbProblem("%s is a %s database, but program=%s needs a %s database%s"
                        % (db, "protein" if dbtype == "prot" else "nucleotide", program,
                           "protein" if want == "prot" else "nucleotide",
                           " (organelle genomes are nucleotide: use program=tblastn)" if dbtype == "nucl" else ""))
    m, t = DB_SEQ_RE.search(text), DB_TITLE_RE.search(text)
    return {"type": dbtype, "title": t.group(1).strip() if t else None,
            "sequences": int(m.group(1).replace(",", "")) if m else None}


def species_inventory(db, dbtype):
    """
    Read every sequence ID of a formatted database (local, with blastdbcmd) and count how many carry each
    '__Species_Name' tag. Returns {'ids': N, 'tagged': M, 'species': {name: count}, 'first_ids': [...]}.
    """
    inv = {"ids": 0, "tagged": 0, "species": {}, "first_ids": []}
    try:
        p = subprocess.Popen(["blastdbcmd", "-db", db, "-dbtype", dbtype, "-entry", "all", "-outfmt", "%a"],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
    except OSError:
        return inv
    for line in p.stdout:
        sid = line.strip()
        if not sid:
            continue
        inv["ids"] += 1
        if len(inv["first_ids"]) < 5:
            inv["first_ids"].append(sid)
        sp = h.split_id(sid)[1]
        if sp:
            inv["tagged"] += 1
            inv["species"][sp] = inv["species"].get(sp, 0) + 1
    p.wait()
    return inv


def build_index(fasta, label, index_dir, program):
    """
    Format a FASTA file as a BLAST database (makeblastdb -parse_seqids, local) and return its path prefix. The index is a
    derived cache: it is rebuilt only if the FASTA changed (name, size or modification time), so a shared --index-dir
    is built once and reused. Raises DbProblem if makeblastdb fails.
    """
    dbtype = "prot" if program == "blastp" else "nucl"
    if not h.find_tool("makeblastdb"):
        raise DbProblem("makeblastdb not found in PATH (install BLAST+); it is needed to index %s" % fasta)
    os.makedirs(index_dir, exist_ok=True)
    key = hashlib.sha1(os.path.abspath(fasta).encode("utf-8")).hexdigest()[:8]
    stem = re.sub(r"[^A-Za-z0-9_.\-]", "_", os.path.splitext(os.path.basename(fasta))[0]) or "db"
    prefix = os.path.join(index_dir, "%s_%s" % (stem, key))
    st = os.stat(fasta)
    stamp = "fasta=%s\nsize=%d\nmtime_ns=%d\ndbtype=%s\n" % (os.path.abspath(fasta), st.st_size, st.st_mtime_ns, dbtype)
    stamp_path = prefix + ".hep_index"
    try:
        with open(stamp_path, "r", encoding="utf-8") as fh:
            if fh.read() == stamp and any(os.path.exists(prefix + ext) for ext in (".pin", ".nin")):
                h.log("%s: reusing the BLAST index built earlier from %s" % (label, fasta))
                return prefix
    except OSError:
        pass
    h.log("%s: building a BLAST index from %s (makeblastdb, local; it is a cache and can be deleted)" % (label, fasta))
    log = prefix + ".makeblastdb.log"
    rc = h.run_command(["makeblastdb", "-in", fasta, "-dbtype", dbtype, "-out", prefix, "-parse_seqids", "-title", label], log)
    if rc != 0:
        try:
            tail = " ".join(open(log, "r", encoding="utf-8", errors="replace").read().split())[-300:]
        except OSError:
            tail = ""
        raise DbProblem("makeblastdb failed on %s (exit %d): %s" % (fasta, rc, tail))
    with open(stamp_path, "w", encoding="utf-8") as fh:
        fh.write(stamp)
    size = sum(os.path.getsize(prefix + ext) for ext in (".pin", ".psq", ".phr", ".nin", ".nsq", ".nhr") if os.path.exists(prefix + ext))
    h.log("%s: BLAST index built (%.0f MB) in %s" % (label, size / 1e6, index_dir))
    return prefix


# --------------------------------------------------------------------------
# Running BLAST+
# --------------------------------------------------------------------------
def blast_command(program, chunk, db, out_tsv, cfg, threads):
    cmd = [program, "-query", chunk, "-db", db, "-evalue", "%g" % float(cfg["evalue"]),
           "-outfmt", "7 " + OUTFMT_FIELDS, "-max_target_seqs", str(int(cfg["max_target_seqs"])),
           "-num_threads", str(int(threads)), "-out", out_tsv, "-seg", "yes" if cfg.get("seg") else "no"]
    if cfg.get("max_hsps") is not None:
        cmd += ["-max_hsps", str(int(cfg["max_hsps"]))]
    if cfg.get("ungapped"):
        cmd += ["-ungapped"]
    if cfg.get("comp_based_stats") is not None:
        cmd += ["-comp_based_stats", str(int(cfg["comp_based_stats"]))]
    if cfg.get("dbsize") is not None:
        cmd += ["-dbsize", str(int(cfg["dbsize"]))]
    return cmd


def check_settings(cfg):
    program = cfg.get("program")
    if program not in PROGRAMS:
        raise h.HepError("program must be blastp or tblastn, not %r" % program)
    if cfg.get("ungapped") and cfg.get("comp_based_stats") not in (0, "0"):
        raise h.HepError("ungapped=yes needs comp_based_stats=0: BLAST+ refuses composition-adjusted ungapped searches "
                         "(set --comp-based-stats 0)")
    if program == "tblastn" and float(cfg.get("min_tcov") or 0) > 0:
        raise h.HepError("min_tcov cannot be used with tblastn: the target coordinates are nucleotides, so target "
                         "coverage is not defined")
    if cfg.get("species_source") not in ("header", "map"):
        raise h.HepError("species_source must be 'header' or 'map', not %r" % cfg.get("species_source"))
    if cfg.get("species_source") == "map" and not cfg.get("taxon_map"):
        raise h.HepError("species_source=map needs taxon_map (a file with accession and species)")
    if cfg.get("match_level") not in h.MATCH_LEVELS:
        raise h.HepError("match_level must be one of %s, not %r" % (", ".join(h.MATCH_LEVELS), cfg.get("match_level")))


def run_search(species, fasta, db, cfg, args, out):
    """
    Run BLAST+ on the chunks of one species against one database. Returns (raw_path, not_run_ids): the merged raw
    result file, and the IDs of genes in chunks that failed or did not finish (NOT_RUN, never NO_HIT).
    """
    program = cfg["program"]
    base = os.path.join(out.workdir, species)
    chunk_dir, res_dir, log_dir = (os.path.join(base, d) for d in ("chunks", "results", "logs"))
    for d in (res_dir, log_dir):
        os.makedirs(d, exist_ok=True)
    chunks = h.chunk_fasta(fasta, chunk_dir, int(cfg["chunk_size"]))
    by_name = {os.path.splitext(os.path.basename(c))[0]: c for c in chunks}

    def job(name):
        tsv = os.path.join(res_dir, name + ".tsv")
        rc = h.run_command(blast_command(program, by_name[name], db, tsv, cfg, args.threads),
                           os.path.join(log_dir, name + ".log"))
        h.check_rc(rc, "%s on %s/%s" % (program, species, name))
        confirmed, unconfirmed, _rows = parse_blast_output(tsv)
        expected = {rid for rid, _h, _s in h.read_fasta(by_name[name])}
        if unconfirmed or not expected <= confirmed:
            raise h.CommandError("%s output for %s/%s is incomplete (no completion footer)" % (program, species, name))

    res = h.run_jobs(sorted(by_name), job, workers=int(cfg["parallel_jobs"]),
                     marker_dir=os.path.join(base, "markers"), resume=args.resume)
    not_run = set()
    for name, msg in res["failed"].items():
        h.warn("%s/%s failed: %s" % (species, name, msg))
        not_run |= {h.canonical_id(rid) for rid, _h, _s in h.read_fasta(by_name[name])}
    raw_dir = os.path.join(out.outdir, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    raw_path = os.path.join(raw_dir, species + ".blast.tsv")
    with open(raw_path, "w", encoding="utf-8") as merged:
        for name in sorted(by_name):
            tsv = os.path.join(res_dir, name + ".tsv")
            if name not in res["failed"] and os.path.exists(tsv):
                with open(tsv, "r", encoding="utf-8") as fh:
                    merged.write(fh.read())
    return raw_path, not_run


def find_existing(parse_only, species):
    if os.path.isfile(parse_only):
        return parse_only
    p = os.path.join(parse_only, species + ".blast.tsv")
    return p if os.path.isfile(p) else None


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def add_args(parser):
    g = parser.add_argument_group("homology-blast options")
    g.add_argument("--db", action="append", default=[], metavar="FASTA_OR_PREFIX",
                   help="the database: a species-tagged FASTA file (a BLAST index is built from it and cached), or the "
                        "path prefix of a database already formatted with makeblastdb -parse_seqids. Give it more than "
                        "once to search several databases, each as its own labelled run. Required unless --parse-only")
    g.add_argument("--index-dir", default=None, metavar="DIR",
                   help="where the BLAST index of a FASTA database is kept (default: db_index inside each run's output "
                        "folder). Point several runs at one folder to build the index once and reuse it")
    g.add_argument("--skip-db-inventory", action="store_true",
                   help="do not read all the database's IDs to count species and look for your analysis species in it "
                        "(saves time on very large databases)")
    g.add_argument("--parse-only", default=None, metavar="PATH",
                   help="do not run BLAST+: read existing results, a folder holding <species>.blast.tsv (or one file "
                        "for a single species), made with outfmt '7 " + OUTFMT_FIELDS + "'")
    g.add_argument("--sequential", action="store_true",
                   help="with several --db: search them in the order given, each one only against the genes that passed "
                        "the previous one (iterative filtering). Without it every database searches every gene. If one "
                        "database fails, the chain stops, because genes cannot pass a database that was never searched")


def check_inventory(inv, db, label, cfg, matcher):
    """
    What the database's own IDs say about it. Raises DbProblem if it cannot be used; warns about the mistakes that make
    results wrong without making anything fail: no species tags, and the analysis species inside the database.
    """
    if not inv["ids"]:
        raise DbProblem("%s contains no sequences" % db)
    if any(i.startswith("gnl|BL_ORD_ID") for i in inv["first_ids"]):
        raise DbProblem("%s was built without -parse_seqids, so BLAST+ reports its sequences as gnl|BL_ORD_ID|n instead of "
                        "their IDs. Rebuild it with: makeblastdb -parse_seqids (or give --db the FASTA file instead)" % db)
    if cfg["species_source"] == "header":
        if inv["tagged"] == 0:
            if cfg["ignore_in_phylogeny"]:
                h.warn("%s: none of the %d sequence IDs in this database carries a '__Species_Name' tag, so no hit can be "
                       "recognised as coming from an analysis species and all such hits will be COUNTED. Tag the IDs "
                       "(make-db, or Stage 0 --add-species), or use species_source=map, or set ignore_in_phylogeny=no."
                       % (label, inv["ids"]))
        elif inv["tagged"] < inv["ids"]:
            h.warn("%s: %d of %d sequence IDs carry no species tag; hits to them are counted as hits of unknown species."
                   % (label, inv["ids"] - inv["tagged"], inv["ids"]))
    present = {sp: n for sp, n in inv["species"].items() if matcher.match(sp)}
    if present:
        shown = ", ".join("%s (%d sequences)" % (sp, n) for sp, n in sorted(present.items())[:8])
        h.warn("%s: species you are analysing, or excluded with ignore_taxa, are INSIDE this database: %s. %s The safest fix "
               "is to remove them from the database."
               % (label, shown, "Their hits are reported but not counted (IN_PHYLOGENY_ONLY)." if cfg["ignore_in_phylogeny"]
                  else "Their hits WILL be counted (ignore_in_phylogeny=no)."))


def run_one_database(ctx, db, label, outdir, matcher, tmap, inputs=None, chain=None):
    """
    One database (or one parse-only read) as a complete labelled run. Raises DbProblem on a database problem, after
    removing the output folder this run created, so no header-only tables are left behind.
    `inputs` is the list of (species, FASTA) to search (default: the command's inputs; in a sequential chain, the
    survivors of the previous database); `chain` describes the position in a sequential chain, for run.json.
    """
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
    program = cfg["program"]
    dbinfo, inventory, prefix = None, None, db
    with out:
        if not args.parse_only:
            if os.path.isfile(db):      # a FASTA file: build (or reuse) a BLAST index for it
                prefix = build_index(db, label, args.index_dir or os.path.join(out.outdir, "db_index"), program)
            dbinfo = inspect_db(prefix, program)
            if not args.skip_db_inventory:
                inventory = species_inventory(prefix, dbinfo["type"])
                check_inventory(inventory, prefix, label, cfg, matcher)
        for species, fasta in inputs:
            h.log("%s [%s]: %s" % (species, label, "parsing existing results" if args.parse_only else "running " + program))
            all_ids = h.load_fasta_ids(fasta)
            if args.parse_only:
                raw = find_existing(args.parse_only, species)
                if raw is None:
                    h.warn("%s: no earlier result file found; all its genes are NOT_RUN" % species)
                    confirmed, rows, failed_ids = set(), [], set()
                else:
                    confirmed, unconfirmed, rows = parse_blast_output(raw)
                    failed_ids = set()
                    if unconfirmed:
                        h.warn("%s: %d gene(s) are in an unfinished part of the results and are NOT_RUN"
                               % (species, len(unconfirmed)))
                    rows = [r for r in rows if r["qid"] in confirmed]
            else:
                raw, failed_ids = run_search(species, fasta, prefix, cfg, args, out)
                confirmed, unconfirmed, rows = parse_blast_output(raw)
            seen = {h.canonical_id(q) for q in confirmed}
            hit_rows, counted, in_phylo = evaluate(rows, cfg, program, matcher, tmap, label)
            for r in hit_rows:
                r["species"] = species
            not_run = (set(all_ids) - seen) | failed_ids
            statuses = h.classify_genes(all_ids, counted=counted, in_phylogeny=in_phylo, not_run=not_run)
            tally = out.add_species(species, fasta, hit_rows, statuses)
            h.log("%s [%s]: %d genes, %d passed, %d with a counted hit, %d in-phylogeny only, %d NOT_RUN"
                  % (species, label, tally["input"], tally["passed"], tally[h.HIT], tally[h.IN_PHYLOGENY_ONLY],
                     tally[h.NOT_RUN]))
        tools = []
        if args.parse_only:
            tools.append({"name": "%s (BLAST+), results parsed from an earlier run" % program.upper(), "version": None})
        else:
            tools.append({"name": "%s (BLAST+)" % program.upper(),
                          "version": h.tool_version([program, "-version"], r"\d+\.\d+\.\d+\+?")})
        info = h.run_info_base(ctx)
        info.update(tools=tools, program=program, database=(os.path.abspath(db) if db else None),
                    database_index=(os.path.abspath(prefix) if db and prefix != db else None), database_info=dbinfo,
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
        tmap = read_taxon_map(cfg["taxon_map"]) if cfg.get("taxon_map") else None
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
        for tool in (cfg["program"], "blastdbcmd"):
            if not h.find_tool(tool):
                h.die("%s not found in PATH (install BLAST+; see docs/setup.md)" % tool)
        multi = len(args.db) > 1
        runs, used = [], set()
        for db in args.db:
            label = db_label(db)
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
            # the next database searches only what passed this one
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
