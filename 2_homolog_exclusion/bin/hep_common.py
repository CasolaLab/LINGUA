#!/usr/bin/env python3
"""
hep_common.py - shared helpers for homolog-exclusion-pipeline (LINGUA Stage 2).

Python 3.11+, standard library only. Imported by the module scripts in this folder
(domain_interpro.py, domain_cdd.py, homology_blastp.py, homology_jackhmmer.py, combine.py).

Contents
--------
  1. logging helpers
  2. configuration: defaults -> config file -> command-line flags, with a record of what changed
  3. FASTA helpers and ID handling
  4. species handling (from input file names, optional species map)
  5. output schema: hits.tsv, gene_status.tsv, passed/<species>.faa, passed_ids.tsv, run.json
  6. runner: external commands, parallel chunk jobs, resume markers
  7. tool discovery and versions
  8. command-line scaffolding shared by every module
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

__version__ = "0.1.0"
PIPELINE_NAME = "homolog-exclusion-pipeline"

# Exit status: 0 = finished and every gene was searched; 1 = the run could not be done (an error, nothing usable);
# 2 = finished and the results are written, but some genes are NOT_RUN (or, for combine, INCOMPLETE): they do not pass.
EXIT_OK, EXIT_ERROR, EXIT_INCOMPLETE = 0, 1, 2
_STATE = {"not_run": 0}


def note_incomplete(n):
    """Record that n genes were not properly searched, so the command exits with EXIT_INCOMPLETE."""
    _STATE["not_run"] += int(n)


def exit_status():
    return EXIT_INCOMPLETE if _STATE["not_run"] else EXIT_OK

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULTS_DIR = os.path.normpath(os.path.join(HERE, os.pardir, "defaults"))
FASTA_EXTS = (".faa", ".fa", ".fasta")

# Gene status values (one per input gene per run).
HIT = "HIT"                              # counted evidence of a domain or homolog
NO_HIT = "NO_HIT"                        # ran, nothing counted
EXCLUDED_ONLY = "EXCLUDED_ONLY"          # only non-domain signals (e.g. MobiDB-Lite, Coils, PROSITE patterns)
SPURIOUS = "SPURIOUS"                    # flagged as a spurious gene model (AntiFam)
IN_PHYLOGENY_ONLY = "IN_PHYLOGENY_ONLY"  # hits only from analysis species or ignored taxa; reported, not counted
NOT_RUN = "NOT_RUN"                      # job failed or gene missing from output; never treated as NO_HIT
STATUSES = (HIT, NO_HIT, EXCLUDED_ONLY, SPURIOUS, IN_PHYLOGENY_ONLY, NOT_RUN)


# --------------------------------------------------------------------------
# 1. logging
# --------------------------------------------------------------------------
def log(msg):
    sys.stderr.write("[%s] %s\n" % (datetime.now().strftime("%H:%M:%S"), msg))
    sys.stderr.flush()


def warn(msg):
    log("WARNING: " + msg)


def die(msg, code=1):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(code)


class HepError(Exception):
    """Raised for user-facing problems (bad config, bad input)."""


# --------------------------------------------------------------------------
# 2. configuration
# --------------------------------------------------------------------------
def coerce(value):
    """Turn a config string into None, bool, int, float, or str."""
    if value is None:
        return None
    if not isinstance(value, str):
        return value
    v = value.strip()
    if v == "":
        return None
    low = v.lower()
    if low in ("true", "yes"):
        return True
    if low in ("false", "no"):
        return False
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(v)
    except ValueError:
        return v


def as_list(value):
    """Comma-separated string (or None) -> list of non-empty stripped items."""
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(x).strip() for x in value if str(x).strip()]
    return [x.strip() for x in str(value).split(",") if x.strip()]


def parse_conf(path):
    """
    Parse a config file. Returns (global_keys, sections).
    Format: key=value lines, '#' starts a comment, optional [module] section headers.
    Keys before any section are global (apply to every module that defines the key).
    """
    global_keys = {}
    sections = {}
    current = None
    with open(path, "r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                current = line[1:-1].strip()
                sections.setdefault(current, {})
                continue
            if "=" not in line:
                raise HepError("%s line %d: expected key=value, got: %s" % (path, lineno, raw.rstrip()))
            key, _, val = line.partition("=")
            key = key.strip()
            if not key:
                raise HepError("%s line %d: empty key" % (path, lineno))
            target = sections[current] if current is not None else global_keys
            target[key] = val.strip()
    return global_keys, sections


def load_defaults(module, defaults_dir=None):
    """Load defaults/common.conf then defaults/<module>.conf. Returns {key: coerced value}."""
    defaults_dir = defaults_dir or DEFAULTS_DIR
    cfg = {}
    for name in ("common", module):
        path = os.path.join(defaults_dir, name + ".conf")
        if os.path.isfile(path):
            g, _ = parse_conf(path)
            for k, v in g.items():
                cfg[k] = coerce(v)
        elif name == module:
            raise HepError("no defaults file for module '%s': %s" % (module, path))
    return cfg


def resolve_config(module, config_path=None, overrides=None, defaults_dir=None):
    """
    Resolve settings for a module. Later sources win:
      defaults/common.conf, defaults/<module>.conf, --config file (global keys, then [module] section),
      overrides (command-line flags and --set).
    Returns (cfg, non_default) where non_default = {key: {"default": d, "value": v}}.
    Unknown keys in a [module] section or in overrides raise HepError (catches typos);
    unknown global keys in a config file are ignored, since they may belong to another module.
    """
    defaults = load_defaults(module, defaults_dir)
    cfg = dict(defaults)

    if config_path:
        if not os.path.isfile(config_path):
            raise HepError("config file not found: %s" % config_path)
        g, sections = parse_conf(config_path)
        for k, v in g.items():
            if k in defaults:
                cfg[k] = coerce(v)
        for k, v in sections.get(module, {}).items():
            if k not in defaults:
                raise HepError("%s: unknown setting '%s' in [%s]. Known: %s"
                               % (config_path, k, module, ", ".join(sorted(defaults))))
            cfg[k] = coerce(v)

    for k, v in (overrides or {}).items():
        if k not in defaults:
            raise HepError("unknown setting '%s' for %s. Known: %s"
                           % (k, module, ", ".join(sorted(defaults))))
        cfg[k] = coerce(v)

    non_default = {}
    for k, v in cfg.items():
        if v != defaults.get(k):
            non_default[k] = {"default": defaults.get(k), "value": v}
    return cfg, non_default


# --------------------------------------------------------------------------
# 3. FASTA helpers and ID handling
# --------------------------------------------------------------------------
def read_fasta(path):
    """Yield (record_id, header, sequence). record_id = first whitespace token; header excludes '>'."""
    header = None
    rid = None
    chunks = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if header is not None:
                    yield rid, header, "".join(chunks)
                header = line[1:].strip()
                rid = header.split()[0] if header else ""
                chunks = []
            else:
                if header is None:
                    raise HepError("%s: sequence data before first '>' header" % path)
                chunks.append(line.strip())
    if header is not None:
        yield rid, header, "".join(chunks)


def write_fasta_record(fh, header, seq, width=60):
    fh.write(">" + header + "\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")
    if not seq:
        fh.write("\n")


def split_id(token):
    """
    Split a FASTA ID token into (gene_id, species_suffix_or_None).
    Stage 0 --add-species appends '__<Species_Name>'; the suffix is split at the LAST '__'.
    """
    if "__" in token:
        gene, _, sp = token.rpartition("__")
        if gene and sp:
            return gene, sp
    return token, None


# Text after '__' that looks like a species name: a capital letter, then lower-case letters, then more name parts
# ('Homo_sapiens', 'Ceratodon_purpureus_R40', 'O_sativa'). Random characters ('FwSs-9Ag') do not match.
SPECIES_TAG_RE = re.compile(r"^[A-Z][a-z\-\.]*(_[A-Za-z0-9\.\-]+)*$")


def looks_like_species(tag):
    """True if `tag` (the text after '__' in an ID) looks like a species name and not like random characters."""
    return bool(tag) and bool(SPECIES_TAG_RE.match(tag))


def canonical_id(token):
    """Gene ID with any '__<Species>' suffix removed. Used on both sides of every ID match."""
    return split_id(token)[0]


def load_fasta_ids(path):
    """
    Return the list of canonical gene IDs in a FASTA file, in order.
    Raises HepError on duplicate or empty IDs (matching by ID would be ambiguous).
    """
    ids = []
    seen = set()
    for rid, _header, _seq in read_fasta(path):
        gid = canonical_id(rid)
        if not gid:
            raise HepError("%s: record with an empty ID" % path)
        if gid in seen:
            raise HepError("%s: duplicate gene ID '%s'" % (path, gid))
        seen.add(gid)
        ids.append(gid)
    return ids


def chunk_fasta(path, outdir, size, prefix="chunk", sort_by_length=True):
    """
    Split a FASTA file into files of at most `size` records. Returns the list of chunk paths.
    By default the records are sorted longest first, so the slowest chunks are numbered first and start first
    (jobs are handed out in chunk order): the long queries overlap with the many short ones instead of being left
    for the end, when most workers would sit idle. Ties keep the input order.
    """
    if size < 1:
        raise HepError("chunk size must be >= 1")
    os.makedirs(outdir, exist_ok=True)
    paths = []
    fh = None
    n = 0
    records = read_fasta(path)
    if sort_by_length:
        records = sorted(records, key=lambda r: -len(r[2]))   # sorted() is stable
    try:
        for _rid, header, seq in records:
            if n % size == 0:
                if fh:
                    fh.close()
                p = os.path.join(outdir, "%s_%05d.faa" % (prefix, len(paths)))
                paths.append(p)
                fh = open(p, "w", encoding="utf-8")
            write_fasta_record(fh, header, seq)
            n += 1
    finally:
        if fh:
            fh.close()
    return paths


# --------------------------------------------------------------------------
# 4. species handling
# --------------------------------------------------------------------------
def fasta_stem(path):
    base = os.path.basename(path)
    for ext in FASTA_EXTS:
        if base.lower().endswith(ext):
            return base[:-len(ext)]
    return base


def species_key(name):
    """Normalised species name for comparisons: lower case, spaces/hyphens -> underscores."""
    return re.sub(r"[\s\-]+", "_", str(name).strip()).lower()


MATCH_LEVELS = ("species", "exact", "genus")


def taxon_key(name, level="species"):
    """
    Comparable form of a species name at a chosen level:
      species - the first two words (Genus species): subspecies, varieties and strains count as the same species,
                so 'Marchantia polymorpha subsp. ruderalis' and 'Marchantia_polymorpha' give the same key;
      exact   - the whole normalised name, subspecies or strain included;
      genus   - the first word.
    """
    if level not in MATCH_LEVELS:
        raise HepError("match_level must be one of %s, not %r" % (", ".join(MATCH_LEVELS), level))
    parts = [p for p in species_key(name).split("_") if p]
    if level == "exact":
        return "_".join(parts)
    if level == "genus":
        return parts[0] if parts else ""
    return "_".join(parts[:2])


class TaxonMatcher(object):
    """
    Decides whether a database hit comes from a species that must not count (an analysis species, or one named in
    ignore_taxa). `analysis_species` and `ignore_taxa` are lists of names; `level` is a match level (see taxon_key).
    An ignore_taxa entry of a single word is a genus and matches every species of it, at any level.
    Returns the matching name, or None. A hit whose species is unknown (None or empty) never matches, so it is counted.
    """

    def __init__(self, analysis_species, ignore_taxa=(), level="species"):
        self.level = level
        self.entries = [(taxon_key(s, level), s) for s in analysis_species if str(s).strip()]
        self.genera = []
        for t in ignore_taxa:
            if not str(t).strip():
                continue
            if len([p for p in species_key(t).split("_") if p]) == 1:
                self.genera.append((taxon_key(t, "genus"), t))
            else:
                self.entries.append((taxon_key(t, level), t))

    def match(self, species):
        if not species or not str(species).strip():
            return None
        key = taxon_key(species, self.level)
        for k, name in self.entries:
            if k and k == key:
                return name
        genus = taxon_key(species, "genus")
        for g, name in self.genera:
            if g and g == genus:
                return name
        return None


def load_species_map(path):
    """
    Read a species map like Stage 0's: two columns (species name, file basename), header row,
    tab- or comma-separated. Returns {file_basename: species_name}.
    """
    if not os.path.isfile(path):
        raise HepError("species map not found: %s" % path)
    with open(path, "r", encoding="utf-8", newline="") as fh:
        text = fh.read()
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise HepError("%s: species map needs a header row and at least one entry" % path)
    delim = "\t" if "\t" in lines[0] else ","
    mapping = {}
    for row in list(csv.reader(lines, delimiter=delim))[1:]:
        if len(row) < 2 or not row[0].strip() or not row[1].strip():
            raise HepError("%s: bad row in species map: %r" % (path, row))
        base = row[1].strip()
        for ext in FASTA_EXTS:
            if base.lower().endswith(ext):
                base = base[:-len(ext)]
        mapping[base] = row[0].strip()
    return mapping


def species_from_path(path, strip_suffixes=("_final",), species_map=None):
    """Species for an input FASTA: species-map entry if present, else the file stem minus known suffixes."""
    stem = fasta_stem(path)
    if species_map and stem in species_map:
        return species_map[stem]
    for suf in strip_suffixes:
        if suf and stem.endswith(suf) and len(stem) > len(suf):
            stem = stem[:-len(suf)]
            break
    return stem


def discover_inputs(path, strip_suffixes=("_final",), species_map=None):
    """
    Input may be one FASTA file or a directory of per-species FASTA files (not recursive).
    Returns a sorted list of (species, fasta_path). Two files mapping to the same species is an error.
    """
    if os.path.isfile(path):
        files = [path]
    elif os.path.isdir(path):
        files = sorted(os.path.join(path, f) for f in os.listdir(path)
                       if f.lower().endswith(FASTA_EXTS) and os.path.isfile(os.path.join(path, f)))
    else:
        raise HepError("input not found: %s" % path)
    if not files:
        raise HepError("no FASTA files (%s) found in %s" % ("/".join(FASTA_EXTS), path))
    out = []
    seen = {}
    for f in files:
        sp = species_from_path(f, strip_suffixes, species_map)
        if sp in seen:
            raise HepError("two input files map to species '%s': %s and %s" % (sp, seen[sp], f))
        seen[sp] = f
        out.append((sp, f))
    return sorted(out)


# --------------------------------------------------------------------------
# 5. output schema
# --------------------------------------------------------------------------
HITS_COLUMNS = ["species", "gene_id", "target", "target_species", "source", "evalue", "bitscore",
                "qcov", "tcov", "pident", "counts_as_hit", "note"]
STATUS_COLUMNS = ["species", "gene_id", "status", "n_counted_hits"]


def classify_genes(all_ids, counted=(), excluded=(), spurious=(), in_phylogeny=(), not_run=()):
    """
    Map every gene ID to one status. Precedence (highest first):
    NOT_RUN > HIT > SPURIOUS > EXCLUDED_ONLY > IN_PHYLOGENY_ONLY > NO_HIT.
    Arguments are collections of gene IDs that fall in each category.
    """
    counted, excluded, spurious = set(counted), set(excluded), set(spurious)
    in_phylogeny, not_run = set(in_phylogeny), set(not_run)
    status = {}
    for gid in all_ids:
        if gid in not_run:
            status[gid] = NOT_RUN
        elif gid in counted:
            status[gid] = HIT
        elif gid in spurious:
            status[gid] = SPURIOUS
        elif gid in excluded:
            status[gid] = EXCLUDED_ONLY
        elif gid in in_phylogeny:
            status[gid] = IN_PHYLOGENY_ONLY
        else:
            status[gid] = NO_HIT
    return status


def passes(status, antifam_action="flag"):
    """Does a gene with this status pass the module (stay in passed/<species>.faa)?"""
    if status in (NO_HIT, EXCLUDED_ONLY, IN_PHYLOGENY_ONLY):
        return True
    if status == SPURIOUS:
        return antifam_action != "remove"
    return False  # HIT and NOT_RUN never pass


def _fmt(value):
    if value is None or value == "":
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


RUN_MARKER = ".hep_run"   # written into every output folder this tool creates


def check_output_folder(outdir, module, label, resume=False, force=False):
    """
    Refuse to write into a folder that already holds something, unless it is safe.
      - a new or empty folder is fine;
      - a folder with files but no marker is not ours: always refused (never mixed with other data);
      - a folder of another module or label is refused (two runs would overwrite each other's files);
      - a folder of this same run is refused unless resume=True (continue an interrupted run) or force=True
        (overwrite its results).
    """
    if not os.path.isdir(outdir) or not os.listdir(outdir):
        return
    marker = os.path.join(outdir, RUN_MARKER)
    if not os.path.isfile(marker):
        raise HepError("output folder %s already contains files that this tool did not create. "
                       "Choose a new or empty folder." % outdir)
    prev = {}
    with open(marker, "r", encoding="utf-8") as fh:
        for line in fh:
            k, _, v = line.strip().partition("=")
            prev[k] = v
    if prev.get("module") != module or prev.get("label") != label:
        raise HepError("output folder %s holds the results of '%s' (label '%s'), not '%s' (label '%s'). "
                       "Choose a different folder or give this run its own --label."
                       % (outdir, prev.get("module"), prev.get("label"), module, label))
    if not (resume or force):
        raise HepError("output folder %s already holds results of this run. Use --resume to continue an "
                       "interrupted run, --force to overwrite the results, or choose a new folder." % outdir)


def summary_table(counts):
    """Rows for summary.tsv: one per species plus a total, from the per-species tallies of add_species()."""
    cols = ["species", "input", "passed"] + list(STATUSES)
    rows = [cols]
    total = {c: 0 for c in cols[1:]}
    for sp in sorted(counts):
        rows.append([sp] + [counts[sp].get(c, 0) for c in cols[1:]])
        for c in cols[1:]:
            total[c] += counts[sp].get(c, 0)
    rows.append(["ALL"] + [total[c] for c in cols[1:]])
    return rows


def citation_note(module, outdir, tools, counts):
    """
    End-of-run note: how many genes passed, which tools this run used, and where the references are.
    The list of tools comes from what the module recorded (run.json), so only tools really used are named.
    Full references live in docs/citing.md; this note never repeats them.
    """
    total = sum(c.get("input", 0) for c in counts.values())
    passed = sum(c.get("passed", 0) for c in counts.values())
    used = ["%s %s" % (t.get("name", "unknown tool"), t.get("version") or "(version unknown)")
            for t in (tools or [])]
    used.append("%s %s" % (PIPELINE_NAME, __version__))
    not_run = sum(c.get(NOT_RUN, 0) for c in counts.values())
    warning = ("WARNING: %d gene(s) are NOT_RUN: their search failed or did not finish, so they do NOT pass and are not "
               "\"no hit\". See gene_status.tsv.\n" % not_run) if not_run else ""
    return (
        "Done ({module}): {passed} of {total} genes passed. Results are in {outdir}\n"
        "{warning}"
        "This run used: {used}.\n"
        "Please cite these tools and the {pipeline}. References: docs/citing.md in the pipeline repository.\n"
    ).format(module=module, passed=passed, total=total, outdir=outdir, used="; ".join(used),
             pipeline=PIPELINE_NAME, warning=warning)


class ModuleOutput(object):
    """
    Writes the files every module produces, inside `outdir`:
      hits.tsv, gene_status.tsv, passed/<species>.faa, passed_ids.tsv, run.json, work/ (scratch).
    Use as a context manager; call add_species() once per input file, then close().
    """

    def __init__(self, outdir, module, label=None, resume=False, force=False):
        self.outdir = outdir
        self.module = module
        self.label = label or module
        check_output_folder(outdir, module, self.label, resume, force)
        # A folder this call creates (or that was empty) is ours to remove if the run is aborted before any result exists.
        self._fresh = not os.path.isdir(outdir) or not os.listdir(outdir)
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, RUN_MARKER), "w", encoding="utf-8") as fh:
            fh.write("module=%s\nlabel=%s\n" % (module, self.label))
        self.workdir = os.path.join(outdir, "work")
        self.passed_dir = os.path.join(outdir, "passed")
        os.makedirs(self.passed_dir, exist_ok=True)
        os.makedirs(self.workdir, exist_ok=True)
        self._files = []
        self._hits = self._open("hits.tsv", HITS_COLUMNS)
        self._status = self._open("gene_status.tsv", STATUS_COLUMNS)
        self._passed_ids = self._open("passed_ids.tsv", ["species", "gene_id"])
        self.counts = {}

    def _open(self, name, columns):
        fh = open(os.path.join(self.outdir, name), "w", encoding="utf-8", newline="")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        self._files.append(fh)
        return w

    def add_species(self, species, fasta_path, hit_rows, statuses, antifam_action="flag"):
        """
        Record one species.
          hit_rows: list of dicts with keys from HITS_COLUMNS (species is filled in if missing;
                    counts_as_hit should be True/False or 'yes'/'no').
          statuses: {gene_id: status} keyed by canonical gene ID. Any gene in the FASTA that is
                    missing here is recorded as NOT_RUN (never NO_HIT).
        """
        statuses = dict(statuses)
        records = [(canonical_id(rid), header, seq) for rid, header, seq in read_fasta(fasta_path)]
        fasta_ids = set(g for g, _h, _s in records)
        if len(fasta_ids) != len(records):
            raise HepError("%s: duplicate gene IDs" % fasta_path)
        missing = [g for g in fasta_ids if g not in statuses]
        if missing:
            warn("%s: %d gene(s) have no result and are marked NOT_RUN (e.g. %s)"
                 % (species, len(missing), ", ".join(sorted(missing)[:3])))
            for g in missing:
                statuses[g] = NOT_RUN
        extra = [g for g in statuses if g not in fasta_ids]
        if extra:
            warn("%s: %d result ID(s) not in the input FASTA were ignored (e.g. %s)"
                 % (species, len(extra), ", ".join(sorted(extra)[:3])))
        for g, s in statuses.items():
            if s not in STATUSES:
                raise HepError("invalid status '%s' for gene %s" % (s, g))

        n_counted = {}
        for row in hit_rows:
            r = dict(row)
            r.setdefault("species", species)
            counted = r.get("counts_as_hit")
            r["counts_as_hit"] = "yes" if counted in (True, "yes") else "no"
            if r["counts_as_hit"] == "yes":
                n_counted[r["gene_id"]] = n_counted.get(r["gene_id"], 0) + 1
            self._hits.writerow([_fmt(r.get(c)) for c in HITS_COLUMNS])

        tally = {s: 0 for s in STATUSES}
        passed_path = os.path.join(self.passed_dir, species + ".faa")
        with open(passed_path, "w", encoding="utf-8") as pfh:
            for gid, header, seq in records:
                st = statuses[gid]
                tally[st] += 1
                self._status.writerow([species, gid, st, n_counted.get(gid, 0)])
                if passes(st, antifam_action):
                    write_fasta_record(pfh, header, seq)
                    self._passed_ids.writerow([species, gid])
        tally["input"] = len(records)
        tally["passed"] = sum(tally[s] for s in (NO_HIT, EXCLUDED_ONLY, IN_PHYLOGENY_ONLY)) + \
            (tally[SPURIOUS] if antifam_action != "remove" else 0)
        self.counts[species] = tally
        return tally

    def close(self, run_info=None, keep_work=False, quiet=False):
        """
        Finish the run: write run.json and citations.txt, and (unless quiet) print the end-of-run note.
        run_info["tools"] is a list of {"name": ..., "version": ...} for the external tools this run used.
        """
        for fh in self._files:
            fh.close()
        self._files = []
        note_incomplete(sum(c.get(NOT_RUN, 0) for c in self.counts.values()))
        info = {
            "pipeline": PIPELINE_NAME,
            "version": __version__,
            "module": self.module,
            "label": self.label,
            "date": datetime.now().isoformat(timespec="seconds"),
            "command": " ".join(sys.argv),
            "counts": self.counts,
        }
        info.update(run_info or {})
        with open(os.path.join(self.outdir, "run.json"), "w", encoding="utf-8") as fh:
            json.dump(info, fh, indent=2, sort_keys=True, default=str)
            fh.write("\n")
        with open(os.path.join(self.outdir, "summary.tsv"), "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            for row in summary_table(self.counts):
                w.writerow(row)
        note = citation_note(self.module, self.outdir, info.get("tools"), self.counts)
        with open(os.path.join(self.outdir, "citations.txt"), "w", encoding="utf-8") as fh:
            fh.write(note)
        if not quiet:
            sys.stderr.write("\n" + note)
            sys.stderr.flush()
        if not keep_work:
            shutil.rmtree(self.workdir, ignore_errors=True)

    def abort(self):
        """
        Give up on this run (for example its database could not be used): close the files and, if the output folder
        was created by this run, remove it, so header-only tables cannot be mistaken for a run on zero genes.
        A folder that already held an earlier run (--resume or --force) is left alone.
        """
        for fh in self._files:
            fh.close()
        self._files = []
        if self._fresh:
            shutil.rmtree(self.outdir, ignore_errors=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if self._files:
            for fh in self._files:
                fh.close()
            self._files = []
        return False


# --------------------------------------------------------------------------
# 6. runner
# --------------------------------------------------------------------------
class CommandError(Exception):
    pass


def run_command(cmd, log_path=None, cwd=None):
    """
    Run an external command (list of strings). stdout and stderr go to log_path if given.
    Returns the return code (127 if the program is not found). Never raises on non-zero exit.
    """
    fh = open(log_path, "ab") if log_path else None
    try:
        if fh:
            fh.write(("$ " + " ".join(cmd) + "\n").encode("utf-8"))
            fh.flush()
        try:
            proc = subprocess.run(cmd, stdout=fh if fh else subprocess.DEVNULL,
                                  stderr=subprocess.STDOUT, cwd=cwd)
        except FileNotFoundError:
            if fh:
                fh.write(("program not found: %s\n" % cmd[0]).encode("utf-8"))
            return 127
        return proc.returncode
    finally:
        if fh:
            fh.close()


def check_rc(rc, name):
    if rc != 0:
        raise CommandError("%s failed with exit code %d" % (name, rc))


def run_jobs(names, fn, workers=1, marker_dir=None, resume=False):
    """
    Run fn(name) for each name, `workers` at a time. fn must raise on failure.
    With resume=True and a marker_dir, jobs that already have '<name>.done' are skipped.
    Returns {"done": [...], "skipped": [...], "failed": {name: error message}}.
    A failed job is never silently dropped: callers must mark its genes NOT_RUN.
    """
    def marker(n):
        return os.path.join(marker_dir, n + ".done") if marker_dir else None

    if marker_dir:
        os.makedirs(marker_dir, exist_ok=True)
    done, skipped, failed = [], [], {}
    pending = []
    for n in names:
        if resume and marker_dir and os.path.exists(marker(n)):
            skipped.append(n)
        else:
            pending.append(n)
    if pending:
        with ThreadPoolExecutor(max_workers=max(1, int(workers))) as ex:
            futs = {ex.submit(fn, n): n for n in pending}
            for f in as_completed(futs):
                n = futs[f]
                try:
                    f.result()
                except Exception as e:  # noqa: BLE001 - report every failure, keep going
                    failed[n] = str(e)
                    continue
                done.append(n)
                if marker_dir:
                    open(marker(n), "w").close()
    return {"done": sorted(done), "skipped": sorted(skipped), "failed": failed}


BLAST_HEADER_RE = re.compile(r"^# [A-Z]*BLAST[A-Z]* \d")        # '# BLASTP 2.17.0+', '# RPSBLAST ...', '# TBLASTN ...'
BLAST_QUERY_RE = re.compile(r"^# Query:\s*(\S+)")
BLAST_FOOTER_RE = re.compile(r"^# BLAST processed (\d+) quer")


def read_blast_tabular(path, ncols, outfmt_hint=""):
    """
    Read BLAST+ tabular output written with '-outfmt 7' (blastp, tblastn, rpsblast, ...).
    Returns (confirmed, unconfirmed, rows): raw query IDs proven complete, raw query IDs NOT proven complete, and the
    alignment lines as lists of column strings.

    How completeness is proven (this is how BLAST+ really writes it): a '# <PROGRAM>' line and a '# Query:' line come
    before every query, and one '# BLAST processed N queries' footer comes at the end of a finished run. A group of
    queries is trusted only if a footer follows it AND the footer's N equals the number of queries seen since the previous
    footer. A file made by joining several runs (one footer each) is fine; if one run was cut short, the counts do not add
    up and every query since the last good footer is left unconfirmed, because the file alone cannot say which were lost.
    Callers must treat unconfirmed genes as NOT_RUN, never as "no hit".
    """
    confirmed, unconfirmed, pending, rows = set(), set(), [], []
    saw_header = False
    with open(path, "r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                if BLAST_HEADER_RE.match(line):
                    saw_header = True
                m = BLAST_QUERY_RE.match(line)
                if m:
                    pending.append(m.group(1))
                f = BLAST_FOOTER_RE.match(line)
                if f:
                    if len(pending) == int(f.group(1)):
                        confirmed.update(pending)
                    else:
                        unconfirmed.update(pending)
                    pending = []
                continue
            cols = line.split("\t")
            if len(cols) != ncols:
                raise HepError("%s line %d: expected %d tab-separated columns (outfmt '7 %s'), got %d"
                               % (path, lineno, ncols, outfmt_hint, len(cols)))
            rows.append(cols)
    if rows and not saw_header:
        raise HepError("%s: no BLAST header line found; is this an outfmt 7 file?" % path)
    unconfirmed.update(pending)               # the file ended before a completion footer
    unconfirmed -= confirmed
    return confirmed, unconfirmed, rows


def union_length(intervals):
    """Total length covered by inclusive (start, end) intervals, overlaps counted once; direction does not matter."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


# --------------------------------------------------------------------------
# 7. tool discovery and versions
# --------------------------------------------------------------------------
def find_tool(*names):
    """Return the path of the first program found on PATH, or None."""
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None


def tool_version(cmd, pattern=None):
    """Run a version command and return the first non-empty output line (or the first regex match)."""
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              universal_newlines=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        if pattern:
            m = re.search(pattern, line)
            if m:
                return m.group(0)
            continue
        return line
    return "unknown"


# --------------------------------------------------------------------------
# 8. command-line scaffolding shared by every module
# --------------------------------------------------------------------------
def _flag(key):
    return "--" + key.replace("_", "-")


def build_parser(module, description, defaults_dir=None, extra=None):
    """
    Argument parser with the shared options plus one flag per config key of the module
    (e.g. --evalue, --min-qcov, --chunk-size). Flags left unset fall back to the config resolution.
    `extra`, if given, is called with the parser to add module-specific options (e.g. --db).
    """
    defaults = load_defaults(module, defaults_dir)
    p = argparse.ArgumentParser(prog=module, description=description,
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version="%s %s" % (PIPELINE_NAME, __version__))
    p.add_argument("-i", "--input", required=True,
                   help="per-species protein FASTA file, or a directory of them")
    p.add_argument("-o", "--output", required=True, help="output directory for this run")
    p.add_argument("--threads", type=int, default=1, help="CPU threads per external job")
    p.add_argument("--label", default=None,
                   help="name for this run in combine (default: the module name), e.g. 'organelle'")
    p.add_argument("--config", default=None, help="config file (see presets/ for an example)")
    p.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                   help="override any setting; repeatable")
    p.add_argument("--species-map", default=None,
                   help="two-column file (species name, file basename) as in Stage 0")
    p.add_argument("--resume", action="store_true",
                   help="continue an interrupted run in an existing output folder, skipping chunks that finished")
    p.add_argument("--force", action="store_true",
                   help="overwrite the results already in the output folder (only if it holds this same run)")
    p.add_argument("--keep-work", action="store_true", help="keep the work/ directory")
    p.add_argument("--print-config", action="store_true",
                   help="print the resolved settings and exit")
    p.add_argument("--quiet", action="store_true",
                   help="do not print the end-of-run note (citations.txt is still written)")
    if extra:
        extra(p)
    grp = p.add_argument_group("settings (each overrides the config file and defaults)")
    for key, val in sorted(defaults.items()):
        grp.add_argument(_flag(key), dest="cfg_" + key, default=None, metavar="VALUE",
                         help="default: %s" % ("(unset)" if val is None else val))
    p.set_defaults(_defaults_keys=sorted(defaults))
    return p


def collect_overrides(args):
    """Merge per-key flags and --set KEY=VALUE into one dict."""
    overrides = {}
    for key in args._defaults_keys:
        v = getattr(args, "cfg_" + key, None)
        if v is not None:
            overrides[key] = v
    for item in args.set:
        if "=" not in item:
            raise HepError("--set expects KEY=VALUE, got: %s" % item)
        k, _, v = item.partition("=")
        overrides[k.strip()] = v.strip()
    return overrides


def start_module(module, description, argv=None, defaults_dir=None, extra=None):
    """
    Common start-up for a module script. Returns a dict with:
      args, cfg, non_default, inputs [(species, path)], species_map.
    Handles --print-config (prints and exits). Warns about every non-default setting.
    `extra` adds module-specific command-line options (see build_parser).
    """
    parser = build_parser(module, description, defaults_dir, extra)
    args = parser.parse_args(argv)
    try:
        cfg, non_default = resolve_config(module, args.config, collect_overrides(args), defaults_dir)
        if args.print_config:
            for k in sorted(cfg):
                print("%s=%s" % (k, "" if cfg[k] is None else cfg[k]))
            sys.exit(0)
        species_map = load_species_map(args.species_map) if args.species_map else None
        inputs = discover_inputs(args.input, as_list(cfg.get("species_strip_suffixes")), species_map)
    except HepError as e:
        die(str(e))
    for k, d in sorted(non_default.items()):
        warn("non-default setting %s=%s (default %s)" % (k, d["value"], d["default"]))
    return {"args": args, "cfg": cfg, "non_default": non_default,
            "inputs": inputs, "species_map": species_map}


def run_info_base(ctx):
    """Fields every module records in run.json (threshold values, non-default flags, inputs)."""
    return {
        "settings": ctx["cfg"],
        "non_default_settings": ctx["non_default"],
        "config_file": ctx["args"].config,
        "inputs": [{"species": sp, "file": os.path.abspath(p)} for sp, p in ctx["inputs"]],
    }


if __name__ == "__main__":
    sys.stderr.write("hep_common.py is a helper module; run one of the module scripts instead.\n")
    sys.exit(2)
