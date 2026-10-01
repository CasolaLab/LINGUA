#!/usr/bin/env python3
"""
combine.py - merge the results of several runs into one per-gene table and a final verdict.

Each run (domain-interpro, domain-cdd, homology-blast, homology-jackhmmer) leaves a folder. combine reads those folders and
answers, for every gene of every species: which runs removed it, which never searched it, and does it pass everything?

Per gene the verdict is
  REMOVED     at least one run removed it (a counted hit)
  INCOMPLETE  no run removed it, but at least one run did not finish or never searched it, so a pass cannot be claimed
  PASS        every run searched it and none removed it

A gene that an earlier run removed is not searched by the runs after it (iterative filtering). Its later columns say
NOT_SEARCHED, and `dropped_by` names the run that stopped it. It is never written as 0. Nothing is ever left blank:
a value that does not exist is NA.

Nothing here uses the network. Python 3.11+, standard library only.

  python3 combine.py --runs out/domain-cdd out/homology-blast out/homology-jackhmmer -o out/combine
  python3 combine.py --runs out/ -o out/combine          # every run folder inside out/, in chain order
"""

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402

MODULE = "combine"
NOT_SEARCHED = "NOT_SEARCHED"
PASS, REMOVED, INCOMPLETE = "PASS", "REMOVED", "INCOMPLETE"
FIXED_COLUMNS = ["species", "gene_id", "verdict", "removed_by", "dropped_by", "n_runs_removed", "n_runs_searched",
                 "n_hit_species", "n_hit_genera"]
PASSING_STATUSES = (h.NO_HIT, h.EXCLUDED_ONLY, h.IN_PHYLOGENY_ONLY)   # always in passed_ids.tsv
FAILING_STATUSES = (h.HIT, h.NOT_RUN)                                  # never in passed_ids.tsv


class Run(object):
    """One finished run folder, read from the files every module writes."""

    def __init__(self, path):
        self.path = path
        marker = os.path.join(path, h.RUN_MARKER)
        if not os.path.isfile(marker):
            raise h.HepError("%s is not a run folder of this tool (no %s file)" % (path, h.RUN_MARKER))
        mk = {}
        with open(marker, "r", encoding="utf-8") as fh:
            for line in fh:
                k, _, v = line.strip().partition("=")
                mk[k] = v
        self.module, self.label = mk.get("module"), mk.get("label")
        if self.module == MODULE:
            raise h.HepError("%s is the output of combine itself, not a search run" % path)
        for name in ("gene_status.tsv", "passed_ids.tsv", "hits.tsv", "run.json"):
            if not os.path.isfile(os.path.join(path, name)):
                raise h.HepError("%s: %s is missing (was the run finished?)" % (path, name))
        with open(os.path.join(path, "run.json"), "r", encoding="utf-8") as fh:
            self.info = json.load(fh)
        chain = self.info.get("chain") or {}
        self.position = chain.get("position")
        self.status = {}
        for r in self._rows("gene_status.tsv", ["species", "gene_id", "status"]):
            key = (r["species"], r["gene_id"])
            if key in self.status:
                raise h.HepError("%s: gene %s of %s is listed twice in gene_status.tsv" % (path, r["gene_id"], r["species"]))
            self.status[key] = r["status"]
        self.passed = {(r["species"], r["gene_id"]) for r in self._rows("passed_ids.tsv", ["species", "gene_id"])}
        self.hit_species = {}       # (species, gene) -> set of target species of COUNTED hits
        for r in self._rows("hits.tsv", ["species", "gene_id", "target_species", "counts_as_hit"]):
            if r["counts_as_hit"] == "yes" and r["target_species"] != "NA":
                self.hit_species.setdefault((r["species"], r["gene_id"]), set()).add(r["target_species"])
        self._check_consistent()

    def _rows(self, name, need):
        p = os.path.join(self.path, name)
        with open(p, "r", encoding="utf-8", newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            head = next(rd, None)
            if head is None or any(c not in head for c in need):
                raise h.HepError("%s: %s must have the columns %s (found: %s)" % (self.path, name, ", ".join(need), head))
            for row in rd:
                if row:
                    yield dict(zip(head, row))

    def _check_consistent(self):
        """The two files of a run must agree; if they do not, nothing can be trusted."""
        for key, st in self.status.items():
            if st in PASSING_STATUSES and key not in self.passed:
                raise h.HepError("%s: %s is %s in gene_status.tsv but missing from passed_ids.tsv (run files disagree)"
                                 % (self.path, key[1], st))
            if st in FAILING_STATUSES and key in self.passed:
                raise h.HepError("%s: %s is %s in gene_status.tsv but listed in passed_ids.tsv (run files disagree)"
                                 % (self.path, key[1], st))
        for key in self.passed:
            if key not in self.status:
                raise h.HepError("%s: %s is in passed_ids.tsv but not in gene_status.tsv (run files disagree)" % (self.path, key[1]))

    def removed(self, key):
        """Did this run remove the gene? True / False, or None if it did not decide (not present, or NOT_RUN)."""
        st = self.status.get(key)
        if st is None or st == h.NOT_RUN:
            return None
        return key not in self.passed

    @property
    def is_homology(self):
        return str(self.module).startswith("homology")


def expand_runs(paths):
    """A path that is a run folder is taken as is; a folder of run folders is expanded: --sequential runs by chain position, the rest by the date they were made."""
    out = []
    for p in paths:
        if not os.path.isdir(p):
            raise h.HepError("not a folder: %s" % p)
        if os.path.isfile(os.path.join(p, h.RUN_MARKER)):
            out.append(p)
            continue
        kids = sorted(os.path.join(p, d) for d in os.listdir(p) if os.path.isfile(os.path.join(p, d, h.RUN_MARKER)))
        if not kids:
            raise h.HepError("%s holds no run folders (no %s files inside)" % (p, h.RUN_MARKER))
        runs = [Run(k) for k in kids]
        # --sequential runs carry their position; the others are ordered by the date they were made (a chained run is always
        # made after the run it follows), then by name
        runs.sort(key=lambda r: (r.position if r.position is not None else 10 ** 6, str(r.info.get("date") or ""), r.path))
        out += [r.path for r in runs]
    return out


def combine(runs, na="NA"):
    """
    Build the per-gene rows. `runs` is the list of Run objects in chain order. Returns (rows, warnings), each row a dict.
    """
    universe = set()
    for r in runs:
        universe |= set(r.status)
    rows, warnings = [], []
    for key in sorted(universe):
        sp, gene = key
        cols, removed_by, stopper, not_finished, searched = {}, [], None, False, 0
        for i, r in enumerate(runs):
            st = r.status.get(key)
            if st is not None:
                cols[r.label] = st
                dec = r.removed(key)
                if dec is None:
                    not_finished = True
                else:
                    searched += 1
                    if dec:
                        removed_by.append(r.label)
                if stopper is None and (dec or dec is None):
                    stopper = r.label
            else:
                cols[r.label] = NOT_SEARCHED
                if stopper is None:
                    # nothing earlier explains why this run never received the gene
                    not_finished = True
                    warnings.append("%s (%s) is missing from run '%s' and no earlier run stopped it" % (gene, sp, r.label))
                else:
                    earlier = next(x for x in runs if x.label == stopper)
                    if earlier.removed(key) is None:
                        not_finished = True     # the gene was stopped by a run that did not finish
        if removed_by:
            verdict = REMOVED
        elif not_finished:
            verdict = INCOMPLETE
        else:
            verdict = PASS
        # breadth of homologs: distinct species and genera among counted hits, from runs that searched it
        homology_searched = [r for r in runs if r.is_homology and r.removed(key) is not None]
        if homology_searched:
            names = set()
            for r in homology_searched:
                names |= r.hit_species.get(key, set())
            n_sp = len({h.taxon_key(n, "species") for n in names})
            n_gen = len({h.taxon_key(n, "genus") for n in names})
        else:
            n_sp = n_gen = na
        row = {"species": sp, "gene_id": gene, "verdict": verdict,
               "removed_by": ";".join(removed_by) if removed_by else na,
               "dropped_by": stopper if stopper else na, "n_runs_removed": len(removed_by), "n_runs_searched": searched,
               "n_hit_species": n_sp, "n_hit_genera": n_gen, "_cols": cols, "_removed_by": removed_by}
        rows.append(row)
    return rows, warnings


def run_summary(runs, rows):
    out = [["run", "module", "position", "genes_received", "removed", "passed", "not_run", "not_searched", "removed_first"]]
    for i, r in enumerate(runs, 1):
        recv = [k for k in r.status]
        removed = [k for k in recv if r.removed(k)]
        first = [k for k in removed if not any(x.removed(k) for x in runs[:i - 1])]
        total = len({(x["species"], x["gene_id"]) for x in rows})
        out.append([r.label, r.module, i, len(recv), len(removed), len([k for k in recv if r.removed(k) is False]),
                    len([k for k in recv if r.removed(k) is None]), total - len(recv), len(first)])
    return out


def overlaps(rows):
    counts = {}
    for r in rows:
        if r["verdict"] == REMOVED:
            key = "+".join(r["_removed_by"])
            counts[key] = counts.get(key, 0) + 1
    out = [["removed_by", "n_genes"]]
    for k, v in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        out.append([k, v])
    return out


def verdict_summary(rows):
    per = {}
    for r in rows:
        d = per.setdefault(r["species"], {PASS: 0, REMOVED: 0, INCOMPLETE: 0})
        d[r["verdict"]] += 1
    out = [["species", "genes", PASS, REMOVED, INCOMPLETE]]
    tot = {PASS: 0, REMOVED: 0, INCOMPLETE: 0}
    for sp in sorted(per):
        d = per[sp]
        out.append([sp, sum(d.values()), d[PASS], d[REMOVED], d[INCOMPLETE]])
        for k in tot:
            tot[k] += d[k]
    out.append(["ALL", sum(tot.values()), tot[PASS], tot[REMOVED], tot[INCOMPLETE]])
    return out


def write_tsv(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for row in rows:
            w.writerow(["NA" if (c is None or c == "") else c for c in row])


def collect_pass_fasta(runs, rows):
    """
    Sequences of the PASS genes, with their original headers, taken from the runs' passed/ folders (a gene that passed every
    run is in every run's passed/ folder that received it). Returns {species: [(header, seq)]}. A gene whose sequence cannot
    be found is an error: nothing is silently left out.
    """
    want = {}
    for r in rows:
        if r["verdict"] == PASS:
            want.setdefault(r["species"], set()).add(r["gene_id"])
    found = {}
    for sp, genes in want.items():
        got = {}
        for run in reversed(runs):
            p = os.path.join(run.path, "passed", sp + ".faa")
            if not os.path.isfile(p):
                continue
            for rid, header, seq in h.read_fasta(p):
                g = h.canonical_id(rid)
                if g in genes and g not in got:
                    got[g] = (header, seq)
            if len(got) == len(genes):
                break
        missing = sorted(genes - set(got))
        if missing:
            raise h.HepError("the sequences of %d passed gene(s) of %s were not found in any run's passed/ folder (e.g. %s)"
                             % (len(missing), sp, missing[0]))
        found[sp] = [got[g] for g in sorted(genes)]
    return found


def main(argv=None):
    ap = argparse.ArgumentParser(prog="combine", description=__doc__.strip().split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version="%s %s" % (h.PIPELINE_NAME, h.__version__))
    ap.add_argument("--runs", nargs="+", required=True, metavar="DIR",
                    help="run folders, in the order the runs were made (the chain order for iterative filtering), or a "
                         "folder that contains run folders")
    ap.add_argument("-o", "--output", required=True, help="output folder")
    ap.add_argument("--binary", action="store_true",
                    help="also write gene_matrix_binary.tsv: 1 = the run removed the gene, 0 = the run searched it and did "
                         "not remove it, and the --na value where the run did not search it or did not finish")
    ap.add_argument("--na", default="NA", help="what to write in place of a missing value in the tables (default NA; never blank)")
    ap.add_argument("--force", action="store_true", help="overwrite an earlier combine result in the output folder")
    ap.add_argument("--quiet", action="store_true", help="do not print the end-of-run note")
    a = ap.parse_args(argv)
    if not a.na.strip():
        h.die("--na must not be blank: a blank cell is exactly what this tool avoids")
    try:
        paths = expand_runs(a.runs)
        runs = [Run(p) for p in paths]
        labels = [r.label for r in runs]
        if len(set(labels)) != len(labels):
            dup = sorted({l for l in labels if labels.count(l) > 1})
            raise h.HepError("two runs have the same label (%s); give each run its own --label when you make it" % ", ".join(dup))
        bad = [l for l in labels if l in FIXED_COLUMNS]
        if bad:
            raise h.HepError("a run label cannot be %s: it is a column name in the combined table" % ", ".join(bad))
        h.check_output_folder(a.output, MODULE, MODULE, resume=False, force=a.force)
        rows, warnings = combine(runs, a.na)
        pass_fasta = collect_pass_fasta(runs, rows)
    except h.HepError as e:
        h.die(str(e))
    os.makedirs(a.output, exist_ok=True)
    with open(os.path.join(a.output, h.RUN_MARKER), "w", encoding="utf-8") as fh:
        fh.write("module=%s\nlabel=%s\n" % (MODULE, MODULE))
    write_tsv(os.path.join(a.output, "gene_matrix.tsv"),
              [["species", "gene_id"] + labels + FIXED_COLUMNS[2:]] +
              [[r["species"], r["gene_id"]] + [r["_cols"][l] for l in labels] + [r[c] for c in FIXED_COLUMNS[2:]] for r in rows])
    if a.binary:
        by_label = {r.label: r for r in runs}
        brow = [["species", "gene_id"] + labels]
        for r in rows:
            key = (r["species"], r["gene_id"])
            brow.append([r["species"], r["gene_id"]] + [
                {True: 1, False: 0, None: a.na}[by_label[l].removed(key)] for l in labels])
        write_tsv(os.path.join(a.output, "gene_matrix_binary.tsv"), brow)
    write_tsv(os.path.join(a.output, "run_summary.tsv"), run_summary(runs, rows))
    write_tsv(os.path.join(a.output, "overlaps.tsv"), overlaps(rows))
    summ = verdict_summary(rows)
    write_tsv(os.path.join(a.output, "summary.tsv"), summ)
    write_tsv(os.path.join(a.output, "passed_all_ids.tsv"),
              [["species", "gene_id"]] + [[r["species"], r["gene_id"]] for r in rows if r["verdict"] == PASS])
    pdir = os.path.join(a.output, "passed_all")
    os.makedirs(pdir, exist_ok=True)
    for sp, recs in sorted(pass_fasta.items()):
        with open(os.path.join(pdir, sp + ".faa"), "w", encoding="utf-8") as fh:
            for header, seq in recs:
                h.write_fasta_record(fh, header, seq)
    tools = []
    seen = set()
    for r in runs:
        for t in r.info.get("tools") or []:
            k = (t.get("name"), t.get("version"))
            if k not in seen:
                seen.add(k)
                tools.append(t)
    info = {"module": MODULE, "pipeline": h.PIPELINE_NAME, "version": h.__version__,
            "runs": [{"label": r.label, "module": r.module, "path": os.path.abspath(r.path),
                      "tools": r.info.get("tools"), "date": r.info.get("date")} for r in runs],
            "run_order": labels, "na": a.na, "warnings": warnings,
            "counts": {"genes": summ[-1][1], PASS: summ[-1][2], REMOVED: summ[-1][3], INCOMPLETE: summ[-1][4]}}
    with open(os.path.join(a.output, "run.json"), "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2, sort_keys=True, default=str)
        fh.write("\n")
    counts = {row[0]: {"input": row[1], "passed": row[2]} for row in summ[1:-1]}
    note = h.citation_note(MODULE, a.output, tools, counts)
    with open(os.path.join(a.output, "citations.txt"), "w", encoding="utf-8") as fh:
        fh.write(note)
    for w in warnings[:10]:
        h.warn(w)
    if len(warnings) > 10:
        h.warn("... and %d more (all in run.json)" % (len(warnings) - 10))
    h.note_incomplete(summ[-1][4])       # genes whose pass cannot be claimed
    if not a.quiet:
        sys.stderr.write("\n" + note)


if __name__ == "__main__":
    main()
    sys.exit(h.exit_status())
