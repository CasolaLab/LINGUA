#!/usr/bin/env python3
"""
make_homology_db.py - build the species-tagged FASTA that homology-blast and homology-jackhmmer use as a database.

Input: a folder of per-species protein FASTA files (or single files). Output: one FASTA in which every ID is
<ID>__<Species_Name> (IDs that already carry a tag, as Stage 1 --add-species writes, keep it), plus a species table
(<output>.species.tsv). Analysis species can be left out. Nothing is downloaded or installed; no BLAST index is made
here (homology-blast builds and caches its own).
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402

TABLE_COLUMNS = ["species", "source_file", "n_sequences", "n_ids_already_tagged", "n_ids_retagged", "status"]


def clean_species(name):
    return re.sub(r"\s+", "_", str(name).strip())


def gather_files(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            found = sorted(os.path.join(p, f) for f in os.listdir(p)
                           if f.lower().endswith(h.FASTA_EXTS) and os.path.isfile(os.path.join(p, f)))
            if not found:
                raise h.HepError("no FASTA files (%s) in %s" % ("/".join(h.FASTA_EXTS), p))
            files += found
        elif os.path.isfile(p):
            files.append(p)
        else:
            raise h.HepError("input not found: %s" % p)
    return files


def excluded_species(args, species_map=None):
    names = list(args.exclude_species or [])
    for d in args.exclude_inputs or []:
        names += [sp for sp, _f in h.discover_inputs(d, species_map=species_map)]
    return names


def build(files, out_path, exclude, level, species_map=None):
    """Write the database. Returns the table rows. Raises HepError before writing anything if IDs would collide."""
    matcher = h.TaxonMatcher(exclude, (), level)
    rows, seen, plan = [], {}, []
    for f in files:
        sp = clean_species(h.species_from_path(f, species_map=species_map))
        hit = matcher.match(sp)
        if hit:
            rows.append([sp, f, "0", "0", "0", "EXCLUDED (matches %s)" % hit])
            continue
        plan.append((sp, f))
    used = {}
    for sp, f in plan:
        if sp in used:
            raise h.HepError("two files give species '%s': %s and %s" % (sp, used[sp], f))
        used[sp] = f
    tmp = out_path + ".partial"
    n_total = 0
    try:
        with open(tmp, "w", encoding="utf-8") as out:
            for sp, f in plan:
                n = tagged = odd = 0
                for rid, header, seq in h.read_fasta(f):
                    gene, tag = h.split_id(rid)
                    if not gene:
                        raise h.HepError("%s: record with an empty ID" % f)
                    if tag is not None and h.looks_like_species(tag):
                        new_id, tagged = rid, tagged + 1
                    else:
                        # no tag, or a '__' inside a random-looking ID (not a species name): add our own tag
                        new_id = "%s__%s" % (rid, sp)
                        if tag is not None:
                            odd += 1
                    if new_id in seen:
                        raise h.HepError("duplicate ID '%s' in %s and %s" % (new_id, seen[new_id], f))
                    seen[new_id] = f
                    seq = seq.replace("*", "")
                    if not seq:
                        raise h.HepError("%s: '%s' has no sequence" % (f, rid))
                    rest = header[len(rid):]
                    h.write_fasta_record(out, new_id + rest, seq)
                    n += 1
                if n == 0:
                    raise h.HepError("%s: no sequences" % f)
                rows.append([sp, f, str(n), str(tagged), str(odd), "INCLUDED"])
                n_total += n
        if n_total == 0:
            raise h.HepError("nothing left to write: every input species was excluded")
        os.replace(tmp, out_path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return sorted(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="make-db", description=__doc__.strip().split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version="%s %s" % (h.PIPELINE_NAME, h.__version__))
    ap.add_argument("-i", "--input", nargs="+", required=True, help="folder(s) or file(s) of per-species protein FASTA")
    ap.add_argument("-o", "--output", required=True, help="tagged FASTA to write")
    ap.add_argument("--exclude-inputs", nargs="+", metavar="DIR",
                    help="leave out the species of the FASTA files in these folders (your analysis species)")
    ap.add_argument("--exclude-species", nargs="+", metavar="NAME", help="leave out these species")
    ap.add_argument("--match-level", default="species", choices=h.MATCH_LEVELS,
                    help="how species are compared for exclusion (default species: subspecies count as the same)")
    ap.add_argument("--species-map", help="two columns (species name, file basename), a header row, tab- or comma-separated: the same map Stage 1 writes; applies to the database files and to --exclude-inputs")
    ap.add_argument("--force", action="store_true", help="overwrite an existing output")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    try:
        if os.path.exists(a.output) and not a.force:
            raise h.HepError("%s exists; use --force to overwrite" % a.output)
        d = os.path.dirname(os.path.abspath(a.output))
        os.makedirs(d, exist_ok=True)
        smap = h.load_species_map(a.species_map) if a.species_map else None
        rows = build(gather_files(a.input), a.output, excluded_species(a, smap), a.match_level, smap)
        with open(a.output + ".species.tsv", "w", encoding="utf-8") as fh:
            fh.write("\t".join(TABLE_COLUMNS) + "\n")
            for r in rows:
                fh.write("\t".join(x if x != "" else "NA" for x in r) + "\n")
    except h.HepError as e:
        sys.stderr.write("ERROR: %s\n" % e)
        return 1
    if not a.quiet:
        inc = [r for r in rows if r[5] == "INCLUDED"]
        h.log("wrote %s: %d species, %d sequences (%d species excluded)" % (
            a.output, len(inc), sum(int(r[2]) for r in inc), len(rows) - len(inc)))
        h.log("species table: %s.species.tsv" % a.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
