#!/usr/bin/env python3
"""
make_cds_bed.py - given a GFF3 annotation and a list of gene/transcript/protein IDs, write a BED of their CDS
coordinates and a summary table (exon and CDS segment counts, gene/mRNA/CDS lengths).

This is a general-purpose coordinate lookup, not tied to any one part of this tool: the ID list can be
`passed_all_ids.tsv` from `combine`, a full candidate set, or any list of IDs you have, from any source. It only
needs a GFF3 file and an ID list; no FASTA, no search tools.

Matching: an ID is looked up against each mRNA's own ID, Name, parent gene ID and protein_id attributes; if that
finds nothing, each CDS line is matched directly by its own protein_id, locus_tag, Name or ID attribute (some GFFs,
especially from smaller or non-model genomes, do not carry a clean Name on the mRNA line). An ID that matches
nothing is reported, not silently dropped; a batch run counts as INCOMPLETE (exit status 2) if any ID across any
species was not found, matching the rest of this tool's rule that a missing result is never mistaken for "no match".

Nothing here uses the network. Python 3.11+, standard library only.

  single:  python3 make_cds_bed.py --ids passed_all_ids.tsv --gff annotation.gff3 -o out/
  batch:   python3 make_cds_bed.py --ids-dir ids/ --gff-dir gff/ -o out/ --species-map species.tsv
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hep_common as h  # noqa: E402

MODULE = "make-bed"
ID_EXTS = (".txt", ".csv", ".tsv", ".list")
GFF_EXTS = (".gff", ".gff3")
BED_COLUMNS = ["#chrom", "chromStart", "chromEnd", "name", "score", "strand", "id"]
MASTER_COLUMNS = ["species", "gene_id", "matched_id", "n_exons", "n_cds_segments", "gene_length", "mrna_length", "cds_length"]
NOT_FOUND_COLUMNS = ["species", "gene_id"]


# --------------------------------------------------------------------------
# 1. ID list loading
# --------------------------------------------------------------------------
def normalize_id(identifier):
    """Strip a trailing '.p' (some pipelines mark protein-translation IDs this way). Never returns None."""
    if not identifier:
        return ""
    identifier = identifier.strip()
    return identifier[:-2] if identifier.endswith(".p") else identifier


def load_target_ids(path, id_column=None):
    """
    Read an ID list: a plain one-ID-per-line file, or a delimited table with a header (comma or tab detected from
    the first line). Returns (ids, raw_count): the set of IDs to look for (raw plus normalized), and how many rows
    were read. Raises HepError on a file that cannot be read at all.
    """
    if not os.path.isfile(path):
        raise h.HepError("ID list not found: %s" % path)
    try:
        with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
            text = fh.read()
    except OSError as e:
        raise h.HepError("cannot read %s: %s" % (path, e))
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return set(), 0
    delim = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    ids, raw_count = set(), 0
    if delim:
        reader = csv.DictReader(lines, delimiter=delim)
        fields = reader.fieldnames or []
        col = id_column
        if col and col not in fields:
            raise h.HepError("%s has no column '%s'. Columns found: %s" % (path, col, ", ".join(fields)))
        if not col:
            candidates = [c for c in fields if re.search(r"gene|protein|transcript|\bid\b", c, re.I)]
            col = candidates[0] if candidates else fields[0]
        for row in reader:
            val = (row.get(col) or "").strip()
            if not val:
                continue
            raw_count += 1
            clean = val.lstrip(">").split()[0]
            ids.add(clean)
            ids.add(normalize_id(clean))
    else:
        for line in lines:
            clean = line.strip().lstrip(">").split()[0]
            if clean:
                raw_count += 1
                ids.add(clean)
                ids.add(normalize_id(clean))
    return ids, raw_count


# --------------------------------------------------------------------------
# 2. GFF pass 1: gene lengths, mRNA -> gene links, and which mRNAs match a target
# --------------------------------------------------------------------------
RE_ID = re.compile(r"ID=([^;]+)")
RE_PARENT = re.compile(r"Parent=([^;]+)")
RE_NAME = re.compile(r"Name=([^;]+)")
RE_PROTEIN = re.compile(r"protein_id=([^;]+)")
RE_LOCUS = re.compile(r"locus_tag=([^;]+)")


def gff_lines(gff_path):
    with open(gff_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            yield parts


def build_mrna_map(gff_path, target_ids):
    """
    Pass 1 over the GFF. Returns:
      approved   {mrna_id: {"matched_id", "gene_id", "mrna_length"}}  for mRNAs whose own ID/Name/gene/protein_id
                 attribute is one of the targets
      parent_of  {mrna_id: gene_id}  for every mRNA, matched or not (used to resolve Gene_ID for CDS-fallback matches)
      gene_len   {gene_id: length}   the real span of each 'gene' feature
    """
    approved, parent_of, gene_len = {}, {}, {}
    for parts in gff_lines(gff_path):
        feat, attrs = parts[2], parts[8]
        if feat == "gene":
            m = RE_ID.search(attrs)
            if m:
                try:
                    gene_len[m.group(1)] = int(parts[4]) - int(parts[3]) + 1
                except ValueError:
                    pass
            continue
        if feat not in ("mRNA", "transcript", "mrna"):
            continue
        m_id = RE_ID.search(attrs)
        if not m_id:
            continue
        mrna_id = m_id.group(1)
        m_parent = RE_PARENT.search(attrs)
        gene_id = m_parent.group(1) if m_parent else "NA"
        if gene_id != "NA":
            parent_of[mrna_id] = gene_id
        try:
            mrna_len = int(parts[4]) - int(parts[3]) + 1
        except ValueError:
            mrna_len = 0
        keys = [mrna_id]
        m_name = RE_NAME.search(attrs)
        if m_name:
            keys.append(m_name.group(1))
        m_prot = RE_PROTEIN.search(attrs)
        if m_prot:
            keys.append(m_prot.group(1))
        if gene_id != "NA":
            keys.append(gene_id)
        matched = None
        for k in keys:
            if k in target_ids or normalize_id(k) in target_ids:
                matched = k if k in target_ids else normalize_id(k)
                break
        if matched:
            approved[mrna_id] = {"matched_id": matched, "gene_id": gene_id, "mrna_length": mrna_len}
    return approved, parent_of, gene_len


# --------------------------------------------------------------------------
# 3. GFF pass 2: write the BED and tally per-target metrics
# --------------------------------------------------------------------------
def extract(gff_path, approved, parent_of, gene_len, target_ids, out_bed):
    """
    Pass 2 over the GFF. Writes `out_bed` (one row per CDS segment, deduplicated on exact coordinates).
    Returns (n_bed_rows, matched_ids, rows): the BED row count, the set of target IDs that were found, and the
    per-matched-ID metric rows (dicts with keys matching MASTER_COLUMNS minus 'species').
    """
    seen_coords = set()
    n_bed_rows = 0
    metrics = defaultdict(lambda: {"gene_id": "NA", "n_exons": 0, "n_cds_segments": 0, "mrna_length": 0,
                                   "cds_length": 0, "min_start": None, "max_end": None})

    with open(out_bed, "w", encoding="utf-8", newline="") as bed_fh:
        w = csv.writer(bed_fh, delimiter="\t", lineterminator="\n")
        w.writerow(BED_COLUMNS)
        for parts in gff_lines(gff_path):
            feat = parts[2]
            if feat not in ("CDS", "exon"):
                continue
            attrs = parts[8]
            matched_id, gene_id, mrna_len = None, "NA", 0

            m_parent = RE_PARENT.search(attrs)
            parent_id = m_parent.group(1) if m_parent else None
            if parent_id:
                for par in parent_id.split(","):
                    info = approved.get(par)
                    if info:
                        matched_id, gene_id, mrna_len = info["matched_id"], info["gene_id"], info["mrna_length"]
                        break

            if not matched_id and feat == "CDS":
                # Fallback: some GFFs (smaller or non-model genomes) don't carry a usable Name on the mRNA line, so
                # the CDS line's own identifying attributes are tried directly against the target IDs.
                for m in (RE_PROTEIN.search(attrs), RE_LOCUS.search(attrs), RE_NAME.search(attrs), RE_ID.search(attrs)):
                    if not m:
                        continue
                    cand = m.group(1)
                    if cand in target_ids:
                        matched_id = cand
                    elif normalize_id(cand) in target_ids:
                        matched_id = normalize_id(cand)
                    if matched_id:
                        break
                if matched_id and parent_id:
                    gene_id = parent_of.get(parent_id, "NA")
                    if gene_id == "NA":
                        h.warn("%s: found '%s' by its own attributes, but its parent mRNA '%s' names no gene"
                              % (os.path.basename(gff_path), matched_id, parent_id))

            if not matched_id:
                continue
            m_data = metrics[matched_id]
            if gene_id != "NA":
                m_data["gene_id"] = gene_id
            if mrna_len > 0:
                m_data["mrna_length"] = mrna_len

            if feat == "exon":
                m_data["n_exons"] += 1
                continue

            # feat == "CDS"
            try:
                gff_start, gff_end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            key = (parts[0], gff_start, gff_end, parts[6], matched_id)
            if key not in seen_coords:
                seen_coords.add(key)
                bed_start, bed_end = gff_start - 1, gff_end
                name = "%s__CDS_%d_%d" % (matched_id, bed_start, bed_end)
                w.writerow([parts[0], bed_start, bed_end, name, "0", parts[6], matched_id])
                n_bed_rows += 1
                # Length, segment count and span are per DISTINCT coordinate, so a duplicate CDS line in the GFF
                # (seen in real annotation exports) is not double-counted here even though it is skipped above.
                m_data["n_cds_segments"] += 1
                m_data["cds_length"] += gff_end - gff_start + 1
                m_data["min_start"] = gff_start if m_data["min_start"] is None else min(m_data["min_start"], gff_start)
                m_data["max_end"] = gff_end if m_data["max_end"] is None else max(m_data["max_end"], gff_end)

    rows = []
    for matched_id, d in metrics.items():
        mrna_length = d["mrna_length"]
        if mrna_length == 0 and d["min_start"] is not None:
            mrna_length = d["max_end"] - d["min_start"] + 1
        rows.append({
            "gene_id": d["gene_id"], "matched_id": matched_id, "n_exons": d["n_exons"],
            "n_cds_segments": d["n_cds_segments"], "gene_length": gene_len.get(d["gene_id"], "NA"),
            "mrna_length": mrna_length, "cds_length": d["cds_length"],
        })
    return n_bed_rows, set(metrics), rows


# --------------------------------------------------------------------------
# 4. One species
# --------------------------------------------------------------------------
def process_one(species, id_path, gff_path, out_bed, id_column):
    """
    Returns (stats, rows, not_found_rows). A species this cannot even attempt (its ID list could not be read, or its
    GFF is missing) is reported with h.warn() immediately -- never only as a status string a --quiet run could miss
    -- and its `status` starts with 'error:', which the caller treats as incomplete for the exit status, the same
    as an ID that was searched for and not found.
    """
    stats = {"species": species, "n_input_ids": 0, "n_found": 0, "n_cds_rows": 0, "status": "error"}
    try:
        targets, raw_count = load_target_ids(id_path, id_column)
    except h.HepError as e:
        stats["status"] = "error: %s" % e
        h.warn("%s: %s" % (species, e))
        return stats, [], []
    stats["n_input_ids"] = raw_count
    if not targets:
        stats["status"] = "error: no IDs read from %s" % id_path
        h.warn("%s: no IDs read from %s" % (species, id_path))
        return stats, [], []
    if not os.path.isfile(gff_path):
        stats["status"] = "error: GFF not found: %s" % gff_path
        h.warn("%s: GFF not found: %s" % (species, gff_path))
        return stats, [], []

    approved, parent_of, gene_len = build_mrna_map(gff_path, targets)
    n_bed_rows, matched_ids, rows = extract(gff_path, approved, parent_of, gene_len, targets, out_bed)
    for r in rows:
        r["species"] = species

    # a target may have been added twice (raw + normalized form); count distinct genes actually resolved
    canonical_targets = {normalize_id(t) for t in targets}
    canonical_matched = {normalize_id(m) for m in matched_ids}
    not_found = sorted(canonical_targets - canonical_matched)

    stats["n_found"] = len(canonical_targets) - len(not_found)
    stats["n_cds_rows"] = n_bed_rows
    stats["status"] = "ok" if not not_found else "ok (%d not found)" % len(not_found)
    return stats, rows, [{"species": species, "gene_id": g} for g in not_found]


# --------------------------------------------------------------------------
# 5. Matching an IDs folder to a GFF folder (batch mode)
# --------------------------------------------------------------------------
def match_files(ids_dir, gff_dir, gff_map=None):
    """
    Pair each ID-list file with a GFF file. The species name is the ID file's own stem with underscores turned back
    into spaces (matching the 'Species' column Stage 2 already writes, e.g. 'Arabidopsis_thaliana.tsv' ->
    'Arabidopsis thaliana').

    `gff_map` (from h.load_species_map: {gff_basename: species_name}, i.e. a species-map file naming the GFF files)
    gives an exact pairing by species name and is the reliable way to do this. Without it, a GFF is guessed by
    prefix from the ID file's name (its full stem, its first '_'-separated word, and, for a two-word
    'Genus_species' stem, the common 'Gspecies' abbreviation used by some annotation sources). Raises HepError on
    an ambiguous match (more than one GFF fits, or two ID files would claim the same GFF) rather than silently
    picking one.
    """
    id_files = sorted(f for f in os.listdir(ids_dir) if f.lower().endswith(ID_EXTS))
    gff_files = sorted(f for f in os.listdir(gff_dir) if f.lower().endswith(GFF_EXTS))
    gff_stems = {g.rsplit(".", 1)[0]: g for g in gff_files}
    species_to_gff_stem = {sp: base for base, sp in (gff_map or {}).items()}

    tasks, used_gff = {}, {}
    for f_id in id_files:
        stem = f_id.rsplit(".", 1)[0]
        species = stem.replace("_", " ")

        if gff_map:
            gff_stem = species_to_gff_stem.get(species)
            if not gff_stem or gff_stem not in gff_stems:
                h.warn("%s: no GFF for '%s' in --species-map" % (f_id, species))
                continue
            matches = [gff_stems[gff_stem]]
        else:
            parts = stem.split("_")
            candidates = [stem, parts[0]]
            if len(parts) >= 2:
                candidates.append(parts[0][0] + parts[1])
            matches = [g for g in gff_files if any(g.startswith(c) for c in candidates)]
            if not matches:
                h.warn("no GFF matches '%s' (tried: %s)" % (f_id, ", ".join(candidates)))
                continue
            if len(matches) > 1:
                raise h.HepError("'%s' matches more than one GFF file (%s); use --species-map to say which one"
                                 % (f_id, ", ".join(matches)))

        gff = matches[0]
        if gff in used_gff:
            raise h.HepError("both '%s' and '%s' matched the same GFF file '%s'; use --species-map to tell them apart"
                             % (used_gff[gff], f_id, gff))
        used_gff[gff] = f_id
        tasks[species] = (os.path.join(ids_dir, f_id), os.path.join(gff_dir, gff))
    return tasks


# --------------------------------------------------------------------------
# 6. Main
# --------------------------------------------------------------------------
def add_args(ap):
    single = ap.add_argument_group("single species")
    single.add_argument("--ids", metavar="FILE", help="one ID list file")
    single.add_argument("--gff", metavar="FILE", help="one GFF3 file")
    single.add_argument("--species-name", metavar="NAME", help="label for this species in the output (default: the "
                        "--ids file's name)")
    batch = ap.add_argument_group("several species")
    batch.add_argument("--ids-dir", metavar="DIR", help="folder of ID list files, one per species (species name = "
                       "file name, underscores read as spaces)")
    batch.add_argument("--gff-dir", metavar="DIR", help="folder of GFF3 files")
    ap.add_argument("--species-map", metavar="FILE", help="species name<TAB>GFF file basename (the same two-column "
                    "format used elsewhere in this tool); gives an exact ID-file-to-GFF pairing instead of one "
                    "guessed from file names")
    ap.add_argument("--id-column", metavar="NAME", help="column in the ID list to read (default: the first column "
                    "whose name contains 'gene', 'protein', 'transcript' or 'id')")
    ap.add_argument("--threads", type=int, default=4, help="GFF files processed at once in batch mode")
    ap.add_argument("-o", "--output", required=True, help="output folder")
    ap.add_argument("--force", action="store_true", help="overwrite an earlier run in the output folder")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--version", action="version", version="%s %s" % (h.PIPELINE_NAME, h.__version__))


def write_tsv(path, columns, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for r in rows:
            w.writerow(["NA" if r.get(c) in (None, "") else r.get(c) for c in columns])


def main(argv=None):
    ap = argparse.ArgumentParser(prog=MODULE, description=__doc__.strip().split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    add_args(ap)
    a = ap.parse_args(argv)

    single_given = bool(a.ids or a.gff)
    batch_given = bool(a.ids_dir or a.gff_dir)
    if single_given == batch_given:
        h.die("give --ids and --gff for one species, or --ids-dir and --gff-dir for several, not a mix of both")
    if single_given and not (a.ids and a.gff):
        h.die("single-species mode needs both --ids and --gff")
    if batch_given and not (a.ids_dir and a.gff_dir):
        h.die("batch mode needs both --ids-dir and --gff-dir")

    try:
        h.check_output_folder(a.output, MODULE, MODULE, resume=False, force=a.force)
        gff_map = h.load_species_map(a.species_map) if a.species_map else None
        if single_given:
            label = a.species_name or os.path.splitext(os.path.basename(a.ids))[0]
            tasks = {label: (a.ids, a.gff)}
        else:
            tasks = match_files(a.ids_dir, a.gff_dir, gff_map)
            if not tasks:
                raise h.HepError("no ID file could be matched to a GFF file in %s / %s" % (a.ids_dir, a.gff_dir))
    except h.HepError as e:
        h.die(str(e))

    os.makedirs(a.output, exist_ok=True)
    with open(os.path.join(a.output, h.RUN_MARKER), "w", encoding="utf-8") as fh:
        fh.write("module=%s\nlabel=%s\n" % (MODULE, MODULE))

    def bed_path(sp):
        return os.path.join(a.output, re.sub(r"[^A-Za-z0-9_.\-]", "_", sp) + ".bed")

    all_stats, all_rows, all_not_found = [], [], []
    if len(tasks) == 1 or a.threads <= 1:
        for sp, (idf, gff) in tasks.items():
            st, rows, nf = process_one(sp, idf, gff, bed_path(sp), a.id_column)
            all_stats.append(st); all_rows += rows; all_not_found += nf
            if not a.quiet:
                h.log("%s: %s" % (sp, st["status"]))
    else:
        with ProcessPoolExecutor(max_workers=a.threads) as ex:
            futs = {ex.submit(process_one, sp, idf, gff, bed_path(sp), a.id_column): sp
                   for sp, (idf, gff) in tasks.items()}
            for f in as_completed(futs):
                st, rows, nf = f.result()
                all_stats.append(st); all_rows += rows; all_not_found += nf
                if not a.quiet:
                    h.log("%s: %s" % (st["species"], st["status"]))

    write_tsv(os.path.join(a.output, "master.tsv"), MASTER_COLUMNS, all_rows)
    write_tsv(os.path.join(a.output, "not_found.tsv"), NOT_FOUND_COLUMNS, all_not_found)
    write_tsv(os.path.join(a.output, "stats.tsv"),
             ["species", "n_input_ids", "n_found", "n_cds_rows", "status"], all_stats)

    total_in = sum(s["n_input_ids"] for s in all_stats)
    total_found = sum(s["n_found"] for s in all_stats)
    failed_species = [s["species"] for s in all_stats if s["status"].startswith("error:")]
    h.note_incomplete(len(all_not_found) + len(failed_species))
    if not a.quiet:
        bits = []
        if all_not_found:
            bits.append("%d ID(s) not found: see not_found.tsv" % len(all_not_found))
        if failed_species:
            bits.append("%d species could not be processed at all: see stats.tsv" % len(failed_species))
        h.log("done: %d of %d IDs found across %d species. %s"
             % (total_found, total_in, len(tasks), " ".join(bits) if bits else "all found."))


if __name__ == "__main__":
    main()
    sys.exit(h.exit_status())
