#!/usr/bin/env python3
"""
stage3_to_stage4.py

Copies the output of one stage 3 run (one focal species) into the stage 4 work folder, in the layout run_pipeline.sh
reads from step 2 on. Run it once per focal species, into the same WORKDIR, then run stage 4 from step 2 with ONLY
listing the species copied in.

  WORKDIR/brassicaceae_23species_LINGUA_results/<dir>_stage3_alignments/*.aln   stage 3 step9/*.aln
  WORKDIR/stops_dir/<Abbr>_DNG_45k_cds_stop_codon_summary.tsv                    stage 3 step13/master_stop.tsv
  WORKDIR/indels_dir/<Abbr>_DNG_45k_cds_indels_summary.tsv                       stage 3 step14/master_indels.tsv
  WORKDIR/coordinates_alignments_all_species/<Abbr>_DNG_45k_cds_all_species_summary.tsv
                                                                                 stage 3 step17_all_species.tsv,
                                                                                 first five columns
  WORKDIR/clsgs_45k_1line/<Abbr>.faa                                             the candidate proteins, as
                                                                                 >ID__<Abbr>, ID being the transcript
                                                                                 ID of the alignment (a ".p" suffix
                                                                                 is matched and dropped)

The species is given by its query_prefix in species.tsv (SPECIES_TSV, default data/species.tsv). The names are
derived, not typed, because three steps read them. Step 2 names the gap-free folder after the alignment folder's
prefix, which must give the species' nogaps_dir (<dir> = nogaps_dir without "_stitched_nogaps"). Step 3 names the
query file after the header suffix and step 4 strips that suffix only when it is an alias of the tree configuration,
so <Abbr> is the alias whose lowercase form is the query_prefix. Step 9 reads the focal species from the coordinate
table name, so <Abbr> must also map to the focal species there. The stop and indel tables are found by content.

Refuses, writing nothing, when the stage 3 run has no SUCCEEDED status, when the species is unknown or its names
cannot be derived, when the alignments belong to another species, when a stop, indel or coordinate table is empty or
names a locus with no step 9 alignment, when an alignment has no protein, or when a target exists. A write that
fails part-way removes what it wrote.

Usage: python3 stage3_to_stage4.py WORKDIR STAGE3_DIR QUERY_PREFIX PROTEINS_FASTA [--species-tsv FILE]
Example: python3 stage3_to_stage4.py work stage3_athaliana athaliana athaliana_candidates.faa
"""
import argparse
import csv
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lingua_config as cfg  # noqa: E402
from nodes_by_coordinates import canonicalize_species as step9_species  # noqa: E402

csv.field_size_limit(10 ** 9)
NOGAPS_SUFFIX = "_stitched_nogaps"


def die(msg):
    sys.exit(f"ERROR: {msg}")


def need(path):
    if not os.path.exists(path):
        die(f"missing {path}")
    return path


def read_fasta(path):
    recs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name, seq = line[1:].split()[0].split("__", 1)[0], []
            elif line.strip():
                seq.append(line.strip())
    if name is not None:
        recs.append((name, "".join(seq)))
    return recs


def bare(gene_id):
    return gene_id[:-2] if gene_id.endswith(".p") else gene_id


def names_for(species_tsv, query_prefix):
    with open(need(species_tsv), newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("query_prefix") == query_prefix]
    if len(rows) != 1:
        die(f"{len(rows)} rows of {species_tsv} have query_prefix {query_prefix} (need exactly one)")
    r = rows[0]
    focal = cfg.canonicalize(r["focal_species"])
    nog = r["nogaps_dir"]
    aln_prefix = nog[:-len(NOGAPS_SUFFIX)] if nog.endswith(NOGAPS_SUFFIX) else ""
    if not aln_prefix or "_" in aln_prefix:
        die(f"nogaps_dir {nog} of {query_prefix} is not <name>{NOGAPS_SUFFIX} with no '_' in <name>")
    abbr = [k for k, v in cfg.ALIASES.items()
            if v == focal and k.lower() == query_prefix and step9_species(k) == focal]
    if not abbr:
        die(f"no alias of {focal} in the tree configuration is spelled {query_prefix} (any case) and maps to "
            f"{focal} in step 9; add one to the configuration's aliases")
    abbr.sort(key=lambda k: (k.islower(), k))
    return focal, aln_prefix, abbr[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workdir")
    ap.add_argument("stage3_dir")
    ap.add_argument("query_prefix")
    ap.add_argument("proteins")
    ap.add_argument("--species-tsv", default=os.environ.get("SPECIES_TSV") or
                    os.path.join(os.path.dirname(HERE), "data", "species.tsv"))
    a = ap.parse_args()
    qp = a.query_prefix.lower()
    focal, aln_prefix, abbr = names_for(a.species_tsv, qp)

    s3 = a.stage3_dir
    status = open(need(os.path.join(s3, "STATUS"))).read().split()
    if not status or status[0] != "SUCCEEDED":
        die(f"{s3}/STATUS does not say SUCCEEDED ({' '.join(status[:2]) or 'empty'})")
    step9 = need(os.path.join(s3, "step9"))
    stops = need(os.path.join(s3, "step13", "master_stop.tsv"))
    indels = need(os.path.join(s3, "step14", "master_indels.tsv"))
    coords = need(os.path.join(s3, "step17_all_species.tsv"))
    report = need(os.path.join(step9, "step9_trimming_report.tsv"))
    need(a.proteins)

    alns = sorted(f for f in os.listdir(step9) if f.endswith(".aln"))
    if not alns:
        die(f"no .aln files in {step9}")
    names = {f[:-4] for f in alns}
    off = [f for f in alns if not f.startswith(focal + "_")]
    if off:
        die(f"{len(off)} alignments do not start with {focal}_ (the species of {qp}), e.g. {off[0]}")

    def check_table(path, col):
        with open(path, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if col not in (rd.fieldnames or []):
                die(f"{path} has no {col} column")
            seen = {r[col][:-4] if r[col].endswith(".aln") else r[col] for r in rd}
        if not seen:
            die(f"{path} has no rows; stage 4 needs at least one row per species in this table")
        extra = seen - names
        if extra:
            die(f"{len(extra)} loci in {path} have no step 9 alignment, e.g. {sorted(extra)[0]}")
        return len(seen)

    n_stop = check_table(stops, "file_name")
    n_indel = check_table(indels, "file_name")
    n_coord = check_table(coords, "file_name")

    with open(report, newline="") as fh:
        rd = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        tx = {r["file"][:-4]: r["transcript"] for r in rd if r["file"][:-4] in names}
    if set(tx) != names:
        die(f"{len(names - set(tx))} alignments are missing from {report}")
    prots = read_fasta(a.proteins)
    by_id = {}
    for pid, _ in prots:
        if bare(pid) in by_id:
            die(f"protein {pid} appears twice in {a.proteins}")
        by_id[bare(pid)] = pid
    no_prot = sorted(n for n, t in tx.items() if bare(t) not in by_id)
    if no_prot:
        die(f"{len(no_prot)} alignments have no protein in {a.proteins}, e.g. {no_prot[0]}")
    as_aln = {bare(t): t for t in tx.values()}
    no_aln = sorted(set(by_id) - set(as_aln))

    w = a.workdir
    res = os.path.join(w, "brassicaceae_23species_LINGUA_results")
    aln_dir = os.path.join(res, f"{aln_prefix}_stage3_alignments")
    targets = {
        "stops": os.path.join(w, "stops_dir", f"{abbr}_DNG_45k_cds_stop_codon_summary.tsv"),
        "indels": os.path.join(w, "indels_dir", f"{abbr}_DNG_45k_cds_indels_summary.tsv"),
        "coords": os.path.join(w, "coordinates_alignments_all_species", f"{abbr}_DNG_45k_cds_all_species_summary.tsv"),
        "prot": os.path.join(w, "clsgs_45k_1line", f"{abbr}.faa"),
    }
    for p in [aln_dir] + list(targets.values()):
        if os.path.exists(p):
            die(f"{p} exists; use a new WORKDIR or move it away")
    if os.path.isdir(res):
        for d in os.listdir(res):
            if d.split("_", 1)[0].lower() == aln_prefix:
                die(f"{res}/{d} already holds alignments for {qp} (step 2 would merge them)")

    written = []
    try:
        os.makedirs(res, exist_ok=True)
        os.makedirs(aln_dir)
        written.append(aln_dir)
        for d in ("stops_dir", "indels_dir", "coordinates_alignments_all_species", "clsgs_45k_1line"):
            os.makedirs(os.path.join(w, d), exist_ok=True)
        for f in alns:
            shutil.copy2(os.path.join(step9, f), os.path.join(aln_dir, f))
        written.append(targets["stops"])
        shutil.copy2(stops, targets["stops"])
        written.append(targets["indels"])
        shutil.copy2(indels, targets["indels"])
        written.append(targets["coords"])
        with open(coords, newline="") as fh, open(targets["coords"], "w") as o:
            for line in fh:
                o.write("\t".join(line.rstrip("\n").split("\t")[:5]) + "\n")
        written.append(targets["prot"])
        with open(targets["prot"], "w") as o:
            for pid, seq in prots:
                o.write(f">{as_aln.get(bare(pid), pid)}__{abbr}\n{seq}\n")
    except BaseException:
        for p in reversed(written):
            if os.path.isdir(p):
                shutil.rmtree(p, ignore_errors=True)
            elif os.path.exists(p):
                os.remove(p)
        print("ERROR: the copy failed; what this run wrote was removed", file=sys.stderr)
        raise

    print(f"{focal} ({qp}): {len(alns)} alignments; loci in the stop table {n_stop}, indel table {n_indel}, "
          f"coordinate table {n_coord}; proteins {len(prots)} ({len(no_aln)} without an alignment)", file=sys.stderr)
    if no_aln:
        print(f"WARNING: {len(no_aln)} proteins have no step 9 alignment (stage 4 searches them anyway), "
              f"e.g. {no_aln[0]}", file=sys.stderr)
    print(f'written into {w}; after the last species: ONLY="{qp} ..." bash 4_enablers/run_pipeline.sh {w} 2 9 '
          f"(ONLY lists every species copied in; leave it out only when all species of species.tsv are there)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
