"""
step1_coordinates.py

Stage 5 step 1. One table of where every classified locus lies in every genome in which its ORF is found.
For each focal species it joins three stage 3-4 outputs:
  - the stage 3 step 17 coordinate table (<prefix>_*_all_species_summary.tsv: file_name, species_name,
    chromosome, start, end), one row per genome in which the locus alignment has a block;
  - the stage 4 classification (<prefix>_nodes_summary_filtered.tsv: geneID, ORF_Node, Synteny_Node,
    Enablers_Nodes, Gene_type); rows of loci absent from it, or without an ORF node, are dropped;
  - the stage 4 tBLASTn summary (<prefix>_cdngs_tblastn_*_summary.tsv: one 0/1 column per species); a row is
    kept only when the candidate protein is found by tBLASTn in that genome.
This is stage 4 step 9 (nodes_by_coordinates.py) without its Brassicaceae species lists: species come from
tree_config.json, and the gene ID is recovered from file_name by matching the IDs of the classification table
rather than by pattern, so IDs containing underscores (XP_010545484.1) and chromosome names containing
underscores (scaffold_1) are both handled. A trailing ".p" is removed from gene IDs, as stage 4 does. The run stops when
fewer than 95% of a species' classified loci are named in its coordinate table (--min-coord-match), the sign that the
two files write gene IDs differently.

Usage:
  python step1_coordinates.py --tree-config tree_config.json --coords DIR --nodes DIR --tblastn DIR --out coords.tsv
      [--coords-glob '*_all_species_summary.tsv'] [--nodes-glob '*_nodes_summary_filtered.tsv']
      [--tblastn-glob '*_cdngs_tblastn_*_summary.tsv']
  The species of each file is the longest leading run of its name parts that is a species or alias of tree_config.
Output: coords.tsv (geneID, focal_species, file_name, species_name, chromosome, start, end, ORF_Node,
        Synteny_Node, Enablers_Nodes, Gene_type) and a per-species row count on stderr
"""
import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s5_common import Tree, die, log, need, read_tsv, strip_p, write_tsv  # noqa: E402

OUT_COLS = ["geneID", "focal_species", "file_name", "species_name", "chromosome", "start", "end",
            "ORF_Node", "Synteny_Node", "Enablers_Nodes", "Gene_type"]
NA_VALUES = {"", "na", "Na", "NA", "NaN", "nan", "None"}


def species_of_file(tree, fname, suffix_re):
    """Focal species of a file: the longest leading run of underscore-separated name parts that is a
    species name or alias in tree_config (Brapa_FPsc_DNG_... -> Brassica_rapa_FPsc)."""
    stem = re.split(suffix_re, os.path.basename(fname))[0]
    parts = stem.split("_")
    for k in range(len(parts), 0, -1):
        s = tree.canon("_".join(parts[:k]))
        if s in tree.species:
            return s
    return None


def files_by_species(tree, folder, pattern, suffix_re, what):
    out = {}
    for f in sorted(glob.glob(os.path.join(folder, pattern))):
        if os.path.basename(f).startswith("._"):
            continue
        s = species_of_file(tree, f, suffix_re)
        if s is None:
            die(f"{what}: cannot tell the species of {f} from its name (add an alias to tree_config.json)")
        if s in out:
            die(f"{what}: two files for {s}: {out[s]} and {f}")
        out[s] = f
    return out


def gene_from_file_name(fn, ids, cache):
    """The geneID inside an alignment file name (<Species>_<geneID>_<chrom>_<start>_<end>.aln): the longest run of
    underscore-separated parts that is a known geneID of the focal species."""
    if fn in cache:
        return cache[fn]
    base = re.sub(r"\.(aln|fa|fasta|fna)$", "", os.path.basename(fn))
    parts = base.split("_")
    best = None
    for i in range(len(parts)):
        for j in range(len(parts), i, -1):
            g = "_".join(parts[i:j])
            if g in ids and (best is None or len(g) > len(best)):
                best = g
                break
    cache[fn] = best
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree-config", required=True)
    ap.add_argument("--coords", required=True, help="stage 3 step 17 tables, *_all_species_summary.tsv")
    ap.add_argument("--nodes", required=True, help="stage 4 *_nodes_summary_filtered.tsv")
    ap.add_argument("--tblastn", required=True, help="stage 4 *_cdngs_tblastn_*_summary.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--coords-glob", default="*_all_species_summary.tsv")
    ap.add_argument("--nodes-glob", default="*_nodes_summary_filtered.tsv")
    ap.add_argument("--tblastn-glob", default="*_cdngs_tblastn_*_summary.tsv")
    ap.add_argument("--min-coord-match", type=float, default=0.95,
                    help="stop when fewer than this share of a species' classified loci is named in its coordinate table")
    a = ap.parse_args()
    need(a.tree_config, a.coords, a.nodes, a.tblastn)
    if os.path.exists(a.out):
        die(f"output exists: {a.out}")
    tree = Tree(a.tree_config)
    cf = files_by_species(tree, a.coords, a.coords_glob, r"_all_species_summary\.tsv$", "coords")
    nf = files_by_species(tree, a.nodes, a.nodes_glob, r"_nodes_summary_filtered\.tsv$", "nodes")
    tf = files_by_species(tree, a.tblastn, a.tblastn_glob, r"_cdngs_tblastn_", "tblastn")
    missing = {s: [w for w, d in (("coords", cf), ("nodes", nf), ("tblastn", tf)) if s not in d]
               for s in set(cf) | set(nf) | set(tf)}
    missing = {s: m for s, m in missing.items() if m}
    if missing:
        die(f"species without all three inputs: {missing}")
    out, report = [], []
    for focal in [s for s in tree.species if s in cf]:
        nodes = {}
        for r in read_tsv(nf[focal]):
            g = strip_p(r["geneID"])
            if g in nodes:
                die(f"{nf[focal]}: geneID {g} twice")
            nodes[g] = r
        tb_rows = read_tsv(tf[focal])
        if not tb_rows:
            die(f"{tf[focal]}: no rows")
        tb = {}
        for r in tb_rows:
            g = strip_p(r["geneID"] or "")
            if g in tb:
                die(f"{tf[focal]}: geneID {g} twice")
            tb[g] = r
        unknown_cols = [c for c in tb_rows[0] if c not in ("geneID", "Node_tblastn")
                        and tree.canon(c) not in tree.species]
        if unknown_cols:
            die(f"{tf[focal]}: columns that are not species of tree_config: {unknown_cols}")
        cache = {}
        rows = read_tsv(cf[focal])
        n_in = len(rows)
        n_nogene = n_noorf = n_notb = 0
        kept = []
        for r in rows:
            g = gene_from_file_name(r["file_name"], nodes, cache)
            if g is None:
                n_nogene += 1
                continue
            nr = nodes[g]
            if nr.get("ORF_Node", "") in NA_VALUES:
                n_noorf += 1
                continue
            sp = tree.canon(r["species_name"])
            t = tb.get(g)
            col = next((c for c in (t or {}) if tree.canon(c) == sp and c not in ("geneID", "Node_tblastn")), None)
            try:
                present = t is not None and col is not None and float(t[col]) > 0
            except (ValueError, TypeError):
                present = False
            if not present:
                n_notb += 1
                continue
            enab = nr.get("Enablers_Nodes", "NA")
            kept.append([g, focal, r["file_name"], r["species_name"], r["chromosome"], r["start"], r["end"],
                         nr["ORF_Node"], nr.get("Synteny_Node", "NA"), "NA" if enab in NA_VALUES else enab,
                         nr["Gene_type"]])
        assert n_in == len(kept) + n_nogene + n_noorf + n_notb
        named = {v for v in cache.values() if v is not None}
        if len(named) < a.min_coord_match * len(nodes):
            die(f"{focal}: only {len(named)} of {len(nodes)} classified loci are named in {cf[focal]}; the gene IDs of "
                f"the classification (e.g. {sorted(set(nodes) - named)[:3]}) and of the alignment file names (e.g. "
                f"{[r['file_name'] for r in rows[:2]]}) probably differ in form")
        loci_with_self = len({k[0] for k in kept if tree.canon(k[3]) == focal})
        report.append((focal, n_in, n_nogene, n_noorf, n_notb, len(kept), len(nodes), loci_with_self))
        out.extend(kept)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    write_tsv(a.out, OUT_COLS, out)
    log("focal_species\trows_in\tno_classified_locus\tno_ORF_node\tnot_found_by_tBLASTn\trows_kept"
        "\tloci_classified\tloci_with_own_coordinates")
    for r in report:
        log("\t".join(map(str, r)))
    log(f"wrote {a.out}: {len(out):,} rows")


if __name__ == "__main__":
    main()
