#!/usr/bin/env python3
"""
python3 summarize_nodes_by_species_split.py hyper_summary/nodes_summary_by_species brassicaceae \
      --tree SpeciesTree_rooted_node_labels_oct_16.txt

python3 summarize_nodes_by_species_filtered_three_outputs.py input_folder out_prefix \
      --tree SpeciesTree_rooted_node_labels_oct_16.txt

summarize_nodes_by_species_split.py

Summarize all <species>_nodes_summary_filtered.tsv files in an input folder.

This script generates THREE output files:

1) DNG age-summary table only
   <out_prefix>_DNGs_by_age.tsv

   Rows: species names
   Columns: Total, Species-specific, N21, N20, ..., N4
   Counts only rows where Gene_type == DNG.

2) Per-species totals table
   <out_prefix>_species_totals.tsv

   Columns:
     Species_name
     Total_genes
     DNGs
     PDNGs
     %DNGs
     Species_specific_DNGs
     %Spsp_DNGs

3) Long-format gene type table
   <out_prefix>_all_gene_types.tsv

   Columns:
     Species_name
     geneID
     Gene_type

Usage:
  python3 summarize_nodes_by_species_filtered_three_outputs.py input_folder out_prefix \
      --tree SpeciesTree_rooted_node_labels_oct_16.txt

Example:
  python3 summarize_nodes_by_species_filtered_three_outputs.py filtered_nodes all_species_filtered \
      --tree SpeciesTree_rooted_node_labels_oct_16.txt

Input files:
  The script reads only files matching:
      *_nodes_summary_filtered.tsv

Required columns in each input file:
  geneID
  ORF_Node
  Gene_type

Filtering rule:
  Rows with ORF_Node equal to N0, N1, N2, or N3 are excluded from all
  output counts and from the long-format gene list.

Species names are inferred from filenames by removing:
  _nodes_summary_filtered.tsv

Species are ordered by relatedness to Arabidopsis_thaliana based on the tree.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


NODE_RE = re.compile(r"^N\d+$", flags=re.IGNORECASE)
NODE_COLS = ["Species-specific"] + [f"N{i}" for i in range(21, 3, -1)]
EXCLUDED_ORF_NODES = {"N0", "N1", "N2", "N3"}
AGE_HEADER = ["Species_name", "Total"] + NODE_COLS
TOTALS_HEADER = [
    "Species_name",
    "Total_genes",
    "DNGs",
    "PDNGs",
    "%DNGs",
    "Species_specific_DNGs",
    "%Spsp_DNGs",
]
LONG_HEADER = ["Species_name", "geneID", "Gene_type", "ORF_Node"]


# -----------------------------
# Species aliases
# -----------------------------
def canonicalize_species(s: str) -> str:
    aliases = {
        "Athaliana": "Arabidopsis_thaliana",
        "Alyrata": "Arabidopsis_lyrata",
        "Ahalleri": "Arabidopsis_halleri",
        "Bstricta": "Boechera_stricta",
        "Cgrandiflora": "Capsella_grandiflora",
        "Crubella": "Capsella_rubella",
        "Cneglecta": "Camelina_neglecta",
        "Csativa": "Camelina_sativa",
        "Camaraamara": "Cardamine_amara_amara",
        "Aalpina": "Arabis_alpina",
        "Anemorensis": "Arabis_nemorensis",
        "Boleracea": "Brassica_oleracea",
        "BrapaFPsc": "Brassica_rapa_FPsc",
        "Esalsugineum": "Eutrema_salsugineum",
        "Iindigotica": "Isatis_indigotica",
        "Merraticum": "Microthlaspi_erraticum",
        "Rsativus": "Raphanus_sativus",
        "Sparvula": "Schrenkiella_parvula",
        "Salba": "Sinapis_alba",
        "TarvensevarMN106": "Thlaspi_arvense_MN106",
        "Aarabicum": "Aethionema_arabicum",
        "Thassleriana": "Tarenaya_hassleriana",
        "Tcacao": "Theobroma_cacao",
        # lowercase short prefixes
        "athaliana": "Arabidopsis_thaliana",
        "alyrata": "Arabidopsis_lyrata",
        "ahalleri": "Arabidopsis_halleri",
        "bstricta": "Boechera_stricta",
        "cgrandiflora": "Capsella_grandiflora",
        "crubella": "Capsella_rubella",
        "cneglecta": "Camelina_neglecta",
        "csativa": "Camelina_sativa",
        "camaraamara": "Cardamine_amara_amara",
        "aalpina": "Arabis_alpina",
        "anemorensis": "Arabis_nemorensis",
        "boleracea": "Brassica_oleracea",
        "brapa": "Brassica_rapa_FPsc",
        "brapafpsc": "Brassica_rapa_FPsc",
        "esalsugineum": "Eutrema_salsugineum",
        "iindigotica": "Isatis_indigotica",
        "merraticum": "Microthlaspi_erraticum",
        "rsativus": "Raphanus_sativus",
        "sparvula": "Schrenkiella_parvula",
        "salba": "Sinapis_alba",
        "tarvensevarmn106": "Thlaspi_arvense_MN106",
        "aarabicum": "Aethionema_arabicum",
        "thassleriana": "Tarenaya_hassleriana",
        "tcacao": "Theobroma_cacao",
    }
    return aliases.get(s, aliases.get(s.lower(), s))


# ============================================================
# Tiny Newick parser
# ============================================================
@dataclass
class Node:
    name: str
    children: List["Node"]
    parent: Optional["Node"] = None

    def is_leaf(self) -> bool:
        return not self.children


def _tokenize_newick(s: str) -> List[str]:
    tokens: List[str] = []
    buf: List[str] = []
    for ch in s.strip():
        if ch in "(),;":
            if buf:
                tokens.append("".join(buf).strip())
                buf = []
            tokens.append(ch)
        else:
            buf.append(ch)
    if buf:
        tokens.append("".join(buf).strip())
    return [t for t in tokens if t != ""]


def _strip_branch_len(label: str) -> str:
    if ":" in label:
        label = label.split(":", 1)[0]
    return label.strip()


def parse_newick(path: Path) -> Node:
    txt = path.read_text(errors="replace").strip()
    txt = re.sub(r"\s+", "", txt)
    tokens = _tokenize_newick(txt)

    stack: List[Optional[Node]] = []
    current: Optional[Node] = None
    last_closed: Optional[Node] = None
    root_candidate: Optional[Node] = None

    i = 0
    while i < len(tokens):
        t = tokens[i]

        if t == "(":
            n = Node(name="", children=[], parent=current)
            if current is not None:
                current.children.append(n)
            else:
                root_candidate = n
            stack.append(current)
            current = n
            i += 1

        elif t == ",":
            i += 1

        elif t == ")":
            last_closed = current
            parent = stack.pop()
            i += 1
            if i < len(tokens) and tokens[i] not in [",", ")", "(", ";"]:
                if last_closed is not None:
                    last_closed.name = _strip_branch_len(tokens[i])
                i += 1
            current = parent

        elif t == ";":
            i += 1

        else:
            lab = _strip_branch_len(t)
            leaf = Node(name=lab, children=[], parent=current)
            if current is not None:
                current.children.append(leaf)
            else:
                root_candidate = leaf
            i += 1

    root = root_candidate if root_candidate is not None else last_closed
    if root is None:
        raise ValueError("Failed to parse Newick tree.")
    while root.parent is not None:
        root = root.parent
    return root


def collect_leaves(root: Node) -> List[Node]:
    leaves: List[Node] = []

    def dfs(n: Node) -> None:
        if n.is_leaf():
            leaves.append(n)
        else:
            for c in n.children:
                dfs(c)

    dfs(root)
    return leaves


def leaf_depth(leaf: Node) -> int:
    d = 0
    cur = leaf.parent
    while cur is not None:
        d += 1
        cur = cur.parent
    return d


def path_to_root(leaf: Node) -> List[Node]:
    path = [leaf]
    cur = leaf.parent
    while cur is not None:
        path.append(cur)
        cur = cur.parent
    return path


def mrca_depth(a: Node, b: Node) -> int:
    a_path = path_to_root(a)
    b_ids = {id(n): n for n in path_to_root(b)}
    best_depth = -1
    for n in a_path:
        if id(n) in b_ids:
            best_depth = max(best_depth, leaf_depth_like_node(n))
    return best_depth


def leaf_depth_like_node(n: Node) -> int:
    d = 0
    cur = n.parent
    while cur is not None:
        d += 1
        cur = cur.parent
    return d


def species_order_by_relatedness(root: Node, reference_species: str) -> List[str]:
    """
    Order tips from closest to reference_species to most distant.
    Uses deepest MRCA with reference as primary key; input tree traversal order as tie-breaker.
    """
    leaves = collect_leaves(root)
    leaf_by_name = {canonicalize_species(l.name): l for l in leaves}
    ref = canonicalize_species(reference_species)
    if ref not in leaf_by_name:
        raise ValueError(f"Reference species not found in tree tips: {ref}")

    ref_leaf = leaf_by_name[ref]
    ordered: List[Tuple[int, int, str]] = []
    for order_i, leaf in enumerate(leaves):
        sp = canonicalize_species(leaf.name)
        if sp == ref:
            score = 10**9
        else:
            score = mrca_depth(ref_leaf, leaf)
        ordered.append((-score, order_i, sp))

    ordered.sort()
    return [sp for _, _, sp in ordered]


# ============================================================
# Summary helpers
# ============================================================
def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def species_from_nodes_summary_filename(path: Path) -> str:
    name = path.name
    suffix = "_nodes_summary_filtered.tsv"
    if name.endswith(suffix):
        return canonicalize_species(name[: -len(suffix)])
    return canonicalize_species(path.stem)


def normalize_orf_node_for_count(orf_node: str, species_name: str) -> str:
    x = (orf_node or "").strip()
    if not x or x.upper() == "NA":
        return "NA"
    if x.endswith("-specific"):
        return "Species-specific"
    if NODE_RE.match(x):
        return x.upper()
    # also accept exact canonical species-specific weirdness
    if x == f"{species_name}-specific":
        return "Species-specific"
    return x


def fmt_pct(num: int, denom: int) -> str:
    if denom == 0:
        return "0.0"
    return f"{(num / denom * 100.0):.1f}"


def init_age_counts() -> Dict[str, int]:
    return {c: 0 for c in NODE_COLS}


def read_species_file(path: Path, species_name: str) -> Tuple[Dict[str, Dict[str, int]], Dict[str, int], List[Dict[str, str]]]:
    """
    Returns:
      age_counts_by_type: {'DNG': counts, 'PDNG': counts}
      totals: {'Total_genes', 'DNGs', 'PDNGs', 'Species_specific_DNGs'}
      long_rows: list of Species_name/geneID/Gene_type rows
    """
    age_counts_by_type = {
        "DNG": init_age_counts(),
        "PDNG": init_age_counts(),
    }
    totals = {
        "Total_genes": 0,
        "DNGs": 0,
        "PDNGs": 0,
        "Species_specific_DNGs": 0,
    }
    long_rows: List[Dict[str, str]] = []

    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            die(f"Could not read header from {path}")
        required = {"geneID", "ORF_Node", "Gene_type"}
        missing = required - set(reader.fieldnames)
        if missing:
            die(f"Missing required columns in {path}: {', '.join(sorted(missing))}")

        for row in reader:
            gid = (row.get("geneID", "") or "").strip()
            gene_type = (row.get("Gene_type", "") or "").strip().upper()
            orf_node = normalize_orf_node_for_count(row.get("ORF_Node", ""), species_name)

            if not gid:
                continue
            if gene_type not in {"DNG", "PDNG"}:
                continue
            if orf_node in EXCLUDED_ORF_NODES:
                continue

            totals["Total_genes"] += 1
            totals[f"{gene_type}s"] += 1

            if gene_type == "DNG" and orf_node == "Species-specific":
                totals["Species_specific_DNGs"] += 1

            if orf_node in age_counts_by_type[gene_type]:
                age_counts_by_type[gene_type][orf_node] += 1

            long_rows.append(
                {
                    "Species_name": species_name,
                    "geneID": gid,
                    "Gene_type": gene_type,
                    "ORF_Node": orf_node,
                }
            )

    return age_counts_by_type, totals, long_rows


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Summarize *_nodes_summary_filtered.tsv files into DNG age, species totals, and long gene-type tables."
    )
    ap.add_argument("input_folder", help="Folder containing *_nodes_summary_filtered.tsv files.")
    ap.add_argument("out_prefix", help="Output prefix. Three TSV files will be written using this prefix.")
    ap.add_argument("--tree", required=True, help="Newick tree file with labeled internal nodes.")
    ap.add_argument(
        "--reference-species",
        default="Arabidopsis_thaliana",
        help="Species used to order tips by relatedness. Default: Arabidopsis_thaliana",
    )
    args = ap.parse_args()

    input_dir = Path(args.input_folder)
    if not input_dir.exists():
        die(f"Input folder not found: {input_dir}")

    tree_path = Path(args.tree)
    if not tree_path.exists():
        die(f"Tree file not found: {tree_path}")

    files = sorted(input_dir.glob("*_nodes_summary_filtered.tsv"))
    if not files:
        die(f"No *_nodes_summary_filtered.tsv files found in {input_dir}")

    root = parse_newick(tree_path)
    tree_species_order = species_order_by_relatedness(root, args.reference_species)
    tree_rank = {sp: i for i, sp in enumerate(tree_species_order)}

    species_to_file: Dict[str, Path] = {}
    for f in files:
        sp = species_from_nodes_summary_filename(f)
        if sp in species_to_file:
            print(f"WARNING: duplicate species file for {sp}; keeping first: {species_to_file[sp]}", file=sys.stderr)
            print(f"WARNING: ignored duplicate: {f}", file=sys.stderr)
            continue
        species_to_file[sp] = f

    species_order = sorted(species_to_file.keys(), key=lambda sp: (tree_rank.get(sp, 10**9), sp))

    out_dng_age = Path(f"{args.out_prefix}_DNGs_by_age.tsv")
    out_totals = Path(f"{args.out_prefix}_species_totals.tsv")
    out_long = Path(f"{args.out_prefix}_all_gene_types.tsv")

    all_dng_age_rows: List[Dict[str, str]] = []
    all_totals_rows: List[Dict[str, str]] = []
    all_long_rows: List[Dict[str, str]] = []

    for sp in species_order:
        path = species_to_file[sp]
        age_counts_by_type, totals, long_rows = read_species_file(path, sp)

        dng_counts = age_counts_by_type["DNG"]
        dng_age_row: Dict[str, str] = {
            "Species_name": sp,
            "Total": str(totals["DNGs"]),
        }
        for col in NODE_COLS:
            dng_age_row[col] = str(dng_counts.get(col, 0))
        all_dng_age_rows.append(dng_age_row)

        total_genes = totals["Total_genes"]
        dngs = totals["DNGs"]
        pdngs = totals["PDNGs"]
        spsp_dngs = totals["Species_specific_DNGs"]
        all_totals_rows.append(
            {
                "Species_name": sp,
                "Total_genes": str(total_genes),
                "DNGs": str(dngs),
                "PDNGs": str(pdngs),
                "%DNGs": fmt_pct(dngs, total_genes),
                "Species_specific_DNGs": str(spsp_dngs),
                "%Spsp_DNGs": fmt_pct(spsp_dngs, dngs),
            }
        )

        all_long_rows.extend(long_rows)

    with out_dng_age.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=AGE_HEADER)
        writer.writeheader()
        writer.writerows(all_dng_age_rows)

    with out_totals.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=TOTALS_HEADER)
        writer.writeheader()
        writer.writerows(all_totals_rows)

    with out_long.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=LONG_HEADER)
        writer.writeheader()
        writer.writerows(all_long_rows)

    print(f"Wrote: {out_dng_age}", file=sys.stderr)
    print(f"Wrote: {out_totals}", file=sys.stderr)
    print(f"Wrote: {out_long}", file=sys.stderr)


if __name__ == "__main__":
    main()
