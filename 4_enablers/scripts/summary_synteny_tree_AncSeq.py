#!/usr/bin/env python3

"""
python3 summary_synteny_tree_AncSeq.py \
  --batch brassicaceae_23species_LINGUA_nogaps \
  --batch-outdir synteny_40perc \
  --tree SpeciesTree_rooted_node_labels_oct_16.txt \
  --synteny 40


python3 summary_synteny_tree_AncSeq.py brassicaceae_23species_LINGUA_nogaps/athaliana_stitched_nogaps synteny_40perc/athaliana_stitched_nogaps-synteny-40-AS.tsv synteny_40perc/athaliana_stitched_nogaps-synteny-40_summary-AS.tsv synteny_40perc/athaliana_stitched_nogaps-synteny-40_AncSeq.tsv --tree SpeciesTree_rooted_node_labels_oct_16.txt --focal-species Arabidopsis_thaliana --synteny 40


BATCH USAGE
-----------
Run all species folders in an input parent directory. The script will search for
subfolders named <species>_stitched_nogaps, infer the focal species from <species>,
and write all three outputs for each folder into --batch-outdir.

python3 summary_synteny_tree_AncSeq.py --batch brassicaceae_23species_LINGUA_nogaps --batch-outdir synteny_40perc --tree SpeciesTree_rooted_node_labels_oct_16.txt --synteny 40

Batch output names:
  <species>_stitched_nogaps-synteny-40.tsv
  <species>_stitched_nogaps-synteny-40_summary.tsv
  <species>_stitched_nogaps-synteny-40_AncSeq.tsv


Tree-aware synteny summary (DNA-only)

This is an updated version of summary_synteny_v1.py that derives node/clade logic
from a Newick tree with labeled internal nodes (e.g., N0, N1, ...).

USAGE
-----
python3 summary_synteny_v2.py <input_folder> <out_per_seq.tsv> <out_by_gene.tsv> \
  --synteny 40 --tree SpeciesTree_rooted_node_labels_oct_16.txt --focal-species Arabidopsis_thaliana

Notes
-----
- Species presence/absence for the gene-level summary is determined from per-sequence
  %Synteny >= --synteny for TIP species only (the 23 species below).
- Ancestor/ancestral sequences in the FASTA headers (e.g., "N9.something") are written
  to the per-seq file with Species="N9" etc., but they do not contribute to the
  per-gene species presence table.
- Node_Synteny-XX is computed using a tree-based rule:
    * requires focal_species == 1
    * chooses the YOUNGEST ancestor node on focal's lineage whose clade contains all positives
    * returns "<focal>-specific" if only focal is positive
    * returns "NA" if focal is not positive

"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

FASTA_EXTS = {".fna", ".fa", ".fasta", ".fas"}

# -----------------------------
# Species column order (tips)
# -----------------------------
SPECIES_COLUMNS: List[str] = [
    "Arabidopsis_thaliana",
    "Arabidopsis_lyrata",
    "Arabidopsis_halleri",
    "Boechera_stricta",
    "Capsella_grandiflora",
    "Capsella_rubella",
    "Camelina_neglecta",
    "Camelina_sativa",
    "Cardamine_amara_amara",
    "Arabis_alpina",
    "Arabis_nemorensis",
    "Brassica_oleracea",
    "Brassica_rapa_FPsc",
    "Eutrema_salsugineum",
    "Isatis_indigotica",
    "Microthlaspi_erraticum",
    "Raphanus_sativus",
    "Schrenkiella_parvula",
    "Sinapis_alba",
    "Thlaspi_arvense_MN106",
    "Aethionema_arabicum",
    "Tarenaya_hassleriana",
    "Theobroma_cacao",
]

TIP_SET: Set[str] = set(SPECIES_COLUMNS)
NODE_RE = re.compile(r"^N\d+$", flags=re.IGNORECASE)


# ============================================================
# FASTA reader
# ============================================================
def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header: Optional[str] = None
    seq_chunks: List[str] = []
    with path.open("r", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                if header is not None:
                    seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


# ============================================================
# Helpers
# ============================================================
def species_from_header(h: str) -> str:
    """
    Species label used in per-seq output.
    Rule: take token before first whitespace, then prefix before first '.'
      - Arabidopsis_thaliana.Chr1 -> Arabidopsis_thaliana
      - N9.N9refChr165 -> N9
    """
    token = h.split()[0]
    return token.split(".", 1)[0]


def parse_geneid_from_filename(fname: str) -> str:
    stem = Path(fname).stem

    # Match the longest species name prefix from SPECIES_COLUMNS
    species_prefix = None
    for sp in sorted(SPECIES_COLUMNS, key=len, reverse=True):
        prefix = sp + "_"
        if stem.startswith(prefix):
            species_prefix = sp
            rest = stem[len(prefix):]
            break

    if species_prefix is None:
        return stem

    parts = rest.split("_")
    if not parts:
        return stem

    # Case 1: protein accessions like XP_010518862.1
    if len(parts) >= 2 and parts[0] in {"XP", "NP", "YP", "WP"}:
        return parts[0] + "_" + parts[1]

    # Case 2: geneIDs like Aa_G1910.h1.t1
    if len(parts) >= 2 and re.fullmatch(r"[A-Za-z]{1,4}", parts[0]) and re.match(r"^G", parts[1]):
        return parts[0] + "_" + parts[1]

    # Default case: geneID is the first token after species name
    return parts[0]


def clean_dna(seq: str) -> str:
    return re.sub(r"\s+", "", seq).upper()


def ungapped_length(seq: str) -> int:
    return len(seq.replace("-", ""))


def starts_with_atg(seq: str) -> int:
    s = seq.replace("-", "")
    return 1 if len(s) >= 3 and s[:3] == "ATG" else 0


def find_fasta_files(root: Path) -> List[Path]:
    files: List[Path] = []
    for p in root.rglob("*"):
        if p.is_file() and p.suffix.lower() in FASTA_EXTS and not p.name.startswith("."):
            files.append(p)
    return sorted(files)


def format_threshold_for_colname(x: float) -> str:
    if float(x).is_integer():
        return str(int(x))
    return f"{x}".rstrip("0").rstrip(".")


# ============================================================
# Newick parser (tiny; supports node labels + branch lengths)
# ============================================================
@dataclass
class Node:
    name: str
    children: List["Node"]
    parent: Optional["Node"] = None

    def is_leaf(self) -> bool:
        return not self.children


def _tokenize_newick(s: str) -> List[str]:
    # Keep delimiters ;,(,)
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
    # remove :0.123 if present
    if ":" in label:
        label = label.split(":", 1)[0]
    return label.strip()


def parse_newick(path: Path) -> Node:
    txt = path.read_text(errors="replace").strip()
    # allow multi-line
    txt = re.sub(r"\s+", "", txt)
    tokens = _tokenize_newick(txt)
    stack: List[Node] = []
    current: Optional[Node] = None

    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t == "(":
            # start a new internal node (name may appear after ')')
            n = Node(name="", children=[], parent=None)
            if current is not None:
                n.parent = current
                current.children.append(n)
            stack.append(current)  # previous current
            current = n
            i += 1
        elif t == ",":
            i += 1
        elif t == ")":
            # finish current internal node; next token might be its name/branchlen
            parent = current.parent if current else None
            prev = stack.pop()  # restore previous current
            # assign internal node name if next token is a label (not delimiter)
            i += 1
            if i < len(tokens) and tokens[i] not in [",", ")", "(", ";"]:
                current.name = _strip_branch_len(tokens[i])
                i += 1
            # move up
            current = prev if prev is not None else current
        elif t == ";":
            i += 1
        else:
            # leaf label
            lab = _strip_branch_len(t)
            leaf = Node(name=lab, children=[], parent=current)
            if current is None:
                # single leaf tree
                current = leaf
            else:
                current.children.append(leaf)
            i += 1

    # Find root: walk up from any node
    root = current
    while root and root.parent:
        root = root.parent
    if root is None:
        raise ValueError("Failed to parse Newick: root is None")
    return root


def derive_node_to_tips(root: Node) -> Dict[str, Set[str]]:
    node_to_tips: Dict[str, Set[str]] = {}

    def dfs(n: Node) -> Set[str]:
        if n.is_leaf():
            return {n.name}
        tips: Set[str] = set()
        for c in n.children:
            tips |= dfs(c)
        if NODE_RE.match(n.name or ""):
            node_to_tips[n.name] = tips.copy()
        return tips

    dfs(root)
    return node_to_tips


def derive_species_to_ancestors(root: Node) -> Dict[str, List[str]]:
    # build parent pointers already exist; build mapping leaf -> ancestors (youngest->oldest)
    leaves: List[Node] = []

    def collect(n: Node):
        if n.is_leaf():
            leaves.append(n)
        else:
            for c in n.children:
                collect(c)

    collect(root)

    out: Dict[str, List[str]] = {}
    for leaf in leaves:
        anc: List[str] = []
        cur = leaf.parent
        while cur is not None:
            if NODE_RE.match(cur.name or ""):
                anc.append(cur.name)
            cur = cur.parent
        out[leaf.name] = anc  # already youngest->oldest as we go up
    return out


def canonicalize_species(s: str) -> str:
    # normalize common abbreviations used in your pipeline
    # (keep minimal; can expand later)
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
    }
    return aliases.get(s, s)


def classify_node_tree_based(
    pres: Dict[str, int],
    focal_species: str,
    node_to_tips: Dict[str, Set[str]],
    species_to_ancestors: Dict[str, List[str]],
) -> str:
    focal = canonicalize_species(focal_species)
    if focal not in pres or pres.get(focal, 0) == 0:
        return "NA"

    positives: Set[str] = {sp for sp, v in pres.items() if v == 1}
    if positives == {focal}:
        return f"{focal}-specific"

    # lineage nodes youngest->oldest
    anc = species_to_ancestors.get(focal, [])
    # choose youngest node N on lineage such that positives ⊆ clade(N)
    for nlab in anc:
        clade = node_to_tips.get(nlab, set())
        # Use canonical tip names for comparison
        if positives.issubset(clade):
            return nlab

    # If not found, root should contain all; fallback:
    return anc[-1] if anc else "NA"


# ============================================================
# Main
# ============================================================
def main() -> int:
    ap = argparse.ArgumentParser(description="DNA-only synteny summarizer (per-seq + per-gene TSV outputs), tree-aware.")

    # Single-folder mode, preserving the original command-line interface:
    #   input_folder out_per_seq_tsv out_by_gene_tsv out_by_ancseq_tsv --focal-species ...
    ap.add_argument("input_folder", nargs="?")
    ap.add_argument("out_per_seq_tsv", nargs="?")
    ap.add_argument("out_by_gene_tsv", nargs="?")
    ap.add_argument("out_by_ancseq_tsv", nargs="?")

    # Batch mode:
    #   --batch <parent_dir> --batch-outdir <outdir>
    ap.add_argument("--batch", help="Parent folder containing subfolders named <species>_stitched_nogaps.")
    ap.add_argument("--batch-outdir", help="Output directory for batch mode.")

    ap.add_argument("--synteny", type=float, default=40.0, help="Threshold for calling 'synteny conserved' (default: 40)")
    ap.add_argument("--tree", required=True, help="Newick file with internal node labels (e.g., N0, N1, ...).")
    ap.add_argument("--focal-species", help="Focal species for single-folder mode (canonical or alias). Not needed in batch mode.")
    ap.add_argument("--debug-nodes", action="store_true", help="Print per-gene debug info about positives/clade checks to stderr.")
    args = ap.parse_args()

    tree_path = Path(args.tree)
    if not tree_path.exists():
        print(f"ERROR: tree file not found: {tree_path}", file=sys.stderr)
        return 2

    root = parse_newick(tree_path)
    node_to_tips_raw = derive_node_to_tips(root)
    species_to_ancestors_raw = derive_species_to_ancestors(root)

    # canonicalize tree tips for mapping (tree already should use canonical names)
    node_to_tips: Dict[str, Set[str]] = {k: {canonicalize_species(x) for x in v} for k, v in node_to_tips_raw.items()}
    species_to_ancestors: Dict[str, List[str]] = {canonicalize_species(k): v for k, v in species_to_ancestors_raw.items()}

    thr_label = format_threshold_for_colname(args.synteny)
    node_synteny_col = f"Node_Synteny-{thr_label}"
    node_synteny_ancseq_col = f"Node_Synteny_AncSeq-{thr_label}"
    all_anc_nodes = sorted(node_to_tips.keys(), key=lambda x: int(x[1:]), reverse=True)

    def run_one(input_folder: Path, out_per_seq_tsv: Path, out_by_gene_tsv: Path, out_by_ancseq_tsv: Path, focal_species: str) -> int:
        in_dir = Path(input_folder)
        if not in_dir.exists():
            print(f"ERROR: input folder not found: {in_dir}", file=sys.stderr)
            return 2

        focal = canonicalize_species(focal_species)
        if focal not in species_to_ancestors:
            print(f"ERROR: focal species not found in tree tips: {focal}", file=sys.stderr)
            return 2

        ancestors_young_to_old = species_to_ancestors.get(focal, [])

        # Dynamic node columns for this focal species (youngest->oldest)
        node_cols = [f"{n}_#species_Synt" for n in ancestors_young_to_old]

        # -----------------------------
        # Headers
        # -----------------------------
        per_seq_header = [
            "geneID",
            "file_name",
            "Species",
            "GeneID+Species",
            "AncSeq",
            "top_CDS_L",
            "other_CDS_L",
            "%Synteny",
            "Start",
        ]

        by_gene_header = ["geneID", node_synteny_col] + node_cols + SPECIES_COLUMNS

        fasta_files = find_fasta_files(in_dir)
        if not fasta_files:
            print(f"ERROR: no FASTA files found in {in_dir}.", file=sys.stderr)
            return 2

        out_per_seq_tsv.parent.mkdir(parents=True, exist_ok=True)
        out_by_gene_tsv.parent.mkdir(parents=True, exist_ok=True)
        out_by_ancseq_tsv.parent.mkdir(parents=True, exist_ok=True)

        # gene -> species -> 0/1 pass
        gene_species_pass: Dict[str, Dict[str, int]] = {}

        # gene -> ancestral node -> 0/1 pass for out_by_ancseq.tsv
        gene_ancseq_pass: Dict[str, Dict[str, int]] = {}

        # ============================================================
        # Per-sequence processing
        # ============================================================
        with open(out_per_seq_tsv, "w", newline="") as out1:
            w1 = csv.writer(out1, delimiter="\t")
            w1.writerow(per_seq_header)

            for fp in fasta_files:
                records = list(iter_fasta(fp))
                if not records:
                    continue

                geneID = parse_geneid_from_filename(fp.name)

                top_seq = clean_dna(records[0][1])
                top_len = ungapped_length(top_seq)

                if geneID not in gene_species_pass:
                    gene_species_pass[geneID] = {sp: 0 for sp in SPECIES_COLUMNS}
                    gene_ancseq_pass[geneID] = {n: 0 for n in all_anc_nodes}

                for h, seq_raw in records:
                    seq = clean_dna(seq_raw)
                    sp_raw = species_from_header(h)
                    sp = canonicalize_species(sp_raw)
                    anc = 1 if NODE_RE.match(sp_raw or "") else 0

                    other_len = ungapped_length(seq)
                    pct = (other_len / top_len * 100.0) if top_len > 0 else 0.0
                    start_flag = starts_with_atg(seq)

                    w1.writerow(
                        [
                            geneID,
                            fp.name,
                            sp_raw,  # keep raw label in per-seq
                            f"{geneID}_{sp_raw}",
                            anc,
                            top_len,
                            other_len,
                            f"{pct:.2f}",
                            start_flag,
                        ]
                    )

                    # Only TIP species contribute to presence table
                    if sp in TIP_SET and pct >= args.synteny:
                        gene_species_pass[geneID][sp] = 1

                    # Only ancestral-node rows contribute to out_by_ancseq.tsv
                    if NODE_RE.match(sp_raw or "") and pct >= args.synteny:
                        if sp_raw in gene_ancseq_pass[geneID]:
                            gene_ancseq_pass[geneID][sp_raw] = 1

        # ============================================================
        # Gene summary
        # ============================================================
        with open(out_by_gene_tsv, "w", newline="") as out2:
            w2 = csv.writer(out2, delimiter="\t")
            w2.writerow(by_gene_header)

            for geneID in sorted(gene_species_pass.keys()):
                spflags = gene_species_pass[geneID]

                node_counts: Dict[str, int] = {}
                positives: Set[str] = {sp for sp, v in spflags.items() if v == 1}

                # Count only the species newly added at each successive node on the focal lineage,
                # matching the behavior of summary_synteny_v1.py. In other words, for each node,
                # count positives present in that node clade but absent from the next younger node clade.
                for i, nlab in enumerate(ancestors_young_to_old):
                    clade = node_to_tips.get(nlab, set())
                    if i == 0:
                        incremental = clade - {focal}
                    else:
                        younger_clade = node_to_tips.get(ancestors_young_to_old[i - 1], set())
                        incremental = clade - younger_clade
                    node_counts[nlab] = sum(1 for sp in positives if sp in incremental)

                node_label = classify_node_tree_based(
                    pres=spflags,
                    focal_species=focal,
                    node_to_tips=node_to_tips,
                    species_to_ancestors=species_to_ancestors,
                )

                if args.debug_nodes:
                    # show why it matched
                    print(f"[DEBUG] {geneID}\tfocal={focal}\tpositives={','.join(sorted(positives)) or 'NONE'}\tNode={node_label}", file=sys.stderr)

                row = (
                    [geneID, node_label]
                    + [node_counts.get(nlab, 0) for nlab in ancestors_young_to_old]
                    + [spflags[sp] for sp in SPECIES_COLUMNS]
                )
                w2.writerow(row)

        # ============================================================
        # Ancestral sequence summary
        # ============================================================
        with open(out_by_ancseq_tsv, "w", newline="") as out3:
            w3 = csv.writer(out3, delimiter="\t")
            w3.writerow(["geneID", node_synteny_ancseq_col] + all_anc_nodes)

            for geneID in sorted(gene_ancseq_pass.keys()):
                node_flags = gene_ancseq_pass[geneID]

                oldest_node = f"{focal}-specific"
                # Only ancestors of the focal species can be its oldest syntenic node.
                lineage = set(ancestors_young_to_old)
                for nlab in sorted(all_anc_nodes, key=lambda x: int(x[1:])):  # N0 oldest -> N21 youngest
                    if nlab not in lineage:
                        continue
                    if node_flags.get(nlab, 0) > 0:
                        oldest_node = nlab
                        break

                row = [geneID, oldest_node] + [node_flags.get(nlab, 0) for nlab in all_anc_nodes]
                w3.writerow(row)

        return 0

    # ============================================================
    # Batch mode
    # ============================================================
    if args.batch:
        if args.input_folder or args.out_per_seq_tsv or args.out_by_gene_tsv or args.out_by_ancseq_tsv:
            print("ERROR: in --batch mode, do not provide the four positional single-folder arguments.", file=sys.stderr)
            return 2
        if not args.batch_outdir:
            print("ERROR: --batch-outdir is required with --batch.", file=sys.stderr)
            return 2

        batch_parent = Path(args.batch)
        if not batch_parent.exists() or not batch_parent.is_dir():
            print(f"ERROR: batch input folder not found or not a directory: {batch_parent}", file=sys.stderr)
            return 2

        outdir = Path(args.batch_outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        species_dir_to_focal = {
            "athaliana": "Arabidopsis_thaliana",
            "alyrata": "Arabidopsis_lyrata",
            "ahalleri": "Arabidopsis_halleri",
            "bstricta": "Boechera_stricta",
            "cgrandiflora": "Capsella_grandiflora",
            "crubella": "Capsella_rubella",
            "cneglecta": "Camelina_neglecta",
            "csativa": "Camelina_sativa",
            "camara": "Cardamine_amara_amara",
            "camaraamara": "Cardamine_amara_amara",
            "camara_amara": "Cardamine_amara_amara",
            "aalpina": "Arabis_alpina",
            "anemorensis": "Arabis_nemorensis",
            "boleracea": "Brassica_oleracea",
            "brapa": "Brassica_rapa_FPsc",
            "esalsugineum": "Eutrema_salsugineum",
            "iindigotica": "Isatis_indigotica",
            "merraticum": "Microthlaspi_erraticum",
            "rsativus": "Raphanus_sativus",
            "sparvula": "Schrenkiella_parvula",
            "salba": "Sinapis_alba",
            "tarvense": "Thlaspi_arvense_MN106",
            "tarvensevarmn106": "Thlaspi_arvense_MN106",
            "tarvensevarMN106": "Thlaspi_arvense_MN106",
            "aarabicum": "Aethionema_arabicum",
            "thassleriana": "Tarenaya_hassleriana",
            "tcacao": "Theobroma_cacao",
        }

        subdirs = sorted(p for p in batch_parent.iterdir() if p.is_dir() and p.name.endswith("_stitched_nogaps"))
        if not subdirs:
            # Fallback: batch_parent itself may be a flat FASTA folder (e.g. athaliana_nestedorfs_nogaps).
            # Infer focal species from the directory name prefix using the same species key table.
            focal_species = None
            dir_name_lower = batch_parent.name.lower()
            for key in sorted(species_dir_to_focal.keys(), key=len, reverse=True):
                if dir_name_lower.startswith(key.lower()):
                    focal_species = species_dir_to_focal[key]
                    break
            if focal_species is None:
                print(
                    f"ERROR: no subfolders matching '*_stitched_nogaps' found in {batch_parent} "
                    f"and cannot infer focal species from directory name '{batch_parent.name}'.",
                    file=sys.stderr,
                )
                return 2
            prefix = batch_parent.name
            out_per  = outdir / f"{prefix}-synteny-{thr_label}.tsv"
            out_gene = outdir / f"{prefix}-synteny-{thr_label}_summary.tsv"
            out_anc  = outdir / f"{prefix}-synteny-{thr_label}_AncSeq.tsv"
            print(f"[BATCH] {batch_parent} (flat folder) -> focal={focal_species}", file=sys.stderr)
            return run_one(batch_parent, out_per, out_gene, out_anc, focal_species)

        status = 0
        for subdir in subdirs:
            species_key = subdir.name[:-len("_stitched_nogaps")]
            focal_species = species_dir_to_focal.get(species_key)
            if focal_species is None:
                print(f"WARNING: skipping {subdir.name}; cannot infer focal species from prefix '{species_key}'.", file=sys.stderr)
                status = 1
                continue

            prefix = subdir.name
            out_per = outdir / f"{prefix}-synteny-{thr_label}.tsv"
            out_gene = outdir / f"{prefix}-synteny-{thr_label}_summary.tsv"
            out_anc = outdir / f"{prefix}-synteny-{thr_label}_AncSeq.tsv"

            print(f"[BATCH] {subdir} -> focal={focal_species}", file=sys.stderr)
            rc = run_one(subdir, out_per, out_gene, out_anc, focal_species)
            if rc != 0:
                status = rc

        return status

    # ============================================================
    # Single-folder mode
    # ============================================================
    missing = [
        name for name, val in [
            ("input_folder", args.input_folder),
            ("out_per_seq_tsv", args.out_per_seq_tsv),
            ("out_by_gene_tsv", args.out_by_gene_tsv),
            ("out_by_ancseq_tsv", args.out_by_ancseq_tsv),
            ("--focal-species", args.focal_species),
        ]
        if not val
    ]
    if missing:
        print(f"ERROR: missing required single-folder arguments: {', '.join(missing)}", file=sys.stderr)
        return 2

    return run_one(
        Path(args.input_folder),
        Path(args.out_per_seq_tsv),
        Path(args.out_by_gene_tsv),
        Path(args.out_by_ancseq_tsv),
        args.focal_species,
    )


if __name__ == "__main__":
    raise SystemExit(main())
