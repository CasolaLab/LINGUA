#!/usr/bin/env python3
"""
Usage:
python3 node_summary_tree_v2.py \
  --batch \
  --tblastn-dir cdngs_tblastn_results_Apr-2026 \
  --synteny-dir synteny_40perc \
  --indels-dir indels_dir \
  --stops-dir Stops_dir \
  --tree SpeciesTree_rooted_node_labels_oct_16.txt \
  --outdir hyper_summary/nodes_summary_by_species


python3 node_summary_tree_v2.py \
    cdngs_tblastn_results_Apr-2026/athaliana_cdngs_tblastn_1078_AncSeq.tsv \
    synteny_40perc/athaliana_stitched_nogaps-synteny-40_AncSeq.tsv \
    indels_dir/all_indels_nodes_Arabidopsis_thaliana.tsv \
    Stops_dir/all_stops_nodes_Arabidopsis_thaliana.tsv \
    --tree SpeciesTree_rooted_node_labels_oct_16.txt \
    --focal-species Arabidopsis_thaliana


node_summary_tree_v2.py

Tree-aware synthesis of four summaries to classify candidate genes as DNG or PDNG
and assign their ORF, synteny, indel/stop enabler, and final enabler nodes.

This script uses the same phylogenetic framework as node_summary_tree.py:
- Newick parser for labeled internal nodes, e.g. N0, N1, N2...
- focal-species lineage derivation from the tree
- species-name canonicalization/alias handling
- tolerant geneID matching, including terminal '.p' variants

INPUT FILES
-----------
1. <focal_species>_cdngs_tblastn_1078_AncSeq.tsv
   Required columns:
     geneID
     Node_AncSeq_ORF

2. <focal_species>_stitched_nogaps-synteny-40_AncSeq.tsv
   Required columns:
     geneID
     Node_Synteny_AncSeq-40
   The script also accepts any column starting with 'Node_Synteny_AncSeq-'
   if the exact -40 column is not present.

3. all_indels_nodes_<focal_species>.tsv
   Required columns:
     geneID or gff_transcriptID
     IndelStart_RefPos
     Nodes_indels
   The file may contain multiple rows per gene. For each gene, Indels_Node is
   selected from the row whose Nodes_indels value contains the node nearest to
   ORF_Node while still older than ORF_Node; rows where Nodes_indels contains
   the ORF_Node are excluded from consideration; the selected row's remaining
   nodes are reported after synteny-age filtering and N0/N1 exclusion.

4. all_stops_nodes_<focal_species>.tsv
   Required columns:
     geneID or gff_transcriptID
     IndelStart_RefPos or StopStart_RefPos or equivalent position column
     Nodes_stops
   The file may contain multiple rows per gene. For each gene, Stops_Node is
   selected from the row whose Nodes_stops value contains the node nearest to
   ORF_Node while still older than ORF_Node; rows where Nodes_stops contains
   the ORF_Node are excluded from consideration; the selected row's remaining
   nodes are reported after synteny-age filtering and N0/N1 exclusion.

OUTPUT
------
By default:
  <canonical_focal_species>_nodes_summary.tsv

Columns:
  geneID
  ORF_Node
  Synteny_Node
  Indels_Node
  Stops_Node
  Enablers_Nodes
  Gene_type
  ORF-Enablers_Node_dist

RULES
-----
- ORF_Node comes from Node_AncSeq_ORF.
- Synteny_Node comes from Node_Synteny_AncSeq-40.
- Indels_Node is selected from all Nodes_indels rows for the gene by choosing
  the shared indel whose node list contains the nearest node to ORF_Node while
  still being older than ORF_Node. After row selection, nodes older than
  Synteny_Node are filtered out and N0/N1 are excluded.
- Stops_Node is selected from all Nodes_stops rows for the gene by choosing
  the shared stop whose node list contains the nearest node to ORF_Node while
  still being older than ORF_Node. After row selection, nodes older than
  Synteny_Node are filtered out and N0/N1 are excluded.
- Enablers_Nodes contains the two nearest nodes to ORF_Node among the filtered
  Indels_Node and Stops_Node values, but only nodes older than ORF_Node are
  considered. N0 and N1 are excluded. The selected two nodes are printed
  older-to-younger.
- Gene_type is DNG if ORF_Node is younger than both nodes in Enablers_Nodes and ORF-Enablers_Node_dist <= 1.
- Gene_type is PDNG otherwise, including cases with fewer than two valid
  enabler nodes or ORF-Enablers_Node_dist > 1.
- Rows with Gene_type = PDNG and ORF_Node equal to N0, N1, N2, or N3 are
  removed from the output file. DNG rows are never removed by this rule.
- ORF-Enablers_Node_dist is the distance in lineage nodes between ORF_Node
  and the youngest node in Enablers_Nodes. It is calculated as the age index
  of the youngest enabler node minus the age index of ORF_Node. NA is reported
  if ORF_Node or Enablers_Nodes cannot be placed on the focal lineage.

Usage:
  python3 node_summary_tree_v2.py \
    cdngs_tblastn_results_Apr-2026/athaliana_cdngs_tblastn_1078_AncSeq.tsv \
    synteny_40perc/athaliana_stitched_nogaps-synteny-40_AncSeq.tsv \
    indels_dir/all_indels_nodes_Arabidopsis_thaliana.tsv \
    Stops_dir/all_stops_nodes_Arabidopsis_thaliana.tsv \
    --tree SpeciesTree_rooted_node_labels_oct_16.txt \
    --focal-species Arabidopsis_thaliana

Optional single-species output:
  --out my_output.tsv

Batch directory mode:
  python3 node_summary_tree_v2_batch.py \
    --batch \
    --tblastn-dir cdngs_tblastn_results_Apr-2026 \
    --synteny-dir synteny_40perc \
    --indels-dir indels_dir \
    --stops-dir Stops_dir \
    --tree SpeciesTree_rooted_node_labels_oct_16.txt \
    --outdir nodes_summary_by_species

Batch mode infers focal species from filenames and writes one output file per species.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any


NODE_RE = re.compile(r"^N\d+$", flags=re.IGNORECASE)
EXCLUDED_ENABLER_NODES = {"N0", "N1"}
EXCLUDED_PDNG_ORF_NODES = {"N0", "N1", "N2", "N3"}


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
        "Camara": "Cardamine_amara_amara",
        "Camara_amara": "Cardamine_amara_amara",
        "Tarvense": "Thlaspi_arvense_MN106",
        "TarvensevarMN106": "Thlaspi_arvense_MN106",
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
        # common lowercase short file prefixes
        "athaliana": "Arabidopsis_thaliana",
        "alyrata": "Arabidopsis_lyrata",
        "ahalleri": "Arabidopsis_halleri",
        "bstricta": "Boechera_stricta",
        "cgrandiflora": "Capsella_grandiflora",
        "crubella": "Capsella_rubella",
        "cneglecta": "Camelina_neglecta",
        "csativa": "Camelina_sativa",
        "camaraamara": "Cardamine_amara_amara",
        "camara": "Cardamine_amara_amara",
        "camara_amara": "Cardamine_amara_amara",
        "tarvense": "Thlaspi_arvense_MN106",
        "tarvensevarmn106": "Thlaspi_arvense_MN106",
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


def derive_node_to_tips(root: Node) -> Dict[str, Set[str]]:
    node_to_tips: Dict[str, Set[str]] = {}

    def dfs(n: Node) -> Set[str]:
        if n.is_leaf():
            return {canonicalize_species(n.name)}
        tips: Set[str] = set()
        for c in n.children:
            tips |= dfs(c)
        if NODE_RE.match(n.name or ""):
            node_to_tips[n.name] = tips.copy()
        return tips

    dfs(root)
    return node_to_tips


def derive_species_to_ancestors(root: Node) -> Dict[str, List[str]]:
    leaves: List[Node] = []

    def collect(n: Node) -> None:
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
                anc.append(cur.name)  # youngest -> oldest
            cur = cur.parent
        out[canonicalize_species(leaf.name)] = anc
    return out


# ============================================================
# General helpers
# ============================================================
def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def normalize_geneid(gid: str) -> str:
    """
    Normalize gene/protein IDs across files so lookups tolerate common suffix variants.

    Current rules:
    - strip whitespace
    - remove one or more final '.p' protein suffixes
    - strip coordinate suffix _start_end_len (e.g. AT1G05055_1450824_1450132_96 -> AT1G05055)
      so headers from nested-ORF FASTAs match clean IDs from filename-derived lookups

    Examples:
      Bostr.0124s0010.1.p          -> Bostr.0124s0010.1
      AT1G01010.1.p                -> AT1G01010.1
      AT1G05055_1450824_1450132_96 -> AT1G05055
      AT1G04120.1_1070340_1070245_32 -> AT1G04120.1
    """
    x = (gid or "").strip()
    x = re.sub(r"(\.p)+$", "", x)
    x = re.sub(r"_\d+_\d+_\d+$", "", x)
    return x


def get_focal_specific_label(focal_species: str) -> str:
    return f"{canonicalize_species(focal_species)}-specific"


def normalize_node_label(label: str, focal_species: str) -> str:
    """
    Normalize node labels and focal-specific labels.

    Examples:
      Athaliana-specific -> Arabidopsis_thaliana-specific
      athaliana-specific -> Arabidopsis_thaliana-specific
      N9                 -> N9
      NA / blank          -> NA
    """
    x = (label or "").strip()
    if not x:
        return "NA"
    if x.lower() in {"na", "nan", "none", "."}:
        return "NA"
    if x.endswith("-specific"):
        prefix = x[:-9]
        prefix = canonicalize_species(prefix)
        if prefix == canonicalize_species(focal_species):
            return get_focal_specific_label(focal_species)
        return f"{prefix}-specific"
    return x


def split_node_list(value: str, focal_species: str) -> List[str]:
    """Parse Nodes_indels / Nodes_stops values separated by |, comma, semicolon, or whitespace."""
    x = (value or "").strip()
    if not x or x.lower() in {"na", "nan", "none", ".", "no-indels", "no-stops"}:
        return []
    parts = re.split(r"[|,;\s]+", x)
    out: List[str] = []
    seen: Set[str] = set()
    for p in parts:
        lab = normalize_node_label(p, focal_species)
        if lab == "NA":
            continue
        if lab not in seen:
            out.append(lab)
            seen.add(lab)
    return out


def node_age_index(label: str, focal_species: str, ancestors_young_to_old: List[str]) -> Optional[int]:
    """
    Convert a node label to a lineage age index.

    Index scale:
      0 = focal-specific / youngest
      1 = youngest internal lineage node
      2 = next older node
      ...

    Higher index = older.
    Returns None for labels not on focal lineage or NA.
    """
    lab = normalize_node_label(label, focal_species)
    if lab == "NA":
        return None
    if lab == get_focal_specific_label(focal_species):
        return 0
    if lab in ancestors_young_to_old:
        return ancestors_young_to_old.index(lab) + 1
    return None


def is_older(node_a: str, node_b: str, focal_species: str, ancestors_young_to_old: List[str]) -> bool:
    """True if node_a is older than node_b on the focal lineage."""
    ia = node_age_index(node_a, focal_species, ancestors_young_to_old)
    ib = node_age_index(node_b, focal_species, ancestors_young_to_old)
    if ia is None or ib is None:
        return False
    return ia > ib


def filter_nodes_not_older_than_synteny(
    nodes_value: str,
    synteny_node: str,
    focal_species: str,
    ancestors_young_to_old: List[str],
) -> str:
    """
    Keep only nodes from nodes_value that are the same age as, or younger than,
    Synteny_Node. Nodes older than Synteny_Node are removed. N0 and N1 are
    also removed and therefore are not reported or used for enabler calls.
    """
    syn_idx = node_age_index(synteny_node, focal_species, ancestors_young_to_old)
    if syn_idx is None:
        return "NA"

    kept: List[Tuple[int, str]] = []
    seen: Set[str] = set()
    for lab in split_node_list(nodes_value, focal_species):
        if lab in EXCLUDED_ENABLER_NODES:
            continue
        if lab in seen:
            continue
        seen.add(lab)
        idx = node_age_index(lab, focal_species, ancestors_young_to_old)
        if idx is None:
            continue
        if idx <= syn_idx:
            kept.append((idx, lab))

    if not kept:
        return "NA"

    kept.sort(key=lambda x: x[0], reverse=True)
    return "|".join(lab for _, lab in kept)


def select_best_nodes_row(
    existing_value: str,
    candidate_value: str,
    focal_species: str,
    ancestors_young_to_old: List[str],
) -> str:
    """
    Select the best Nodes_indels / Nodes_stops row for a gene.

    Rule:
    - compare the oldest node in each row
    - choose the row with the older oldest node
    - if tied, choose the row with the most nodes
    - keep the existing row if still tied
    """
    cand_nodes = split_node_list(candidate_value, focal_species)
    if not cand_nodes:
        return existing_value

    if not existing_value or existing_value == "NA":
        return "|".join(cand_nodes)

    old_nodes = split_node_list(existing_value, focal_species)
    if not old_nodes:
        return "|".join(cand_nodes)

    old_idxs = [node_age_index(n, focal_species, ancestors_young_to_old) for n in old_nodes]
    cand_idxs = [node_age_index(n, focal_species, ancestors_young_to_old) for n in cand_nodes]
    old_idxs = [i for i in old_idxs if i is not None]
    cand_idxs = [i for i in cand_idxs if i is not None]

    if not cand_idxs:
        return existing_value
    if not old_idxs:
        return "|".join(cand_nodes)

    old_oldest = max(old_idxs)
    cand_oldest = max(cand_idxs)

    if cand_oldest > old_oldest:
        return "|".join(cand_nodes)
    if cand_oldest < old_oldest:
        return existing_value

    if len(cand_nodes) > len(old_nodes):
        return "|".join(cand_nodes)
    return existing_value





def select_nodes_row_nearest_to_orf(
    candidate_values: Any,
    orf_node: str,
    focal_species: str,
    ancestors_young_to_old: List[str],
) -> str:
    """
    Select one Nodes_indels / Nodes_stops row for a gene.

    New rule:
    - consider each shared indel/stop row separately
    - skip any row whose node list contains ORF_Node
    - for each remaining row, find the node in that row that is older than ORF_Node
      but nearest to ORF_Node
    - choose the row with the nearest such node
    - if tied, choose the row with the most nodes
    - keep the existing row if still tied

    The returned value is the full selected row's node list. Synteny filtering
    and N0/N1 exclusion are applied later by filter_nodes_not_older_than_synteny().
    """
    if candidate_values is None or candidate_values == "NA":
        return "NA"
    if isinstance(candidate_values, str):
        values = [candidate_values]
    else:
        values = list(candidate_values)

    orf_idx = node_age_index(orf_node, focal_species, ancestors_young_to_old)
    if orf_idx is None:
        return "NA"

    best_value = "NA"
    best_nearest_idx: Optional[int] = None
    best_n_nodes = -1

    for value in values:
        nodes = split_node_list(value, focal_species)
        if not nodes:
            continue
        if normalize_node_label(orf_node, focal_species) in nodes:
            continue

        idxs: List[int] = []
        for n in nodes:
            if n in EXCLUDED_ENABLER_NODES:
                continue
            idx = node_age_index(n, focal_species, ancestors_young_to_old)
            if idx is not None and idx > orf_idx:
                idxs.append(idx)

        if not idxs:
            continue

        # nearest older node = smallest age index greater than ORF_Node index
        nearest_idx = min(idxs)
        n_nodes = len(nodes)

        if best_nearest_idx is None or nearest_idx < best_nearest_idx:
            best_value = "|".join(nodes)
            best_nearest_idx = nearest_idx
            best_n_nodes = n_nodes
        elif nearest_idx == best_nearest_idx and n_nodes > best_n_nodes:
            best_value = "|".join(nodes)
            best_n_nodes = n_nodes

    return best_value
def choose_enabler_nodes(
    indels_node: str,
    stops_node: str,
    orf_node: str,
    focal_species: str,
    ancestors_young_to_old: List[str],
) -> str:
    """
    Select Enablers_Nodes from already-filtered Indels_Node and Stops_Node values.

    N0 and N1 are ignored here as a defensive check, although they should already
    have been removed from Indels_Node and Stops_Node.

    Enablers_Nodes = the two nodes nearest to ORF_Node, but older than ORF_Node.
    The two selected nodes are printed older-to-younger.
    """
    orf_idx = node_age_index(orf_node, focal_species, ancestors_young_to_old)
    if orf_idx is None:
        return "NA"

    candidates = split_node_list(indels_node, focal_species) + split_node_list(stops_node, focal_species)

    valid: List[Tuple[int, str]] = []
    seen: Set[str] = set()
    for lab in candidates:
        if lab in EXCLUDED_ENABLER_NODES:
            continue
        if lab in seen:
            continue
        seen.add(lab)
        idx = node_age_index(lab, focal_species, ancestors_young_to_old)
        if idx is None:
            continue
        if idx > orf_idx:
            valid.append((idx, lab))

    if len(valid) < 2:
        return "NA"

    # First choose the two nearest older nodes: smallest indices greater than ORF index.
    valid.sort(key=lambda x: x[0])
    nearest_two = valid[:2]

    # Print older-to-younger, matching previous enabler-node output convention.
    nearest_two.sort(key=lambda x: x[0], reverse=True)
    return "|".join(lab for _, lab in nearest_two)


def is_orf_younger_than_all_enabler_nodes(
    orf_node: str,
    enablers_nodes: str,
    focal_species: str,
    ancestors_young_to_old: List[str],
) -> bool:
    """True only if ORF_Node is younger than both nodes in Enablers_Nodes."""
    nodes = split_node_list(enablers_nodes, focal_species)
    if len(nodes) < 2:
        return False
    return all(is_older(en_node, orf_node, focal_species, ancestors_young_to_old) for en_node in nodes[:2])


def calculate_orf_enablers_node_dist(
    orf_node: str,
    enablers_nodes: str,
    focal_species: str,
    ancestors_young_to_old: List[str],
) -> str:
    """
    Calculate the distance in lineage nodes between ORF_Node and the youngest
    node in Enablers_Nodes.

    Examples for Arabidopsis_thaliana lineage:
      ORF_Node=N9, Enablers_Nodes=N4|N6 -> 1
      ORF_Node=N9, Enablers_Nodes=N3|N4 -> 2
      ORF_Node=Arabidopsis_thaliana-specific, Enablers_Nodes=N6|N9 -> 1
      ORF_Node=Arabidopsis_thaliana-specific, Enablers_Nodes=N2|N3 -> 4
    """
    orf_idx = node_age_index(orf_node, focal_species, ancestors_young_to_old)
    if orf_idx is None:
        return "NA"

    enabler_idxs: List[int] = []
    for lab in split_node_list(enablers_nodes, focal_species):
        idx = node_age_index(lab, focal_species, ancestors_young_to_old)
        if idx is not None and idx > orf_idx:
            enabler_idxs.append(idx)

    if not enabler_idxs:
        return "NA"

    youngest_enabler_idx = min(enabler_idxs)
    return str(youngest_enabler_idx - orf_idx)


def get_required_col(fieldnames: List[str], exact: str, path: Path) -> str:
    if exact in fieldnames:
        return exact
    die(f"Required column '{exact}' not found in {path}. Available columns: {', '.join(fieldnames)}")
    return exact


def get_first_existing_col(fieldnames: List[str], candidates: List[str], path: Path) -> str:
    for c in candidates:
        if c in fieldnames:
            return c
    die(f"None of these columns found in {path}: {', '.join(candidates)}. Available columns: {', '.join(fieldnames)}")
    return candidates[0]


def detect_synteny_ancseq_col(fieldnames: List[str], path: Path) -> str:
    if "Node_Synteny_AncSeq-40" in fieldnames:
        return "Node_Synteny_AncSeq-40"
    candidates = [c for c in fieldnames if c.startswith("Node_Synteny_AncSeq-")]
    if candidates:
        if len(candidates) > 1:
            print(
                "WARNING: multiple Node_Synteny_AncSeq-* columns found; using the first: "
                f"{candidates[0]}",
                file=sys.stderr,
            )
        return candidates[0]
    die(f"Could not find Node_Synteny_AncSeq-40 or Node_Synteny_AncSeq-* in {path}.")
    return "Node_Synteny_AncSeq-40"


def insert_lookup_value(lookup: Dict[str, Any], gid: str, value: Any) -> None:
    gid_exact = (gid or "").strip()
    gid_norm = normalize_geneid(gid_exact)
    if gid_exact and gid_exact not in lookup:
        lookup[gid_exact] = value
    if gid_norm and gid_norm not in lookup:
        lookup[gid_norm] = value


def lookup_gene(lookup: Dict[str, Any], gid: str, default: Any = "NA") -> Any:
    gid_exact = (gid or "").strip()
    gid_norm = normalize_geneid(gid_exact)
    if gid_exact in lookup:
        return lookup[gid_exact]
    if gid_norm in lookup:
        return lookup[gid_norm]
    return default


def read_orf_table(path: Path, focal_species: str) -> List[Tuple[str, str]]:
    rows: List[Tuple[str, str]] = []
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            die(f"Could not read header from {path}")
        gene_col = get_first_existing_col(reader.fieldnames, ["geneID", "gff_transcriptID"], path)
        node_col = get_required_col(reader.fieldnames, "Node_AncSeq_ORF", path)
        for row in reader:
            gid = (row.get(gene_col, "") or "").strip()
            if not gid:
                continue
            node = normalize_node_label(row.get(node_col, ""), focal_species)
            rows.append((gid, node))
    if not rows:
        die(f"No gene rows read from ORF/tblastn file: {path}")
    return rows


def read_single_node_lookup(path: Path, node_col: str, focal_species: str, *, synteny: bool = False) -> Dict[str, str]:
    lookup: Dict[str, str] = {}
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            die(f"Could not read header from {path}")
        gene_col = get_first_existing_col(reader.fieldnames, ["geneID", "gff_transcriptID"], path)
        use_node_col = detect_synteny_ancseq_col(reader.fieldnames, path) if synteny else get_required_col(reader.fieldnames, node_col, path)
        for row in reader:
            gid = (row.get(gene_col, "") or "").strip()
            if not gid:
                continue
            val = normalize_node_label(row.get(use_node_col, ""), focal_species)
            insert_lookup_value(lookup, gid, val)
    return lookup


def read_multi_node_lookup(
    path: Path,
    nodes_col: str,
    focal_species: str,
    ancestors_young_to_old: List[str],
) -> Dict[str, List[str]]:
    lookup: Dict[str, List[str]] = {}
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            die(f"Could not read header from {path}")
        gene_col = get_first_existing_col(reader.fieldnames, ["geneID", "gff_transcriptID"], path)
        use_nodes_col = get_required_col(reader.fieldnames, nodes_col, path)
        for row in reader:
            gid = (row.get(gene_col, "") or "").strip()
            if not gid:
                continue
            nodes = split_node_list(row.get(use_nodes_col, ""), focal_species)
            val = "|".join(nodes) if nodes else "NA"
            if val == "NA":
                continue

            gid_exact = (gid or "").strip()
            gid_norm = normalize_geneid(gid_exact)

            if gid_exact:
                lookup.setdefault(gid_exact, []).append(val)
            if gid_norm:
                lookup.setdefault(gid_norm, []).append(val)
    return lookup


# ============================================================
# Single-species runner + batch discovery
# ============================================================
def write_one_summary(
    tblastn_path: Path,
    synteny_path: Path,
    indels_path: Path,
    stops_path: Path,
    tree_path: Path,
    focal_species: str,
    out_path: Path,
) -> None:
    """Run the existing single-species synthesis logic for one focal species."""
    focal = canonicalize_species(focal_species)

    if not tree_path.exists():
        die(f"Tree file not found: {tree_path}")

    root = parse_newick(tree_path)
    # Kept because this is part of the same phylogenetic framework and useful for future checks.
    _node_to_tips = derive_node_to_tips(root)
    species_to_ancestors = derive_species_to_ancestors(root)

    if focal not in species_to_ancestors:
        die(f"Focal species not found in tree tips: {focal}")

    ancestors_young_to_old = species_to_ancestors[focal]
    if not ancestors_young_to_old:
        die(f"No labeled ancestor nodes found for focal species: {focal}")

    for p in [tblastn_path, synteny_path, indels_path, stops_path]:
        if not p.exists():
            die(f"Input file not found: {p}")

    orf_rows = read_orf_table(tblastn_path, focal)
    synteny_lookup = read_single_node_lookup(synteny_path, "Node_Synteny_AncSeq-40", focal, synteny=True)
    indels_lookup = read_multi_node_lookup(indels_path, "Nodes_indels", focal, ancestors_young_to_old)
    stops_lookup = read_multi_node_lookup(stops_path, "Nodes_stops", focal, ancestors_young_to_old)

    out_cols = [
        "geneID",
        "ORF_Node",
        "Synteny_Node",
        "Indels_Node",
        "Stops_Node",
        "Enablers_Nodes",
        "Gene_type",
        "ORF-Enablers_Node_dist",
    ]

    n_syn_match = 0
    n_indel_match = 0
    n_stop_match = 0
    n_enabler = 0
    n_dng = 0
    n_pdng = 0
    n_pdng_old_orf_removed = 0

    with out_path.open("w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, delimiter="\t", fieldnames=out_cols)
        writer.writeheader()

        for gid, orf_node in orf_rows:
            synteny_node = lookup_gene(synteny_lookup, gid, "NA")
            indels_candidates = lookup_gene(indels_lookup, gid, "NA")
            stops_candidates = lookup_gene(stops_lookup, gid, "NA")

            indels_node = select_nodes_row_nearest_to_orf(
                candidate_values=indels_candidates,
                orf_node=orf_node,
                focal_species=focal,
                ancestors_young_to_old=ancestors_young_to_old,
            )
            stops_node = select_nodes_row_nearest_to_orf(
                candidate_values=stops_candidates,
                orf_node=orf_node,
                focal_species=focal,
                ancestors_young_to_old=ancestors_young_to_old,
            )

            if synteny_node != "NA":
                n_syn_match += 1
            if indels_node != "NA":
                n_indel_match += 1
            if stops_node != "NA":
                n_stop_match += 1

            indels_node = filter_nodes_not_older_than_synteny(
                nodes_value=indels_node,
                synteny_node=synteny_node,
                focal_species=focal,
                ancestors_young_to_old=ancestors_young_to_old,
            )
            stops_node = filter_nodes_not_older_than_synteny(
                nodes_value=stops_node,
                synteny_node=synteny_node,
                focal_species=focal,
                ancestors_young_to_old=ancestors_young_to_old,
            )

            enablers_nodes = choose_enabler_nodes(
                indels_node=indels_node,
                stops_node=stops_node,
                orf_node=orf_node,
                focal_species=focal,
                ancestors_young_to_old=ancestors_young_to_old,
            )
            if enablers_nodes != "NA":
                n_enabler += 1

            orf_enablers_node_dist = calculate_orf_enablers_node_dist(
                orf_node=orf_node,
                enablers_nodes=enablers_nodes,
                focal_species=focal,
                ancestors_young_to_old=ancestors_young_to_old,
            )

            dist_ok = False
            if orf_enablers_node_dist != "NA":
                try:
                    dist_ok = int(orf_enablers_node_dist) <= 1
                except ValueError:
                    dist_ok = False

            gene_type = "DNG" if (
                is_orf_younger_than_all_enabler_nodes(
                    orf_node,
                    enablers_nodes,
                    focal,
                    ancestors_young_to_old,
                )
                and dist_ok
            ) else "PDNG"
            if gene_type == "DNG":
                n_dng += 1
            else:
                n_pdng += 1

            if gene_type == "PDNG" and normalize_node_label(orf_node, focal) in EXCLUDED_PDNG_ORF_NODES:
                n_pdng_old_orf_removed += 1
                continue

            writer.writerow(
                {
                    "geneID": gid,
                    "ORF_Node": orf_node,
                    "Synteny_Node": synteny_node,
                    "Indels_Node": indels_node,
                    "Stops_Node": stops_node,
                    "Enablers_Nodes": enablers_nodes,
                    "Gene_type": gene_type,
                    "ORF-Enablers_Node_dist": orf_enablers_node_dist,
                }
            )

    total = len(orf_rows)
    print(f"Wrote: {out_path}", file=sys.stderr)
    print(f"DEBUG: focal species = {focal}", file=sys.stderr)
    print(f"DEBUG: focal lineage youngest->oldest = {'|'.join(ancestors_young_to_old)}", file=sys.stderr)
    print(f"DEBUG: synteny matches = {n_syn_match} / {total}", file=sys.stderr)
    print(f"DEBUG: indel matches = {n_indel_match} / {total}", file=sys.stderr)
    print(f"DEBUG: stop matches = {n_stop_match} / {total}", file=sys.stderr)
    print(f"DEBUG: genes with valid Enablers_Nodes = {n_enabler} / {total}", file=sys.stderr)
    print(f"DEBUG: DNG = {n_dng}; PDNG = {n_pdng}", file=sys.stderr)
    print(f"DEBUG: PDNG rows removed because ORF_Node is N0/N1/N2/N3 = {n_pdng_old_orf_removed}", file=sys.stderr)


def species_prefix_for_output(focal_species: str, fallback: str) -> str:
    """Return a compact output prefix when possible; otherwise use the canonical species name."""
    focal = canonicalize_species(focal_species)
    short_aliases = {
        "Arabidopsis_thaliana": "athaliana",
        "Arabidopsis_lyrata": "alyrata",
        "Arabidopsis_halleri": "ahalleri",
        "Boechera_stricta": "bstricta",
        "Capsella_grandiflora": "cgrandiflora",
        "Capsella_rubella": "crubella",
        "Camelina_neglecta": "cneglecta",
        "Camelina_sativa": "csativa",
        "Cardamine_amara_amara": "camaraamara",
        "Arabis_alpina": "aalpina",
        "Arabis_nemorensis": "anemorensis",
        "Brassica_oleracea": "boleracea",
        "Brassica_rapa_FPsc": "brapa",
        "Eutrema_salsugineum": "esalsugineum",
        "Isatis_indigotica": "iindigotica",
        "Microthlaspi_erraticum": "merraticum",
        "Raphanus_sativus": "rsativus",
        "Schrenkiella_parvula": "sparvula",
        "Sinapis_alba": "salba",
        "Thlaspi_arvense_MN106": "tarvensevarmn106",
        "Aethionema_arabicum": "aarabicum",
        "Tarenaya_hassleriana": "thassleriana",
        "Theobroma_cacao": "tcacao",
    }
    return short_aliases.get(focal, fallback or focal)


def infer_focal_from_tblastn_filename(path: Path) -> Optional[Tuple[str, str]]:
    """
    Infer species from files such as:
      athaliana_cdngs_tblastn_1078_AncSeq.tsv
      Brassica_oleracea_cdngs_tblastn_123_AncSeq.tsv
    Returns (canonical_species, file_prefix).
    """
    m = re.match(r"^(.+?)_cdngs_tblastn_.*_AncSeq\.tsv$", path.name)
    if not m:
        return None
    prefix = m.group(1)
    focal = canonicalize_species(prefix)
    return focal, prefix


def build_file_index_by_species(directory: Path, kind: str) -> Dict[str, Path]:
    """Index files in a directory by canonical focal species."""
    out: Dict[str, Path] = {}
    patterns: List[str]
    if kind == "tblastn":
        patterns = ["*_cdngs_tblastn_*_AncSeq.tsv"]
    elif kind == "synteny":
        patterns = ["*_nogaps-synteny-*_AncSeq.tsv"]
    elif kind == "indels":
        patterns = ["all_indels_nodes_*.tsv"]
    elif kind == "stops":
        patterns = ["all_stops_nodes_*.tsv"]
    else:
        die(f"Unknown file-index kind: {kind}")

    files: List[Path] = []
    for pat in patterns:
        files.extend(sorted(directory.glob(pat)))

    for path in files:
        species_token: Optional[str] = None
        if kind == "tblastn":
            m = re.match(r"^(.+?)_cdngs_tblastn_.*_AncSeq\.tsv$", path.name)
            species_token = m.group(1) if m else None
        elif kind == "synteny":
            m = re.match(r"^(.+?)_\w+_nogaps-synteny-[\d.]+_AncSeq\.tsv$", path.name)
            species_token = m.group(1) if m else None
        elif kind == "indels":
            m = re.match(r"^all_indels_nodes_(.+?)\.tsv$", path.name)
            species_token = m.group(1) if m else None
        elif kind == "stops":
            m = re.match(r"^all_stops_nodes_(.+?)\.tsv$", path.name)
            species_token = m.group(1) if m else None

        if not species_token:
            continue
        focal = canonicalize_species(species_token)
        if focal in out:
            print(
                f"WARNING: multiple {kind} files found for {focal}; keeping first: {out[focal]} and ignoring {path}",
                file=sys.stderr,
            )
            continue
        out[focal] = path

    return out


def run_batch_mode(
    tblastn_dir: Path,
    synteny_dir: Path,
    indels_dir: Path,
    stops_dir: Path,
    tree_path: Path,
    outdir: Path,
    only_species: Optional[List[str]] = None,
) -> None:
    """Run all species that can be paired across the four target directories."""
    for d in [tblastn_dir, synteny_dir, indels_dir, stops_dir]:
        if not d.exists() or not d.is_dir():
            die(f"Input directory not found: {d}")

    outdir.mkdir(parents=True, exist_ok=True)

    tblastn_index = build_file_index_by_species(tblastn_dir, "tblastn")
    synteny_index = build_file_index_by_species(synteny_dir, "synteny")
    indels_index = build_file_index_by_species(indels_dir, "indels")
    stops_index = build_file_index_by_species(stops_dir, "stops")

    requested: Optional[Set[str]] = None
    if only_species:
        requested = {canonicalize_species(x) for x in only_species}

    species_to_run = sorted(tblastn_index.keys())
    if requested is not None:
        species_to_run = [sp for sp in species_to_run if sp in requested]

    if not species_to_run:
        die("No tblastn input files found for batch mode, or none matched --only-species.")

    n_run = 0
    n_skip = 0
    for focal in species_to_run:
        missing = []
        if focal not in synteny_index:
            missing.append("synteny")
        if focal not in indels_index:
            missing.append("indels")
        if focal not in stops_index:
            missing.append("stops")
        if missing:
            n_skip += 1
            print(
                f"WARNING: skipping {focal}; missing {', '.join(missing)} file(s).",
                file=sys.stderr,
            )
            continue

        prefix_info = infer_focal_from_tblastn_filename(tblastn_index[focal])
        fallback_prefix = prefix_info[1] if prefix_info else focal
        out_prefix = species_prefix_for_output(focal, fallback_prefix)
        out_path = outdir / f"{out_prefix}_nodes_summary.tsv"

        print(f"\n=== Running {focal} ===", file=sys.stderr)
        print(f"tblastn: {tblastn_index[focal]}", file=sys.stderr)
        print(f"synteny: {synteny_index[focal]}", file=sys.stderr)
        print(f"indels : {indels_index[focal]}", file=sys.stderr)
        print(f"stops  : {stops_index[focal]}", file=sys.stderr)

        write_one_summary(
            tblastn_path=tblastn_index[focal],
            synteny_path=synteny_index[focal],
            indels_path=indels_index[focal],
            stops_path=stops_index[focal],
            tree_path=tree_path,
            focal_species=focal,
            out_path=out_path,
        )
        n_run += 1

    print(f"\nBATCH COMPLETE: wrote {n_run} species summaries; skipped {n_skip} species.", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Tree-aware synthesis of ORF, synteny, indel, and stop nodes to classify DNG vs PDNG. Supports single-species and batch directory modes."
    )

    # Original single-species positional mode. Kept optional so directory batch mode can use flags only.
    ap.add_argument("tblastn_ancseq_tsv", nargs="?", help="Single-species mode: TSV containing Node_AncSeq_ORF")
    ap.add_argument("synteny_ancseq_tsv", nargs="?", help="Single-species mode: TSV containing Node_Synteny_AncSeq-40")
    ap.add_argument("indels_nodes_tsv", nargs="?", help="Single-species mode: TSV containing Nodes_indels")
    ap.add_argument("stops_nodes_tsv", nargs="?", help="Single-species mode: TSV containing Nodes_stops")

    # Batch mode: provide directories once and run all matched species.
    ap.add_argument("--batch", action="store_true", help="Run all species by matching files across target directories.")
    ap.add_argument("--tblastn-dir", default=None, help="Batch mode: directory containing *_cdngs_tblastn_*_AncSeq.tsv files.")
    ap.add_argument("--synteny-dir", default=None, help="Batch mode: directory containing *_stitched_nogaps-synteny-40_AncSeq.tsv files.")
    ap.add_argument("--indels-dir", default=None, help="Batch mode: directory containing all_indels_nodes_<species>.tsv files.")
    ap.add_argument("--stops-dir", default=None, help="Batch mode: directory containing all_stops_nodes_<species>.tsv files.")
    ap.add_argument("--outdir", default="nodes_summary_by_species", help="Batch mode output directory. Default: nodes_summary_by_species")
    ap.add_argument("--only-species", nargs="*", default=None, help="Optional batch-mode species subset, using full names or supported aliases.")

    ap.add_argument("--tree", required=True, help="Newick tree file with labeled internal nodes, e.g. N0, N1...")
    ap.add_argument("--focal-species", default=None, help="Single-species mode: canonical focal species name or supported short alias.")
    ap.add_argument("--out", default=None, help="Single-species mode output TSV. Default: <canonical_focal_species>_nodes_summary.tsv")
    args = ap.parse_args()

    tree_path = Path(args.tree)

    if args.batch:
        needed = [args.tblastn_dir, args.synteny_dir, args.indels_dir, args.stops_dir]
        if any(x is None for x in needed):
            die("Batch mode requires --tblastn-dir, --synteny-dir, --indels-dir, and --stops-dir.")
        run_batch_mode(
            tblastn_dir=Path(args.tblastn_dir),
            synteny_dir=Path(args.synteny_dir),
            indels_dir=Path(args.indels_dir),
            stops_dir=Path(args.stops_dir),
            tree_path=tree_path,
            outdir=Path(args.outdir),
            only_species=args.only_species,
        )
        return

    positional = [args.tblastn_ancseq_tsv, args.synteny_ancseq_tsv, args.indels_nodes_tsv, args.stops_nodes_tsv]
    if any(x is None for x in positional):
        die(
            "Single-species mode requires four input TSV files, or use --batch with the four target directories."
        )
    if not args.focal_species:
        die("Single-species mode requires --focal-species. In batch mode, species are inferred from filenames.")

    focal = canonicalize_species(args.focal_species)
    out_path = Path(args.out) if args.out else Path(f"{focal}_nodes_summary.tsv")

    write_one_summary(
        tblastn_path=Path(args.tblastn_ancseq_tsv),
        synteny_path=Path(args.synteny_ancseq_tsv),
        indels_path=Path(args.indels_nodes_tsv),
        stops_path=Path(args.stops_nodes_tsv),
        tree_path=tree_path,
        focal_species=focal,
        out_path=out_path,
    )


if __name__ == "__main__":
    main()
