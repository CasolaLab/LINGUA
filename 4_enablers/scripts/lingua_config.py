"""
lingua_config.py

Tree-specific settings shared by the LINGUA scripts, read from data/tree_config.json
(written by make_tree_config.py). Set LINGUA_CONFIG to use a different file.
"""

import json
import os
import sys

_default = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "data", "tree_config.json")
_path = os.environ.get("LINGUA_CONFIG", _default)
if not os.path.isfile(_path):
    sys.exit(f"ERROR: tree configuration not found: {_path} (run make_tree_config.py)")
CFG = json.load(open(_path))

SPECIES = list(CFG["species"])
ALIASES = dict(CFG["aliases"])
NODE_ORDER_YOUNG_TO_OLD = list(CFG["node_order_young_to_old"])
NODE_CLADES = {k: list(v) for k, v in CFG["node_clades"].items()}
EXCLUDED_ORF_NODES = set(CFG["excluded_orf_nodes"])
EXCLUDED_ENABLER_NODES = set(CFG["excluded_enabler_nodes"])
CONDITIONAL_ENABLER_NODES = {k: list(v) for k, v in CFG.get("conditional_enabler_nodes", {}).items()}
STOPS_FOCAL_FRAME_ONLY = bool(CFG.get("stops_focal_frame_only", False))
# Method options: data/method_options.json, or LINGUA_OPTIONS. Each key switches one rule on (README.md, "Method
# options"); a key that is absent is off, except conserved_orf_filter, which run_pipeline.sh treats as on.
_opt_default = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "data", "method_options.json")
_opt_path = os.environ.get("LINGUA_OPTIONS", _opt_default)
OPTIONS = json.load(open(_opt_path)) if os.path.isfile(_opt_path) else {}


def opt(name, default=False):
    return OPTIONS.get(name, default)


def orf_node_from_cells(has_orf, lineage_old_to_young, specific_label):
    """ORF node from per-ancestor ORF calls (has_orf(node) -> bool). Default: the oldest focal-lineage ancestor with the
    ORF. With orf_node_contiguous (T2T 2026-11-13): the oldest ancestor of the unbroken run of ORF-bearing ancestors
    above the focal species, so an old ancestor whose ORF is not carried by the younger lineage ancestors cannot date
    the gene; <focal>-specific when the youngest ancestor lacks the ORF."""
    if not opt("orf_node_contiguous"):
        return next((n for n in lineage_old_to_young if has_orf(n)), specific_label)
    node = specific_label
    for n in reversed(list(lineage_old_to_young)):
        if not has_orf(n):
            break
        node = n
    return node


ABBREV = {sp: a for a, sp in ALIASES.items() if a != sp and a != sp.lower() and not a.islower()}


def canonicalize(name):
    s = str(name).strip()
    return ALIASES.get(s, ALIASES.get(s.lower(), s))


def side_branch(node, lineage_old_to_young, focal):
    """Living species on the side branch of a focal-lineage node: its descendants that do not descend from the next
    younger lineage node (the focal species itself when the node is the youngest)."""
    i = list(lineage_old_to_young).index(node)
    younger = set(NODE_CLADES[lineage_old_to_young[i + 1]]) if i + 1 < len(lineage_old_to_young) else {focal}
    return set(NODE_CLADES[node]) - younger


def fill_single_gaps(present, lineage_old_to_young, living_carriers, focal):
    """event_gap_living_inferred (T2T 2026-11-13): a single lineage ancestor that lacks an event between two carrier
    ancestors is counted as a carrier when a living species on its side branch carries the same event (same position;
    indels same length). Returns (carriers old -> young, inferred nodes)."""
    lin = list(lineage_old_to_young)
    idx = sorted(lin.index(n) for n in present)
    out, inferred = [lin[idx[0]]], []
    for a, b in zip(idx, idx[1:]):
        if b == a + 2:
            g = lin[a + 1]
            if side_branch(g, lin, focal) & set(living_carriers):
                out.append(g)
                inferred.append(g)
        out.append(lin[b])
    return out, inferred
