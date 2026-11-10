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
# Method options (lingua_mammals_pipeline): data/method_options.json, or LINGUA_OPTIONS. Absent or
# empty means the baseline behaviour of 2026-10-09; each key switches on one improvement.
_opt_default = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "data", "method_options.json")
_opt_path = os.environ.get("LINGUA_OPTIONS", _opt_default)
OPTIONS = json.load(open(_opt_path)) if os.path.isfile(_opt_path) else {}


def opt(name, default=False):
    return OPTIONS.get(name, default)


ABBREV = {sp: a for a, sp in ALIASES.items() if a != sp and a != sp.lower() and not a.islower()}


def canonicalize(name):
    s = str(name).strip()
    return ALIASES.get(s, ALIASES.get(s.lower(), s))
