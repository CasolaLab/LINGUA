#!/usr/bin/env python3
"""
Usage:

Desktop:
python blast_prts-vs-nogaps_AncSeq_v2.py 45k_renamed/athaliana_cdngs_1078.faa brassicaceae_23species_LINGUA_nogaps/athaliana_stitched_nogaps cdngs_tblastn_results_Apr-2026/athaliana_cdngs_tblastn_1078.tsv cdngs_tblastn_results_Apr-2026/athaliana_cdngs_tblastn_1078_summary.tsv cdngs_tblastn_results_Apr-2026/athaliana_cdngs_tblastn_1078_AncSeq.tsv --tree SpeciesTree_rooted_node_labels_oct_16.txt --focal-species Arabidopsis_thaliana --frac_ORFss 80 --pident 40 --tblastn tblastn --makeblastdb makeblastdb --sanity-check --verbose

python blast_prts-vs-nogaps_AncSeq_v2.py 45k_renamed/aarabicum_cdngs_1143.faa brassicaceae_23species_LINGUA_nogaps/aarabicum_stitched_nogaps cdngs_tblastn_results_Apr-2026/aarabicum_cdngs_tblastn_1143.tsv cdngs_tblastn_results_Apr-2026/aarabicum_cdngs_tblastn_1143_summary.tsv cdngs_tblastn_results_Apr-2026/aarabicum_cdngs_tblastn_1143_AncSeq.tsv --tree SpeciesTree_rooted_node_labels_oct_16.txt --focal-species Aethionema_arabicum --frac_ORFss 80 --pident 40 --tblastn tblastn --makeblastdb makeblastdb --sanity-check --verbose

python blast_prts-vs-nogaps_AncSeq_v2.py 45k_renamed/merraticum_cdngs_4511.faa brassicaceae_23species_LINGUA_nogaps/merraticum_stitched_nogaps cdngs_tblastn_results_Apr-2026/merraticum_cdngs_tblastn_4511.tsv cdngs_tblastn_results_Apr-2026/merraticum_cdngs_tblastn_4511_summary.tsv cdngs_tblastn_results_Apr-2026/merraticum_cdngs_tblastn_4511_AncSeq.tsv --tree SpeciesTree_rooted_node_labels_oct_16.txt --focal-species Microthlaspi_erraticum --frac_ORFss 80 --pident 40 --tblastn tblastn --makeblastdb makeblastdb --sanity-check --verbose

python blast_prts-vs-nogaps_AncSeq_v2.py 45k_renamed/tarvensevarmn106_cdngs_330.faa brassicaceae_23species_LINGUA_nogaps/tarvense_stitched_nogaps cdngs_tblastn_results_Apr-2026/tarvense_cdngs_tblastn_330.tsv cdngs_tblastn_results_Apr-2026/tarvense_cdngs_tblastn_330_summary.tsv cdngs_tblastn_results_Apr-2026/tarvense_cdngs_tblastn_330_AncSeq.tsv --tree SpeciesTree_rooted_node_labels_oct_16.txt --focal-species Thlaspi_arvense --frac_ORFss 80 --pident 40 --tblastn tblastn --makeblastdb makeblastdb --sanity-check --verbose

python blast_prts-vs-nogaps_AncSeq_v2.py 45k_renamed/ahalleri_cdngs_932.faa brassicaceae_23species_LINGUA_nogaps/ahalleri_stitched_nogaps cdngs_tblastn_results_Apr-2026/ahalleri_cdngs_tblastn_932.tsv cdngs_tblastn_results_Apr-2026/ahalleri_cdngs_tblastn_932_summary.tsv cdngs_tblastn_results_Apr-2026/ahalleri_cdngs_tblastn_932_AncSeq.tsv --tree SpeciesTree_rooted_node_labels_oct_16.txt --focal-species Arabidopsis_halleri --frac_ORFss 80 --pident 40 --tblastn tblastn --makeblastdb makeblastdb --sanity-check --verbose


blast_prts-vs-nogaps_AncSeq_v2.py

Per-query tblastn of protein queries (multi-FASTA) against a matched subject DNA FASTA
file (one per query), matched by GeneID contained in subject filename.

This version supports:
- Output #1: per-HSP tblastn hits with derived fields (no-gaps sequences, PTC position, frac_CDS_seq, etc.)
- Output #2: per-geneID summary table by species with HSP-count thresholds, plus Node_tblastn classification.
- Output #3: per-geneID ancestral-sequence ORF presence table by node (N21..N0).

Key updates:
- subject_file column contains only the basename (file name, no path)
- frac_CDS_seq is computed as (longest stop-free segment length in sseq-NO-GAPS / (qlen-1)) * 100
- Species in output #1 is derived from sseqid prefix before the first '.' (e.g., Arabidopsis_thaliana.Chr1 -> Arabidopsis_thaliana; N9.N9ref... -> N9)
- Node_tblastn is computed tree-based and species-specific when --tree and --focal-species are provided.

Example:
  python3 blast_prts-vs-nogaps_v3.py query.faa subject_folder out.tsv summary.tsv ancseq.tsv \
    --tree SpeciesTree_rooted_node_labels_oct_16.txt \
    --focal-species Arabidopsis_thaliana \
    --frac_ORFss 80 --pident 40 \
    --tblastn /path/to/tblastn --makeblastdb /path/to/makeblastdb \
    --sanity-check --debug-nodes
"""
from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple


# -----------------------------
# Species canonicalization
# -----------------------------
# Canonical (full) species columns for summary output (23 taxa).
SPECIES_COLUMNS: List[str] = [
    "Arabidopsis_thaliana",
    "Arabidopsis_halleri",
    "Arabidopsis_lyrata",
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

# Abbrev names seen in your Newick (and sometimes in headers) -> canonical full names.
ABBREV_TO_CANON: Dict[str, str] = {
    "Athaliana": "Arabidopsis_thaliana",
    "Ahalleri": "Arabidopsis_halleri",
    "Alyrata": "Arabidopsis_lyrata",
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

CANON_SET: Set[str] = set(SPECIES_COLUMNS)

def canon_taxon(name: str) -> str:
    """Map known abbreviations to canonical names; otherwise return as-is."""
    return ABBREV_TO_CANON.get(name, name)


def resolve_species_name(name: str, species_list: List[str]) -> str:
    """Resolve user-provided species name to canonical project name.

    Accepts canonical names, abbreviations, and unique partial canonical matches.
    Examples:
      Thlaspi_arvense_MN106 -> Thlaspi_arvense_MN106
      TarvensevarMN106      -> Thlaspi_arvense_MN106
      Thlaspi_arvense       -> Thlaspi_arvense_MN106
    """
    q = (name or '').strip()
    if not q:
        raise ValueError('Empty species name')

    if q in species_list:
        return q

    if q in ABBREV_TO_CANON:
        return ABBREV_TO_CANON[q]

    # Case-insensitive exact canonical match
    canon_lut = {sp.lower(): sp for sp in species_list}
    if q.lower() in canon_lut:
        return canon_lut[q.lower()]

    # Unique prefix/substring match against canonical names
    prefix = [sp for sp in species_list if sp.lower().startswith(q.lower())]
    if len(prefix) == 1:
        return prefix[0]

    contains = [sp for sp in species_list if q.lower() in sp.lower()]
    if len(contains) == 1:
        return contains[0]

    raise ValueError(
        f"Could not resolve species name '{name}'. Valid species: " + ', '.join(sorted(species_list))
    )


# Explicit node membership rules supplied by user.
# These are preferred for Node_tblastn classification because they directly reflect
# the intended node semantics for this project.
NODE_CLADE_ABBREV: Dict[str, List[str]] = {
    "N0":  ["Athaliana","Ahalleri","Alyrata","Bstricta","Cgrandiflora","Crubella","Cneglecta","Csativa","Camaraamara","Aalpina","Anemorensis","Boleracea","BrapaFPsc","Esalsugineum","Iindigotica","Merraticum","Rsativus","Sparvula","Salba","TarvensevarMN106","Aarabicum","Thassleriana","Tcacao"],
    "N1":  ["Athaliana","Ahalleri","Alyrata","Bstricta","Cgrandiflora","Crubella","Cneglecta","Csativa","Camaraamara","Aalpina","Anemorensis","Boleracea","BrapaFPsc","Esalsugineum","Iindigotica","Merraticum","Rsativus","Sparvula","Salba","TarvensevarMN106","Aarabicum","Thassleriana"],
    "N2":  ["Athaliana","Ahalleri","Alyrata","Bstricta","Cgrandiflora","Crubella","Cneglecta","Csativa","Camaraamara","Aalpina","Anemorensis","Boleracea","BrapaFPsc","Esalsugineum","Iindigotica","Merraticum","Rsativus","Sparvula","Salba","TarvensevarMN106","Aarabicum"],
    "N3":  ["Athaliana","Ahalleri","Alyrata","Bstricta","Cgrandiflora","Crubella","Cneglecta","Csativa","Camaraamara","Aalpina","Anemorensis","Boleracea","BrapaFPsc","Esalsugineum","Iindigotica","Merraticum","Rsativus","Sparvula","Salba","TarvensevarMN106"],
    "N4":  ["Camaraamara","Alyrata","Ahalleri","Athaliana","Cneglecta","Csativa","Crubella","Cgrandiflora","Bstricta"],
    "N5":  ["Rsativus","BrapaFPsc","Boleracea","Salba","Iindigotica","Sparvula","TarvensevarMN106","Esalsugineum","Merraticum","Anemorensis","Aalpina"],
    "N6":  ["Alyrata","Ahalleri","Athaliana","Cneglecta","Csativa","Crubella","Cgrandiflora","Bstricta"],
    "N7":  ["Anemorensis","Aalpina"],
    "N8":  ["Rsativus","BrapaFPsc","Boleracea","Salba","Iindigotica","Sparvula","TarvensevarMN106","Esalsugineum","Merraticum"],
    "N9":  ["Alyrata","Ahalleri","Athaliana"],
    "N10": ["Cneglecta","Csativa","Crubella","Cgrandiflora","Bstricta"],
    "N11": ["Rsativus","BrapaFPsc","Boleracea","Salba","Iindigotica","Sparvula","TarvensevarMN106","Esalsugineum"],
    "N12": ["Alyrata","Ahalleri"],
    "N13": ["Cneglecta","Csativa","Crubella","Cgrandiflora"],
    "N14": ["TarvensevarMN106","Esalsugineum"],
    "N15": ["Rsativus","BrapaFPsc","Boleracea","Salba","Iindigotica","Sparvula"],
    "N16": ["Cneglecta","Csativa"],
    "N17": ["Crubella","Cgrandiflora"],
    "N18": ["Rsativus","BrapaFPsc","Boleracea","Salba","Iindigotica"],
    "N19": ["Rsativus","BrapaFPsc","Boleracea","Salba"],
    "N20": ["Rsativus","BrapaFPsc","Boleracea"],
    "N21": ["BrapaFPsc","Boleracea"],
}
NODE_CLADE_CANON: Dict[str, Set[str]] = {k: {canon_taxon(x) for x in v} for k, v in NODE_CLADE_ABBREV.items()}


def explicit_ancestors_young_to_old(focal_species: str) -> List[str]:
    focal = canon_taxon(focal_species)
    nodes = [n for n, clade in NODE_CLADE_CANON.items() if focal in clade]
    nodes.sort(key=lambda n: (len(NODE_CLADE_CANON[n]), -int(n[1:])))
    return nodes


# -----------------------------
# FASTA utils
# -----------------------------
def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    seq_parts: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                if header is None:
                    continue
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts)


def write_single_fasta(path: Path, header: str, seq: str) -> None:
    with path.open("w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")


# -----------------------------
# GeneID / subject matching
# -----------------------------
def strip_species_suffix(token: str) -> str:
    """Strip a trailing species tag from query IDs when present.

    Handles forms like:
      AL1014U10010.t1_Alyrata
      AT1G03395.1_Arabidopsis_thaliana
      TaMN106.1G006500.1.p__TarvensevarMN106

    Also removes any leftover trailing underscore(s).
    """
    token = token.strip()
    if not token:
        return token

    suffixes = []
    for ab in ABBREV_TO_CANON.keys():
        suffixes.append("_" + ab)
        suffixes.append("__" + ab)
    for canon in SPECIES_COLUMNS:
        suffixes.append("_" + canon)
        suffixes.append("__" + canon)
    suffixes = sorted(set(suffixes), key=len, reverse=True)

    for suf in suffixes:
        if token.endswith(suf):
            token = token[:-len(suf)]
            break

    token = re.sub(r"_+$", "", token)
    return token


def geneid_from_query_header(hdr: str) -> str:
    """Extract the geneID used to match subject filenames."""
    token = hdr.split()[0]
    return strip_species_suffix(token)


def geneid_variants(gene_id: str) -> List[str]:
    """Return candidate identifiers to match subject filenames."""
    variants: List[str] = []
    gene_id = gene_id.strip()
    if not gene_id:
        return variants

    def add(v: str):
        if v and v not in variants:
            variants.append(v)

    # Strip coordinate suffix _coord1_coord2_len (e.g. AT1G05055_1450824_1450132_96 ->
    # AT1G05055) so headers carrying genomic coordinates match subject filenames.
    # No-op on IDs that don't end in three consecutive underscore-digits groups.
    coord_stripped = re.sub(r"_\d+_\d+_\d+$", "", gene_id)
    seeds = [gene_id, strip_species_suffix(gene_id), coord_stripped, strip_species_suffix(coord_stripped)]

    for seed in seeds:
        add(seed)
        s = re.sub(r"_+$", "", seed)
        add(s)
        add(re.sub(r"\.[ptm]$", "", s))
        add(re.sub(r"\.\d+$", "", s))
        add(re.sub(r"\.[ptm]\d+$", "", s))
        v = re.sub(r"\.p$", "", s)
        add(v)
        add(re.sub(r"\.\d+$", "", v))
        add(re.sub(r"\.v\d+(?:\.\d+)*$", "", s))
        add(re.sub(r"(\.[ptm]\d+)\.v\d+(?:\.\d+)*$", "", s))

    # Gene ID / transcript ID interop: filenames sometimes use gene-level IDs (AT1G31420)
    # and sometimes transcript-level IDs (AT1G31420.1, .2). Generate both directions so a
    # bare header also tries .1/.2, and a .N header already tries bare via re.sub above.
    bare_ids = [v for v in list(variants) if v and "_" not in v and not re.search(r"\.\d+$", v)]
    for bare in bare_ids:
        add(bare + ".1")
        add(bare + ".2")

    return variants


def build_subject_file_list(subject_folder: Path) -> List[Path]:
    """Return all subject fasta files in the folder."""
    exts = {".fa", ".fna", ".fasta", ".fas", ".faa", ".aln"}
    files: List[Path] = []
    for p in sorted(subject_folder.iterdir()):
        if p.is_file() and p.suffix.lower() in exts:
            files.append(p)
    return files


def _gene_boundary_pattern(gene_id: str) -> re.Pattern:
    gid = re.escape(gene_id)
    return re.compile(rf"(?<![A-Za-z0-9]){gid}(?![A-Za-z0-9])")


def _disambiguate_by_coords(gene_id: str, candidates: List[Path]) -> Optional[List[Path]]:
    """When coordinate stripping causes multiple matches, use header coords to pick the right file.

    Header format:  GENEID_coord1_coord2_len  (e.g. AT1G31420_11250395_11250433_12)
    Filename format: ..._chr?_start_end_nogaps.fna  (0-based start, 1-based end)

    Returns a single-element list if exactly one file matches the header coords,
    otherwise None (caller falls back to the full candidate list).
    """
    m = re.search(r"_(\d+)_(\d+)_\d+$", gene_id)
    if not m:
        return None
    c1, c2 = int(m.group(1)), int(m.group(2))
    lo, hi = min(c1, c2), max(c1, c2)

    hits = []
    for p in candidates:
        fm = re.search(r"_(\d+)_(\d+)_nogaps", p.name)
        if not fm:
            continue
        f_start, f_end = int(fm.group(1)), int(fm.group(2))
        # Accept 0-based start (f_start + 1 == lo) or 1-based start (f_start == lo)
        if f_end == hi and (f_start + 1 == lo or f_start == lo):
            hits.append(p)
    return hits if len(hits) == 1 else None


def find_subjects_for_gene(gene_id: str, subject_files: List[Path]) -> List[Path]:
    """Find subject FASTA files whose basename contains a boundary-matched gene-id variant."""
    matches: List[Path] = []
    variants = geneid_variants(gene_id)
    for var in variants:
        pat = _gene_boundary_pattern(var)
        cur = [p for p in subject_files if pat.search(p.name)]
        if cur:
            matches.extend(cur)
    seen: Set[Path] = set()
    out: List[Path] = []
    for p in matches:
        if p not in seen:
            seen.add(p)
            out.append(p)

    # If coordinate stripping produced multiple candidates, use header coords to disambiguate.
    if len(out) > 1:
        resolved = _disambiguate_by_coords(gene_id, out)
        if resolved:
            return resolved

    return out


# -----------------------------
# BLAST / derived fields
# -----------------------------
BLAST_FIELDS: List[str] = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "qseq",
    "sseq",
    "sframe",
]

DERIVED_FIELDS: List[str] = [
    "qseq_NoGaps",
    "sseq_NoGaps",
    "ORF_START-STOP",
    "frac_ORFss",
    "1st_PTC_position",
    "frac_CDS_seq",
]

SANITY_FIELDS: List[str] = [
    "longest_NoPTC_segment",
]

EXTRA_FIELDS: List[str] = [
    "Species",
    "geneID_species",
]

OUT_HEADER_BASE: List[str] = ["geneID", "subject_file"] + BLAST_FIELDS + DERIVED_FIELDS + EXTRA_FIELDS


def run(cmd: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def nogaps(seq: str) -> str:
    return seq.replace("-", "")


def replace_stops_with_X(seq: str) -> str:
    return seq.replace("*", "X")


def first_ptc_pos_1based(sseq_nogaps: str) -> int:
    i = sseq_nogaps.find("X")
    return 0 if i < 0 else i + 1


def longest_no_ptc_segment_len(sseq_nogaps: str) -> int:
    if not sseq_nogaps:
        return 0
    # split on X and take max segment length
    return max((len(seg) for seg in sseq_nogaps.split("X")), default=0)


def longest_orf_start_stop_len(sseq_nogaps: str) -> int:
    """Length used for ORF_START-STOP.

    Strict ORF definition:
      - If no M is present, return 0.
      - If at least one M is present, use the stretch from the first M to the
        first downstream X, excluding the X; if no downstream X exists, use the
        stretch from the first M to the end of sseq_nogaps.
    """
    if not sseq_nogaps:
        return 0

    first_m = sseq_nogaps.find("M")
    if first_m < 0:
        return 0

    first_x_after_m = sseq_nogaps.find("X", first_m)
    if first_x_after_m >= 0:
        return first_x_after_m - first_m
    return len(sseq_nogaps) - first_m

def frac_orfss(qlen: int, orf_start_stop_len: int) -> float:
    """frac_ORFss = (ORF_START-STOP / (qlen-1)) * 100."""
    denom = qlen - 1
    if denom <= 0:
        return 0.0
    return (orf_start_stop_len / denom) * 100.0


def frac_cds_seq(qlen: int, sseq_nogaps: str) -> float:
    """
    frac_CDS_seq = (longest stop-free segment length in sseq-NoGaps / (qlen-1)) * 100
    (If no PTCs, this reduces to len(sseq-NoGaps)/(qlen-1)*100.)
    """
    denom = qlen - 1
    if denom <= 0:
        return 0.0
    longest = longest_no_ptc_segment_len(sseq_nogaps)
    return (longest / denom) * 100.0


def species_from_sseqid(sseqid: str) -> str:
    """
    Species label for output #1:
      - take the substring before first '.'
      - if empty, return 'NA'
    Examples:
      Arabidopsis_thaliana.Chr1 -> Arabidopsis_thaliana
      N9.N9refChr165 -> N9
    """
    if not sseqid:
        return "NA"
    pref = sseqid.split(".", 1)[0].strip()
    return pref if pref else "NA"


# -----------------------------
# Newick parsing and tree-derived lookups
# -----------------------------
@dataclass
class Node:
    name: Optional[str]  # taxon name for leaves, node label for internal nodes (e.g., N9), or None
    children: List["Node"]

    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, s: str):
        self.s = s.strip()
        self.i = 0
        self.n = len(self.s)

    def peek(self) -> str:
        return self.s[self.i] if self.i < self.n else ""

    def consume(self, ch: str) -> None:
        if self.peek() != ch:
            raise ValueError(f"Expected '{ch}' at pos {self.i}, got '{self.peek()}'")
        self.i += 1

    def skip_ws(self) -> None:
        while self.i < self.n and self.s[self.i].isspace():
            self.i += 1

    def parse_name(self) -> str:
        """
        Parse a node/leaf label token (until :,),; or whitespace.
        """
        self.skip_ws()
        start = self.i
        while self.i < self.n and self.s[self.i] not in [":", ",", ")", "(", ";"] and not self.s[self.i].isspace():
            self.i += 1
        return self.s[start:self.i].strip()

    def skip_branchlen(self) -> None:
        self.skip_ws()
        if self.peek() == ":":
            self.i += 1
            # consume number (and possible scientific notation)
            while self.i < self.n and self.s[self.i] not in [",", ")", ";"]:
                if self.s[self.i] == "[":
                    break
                self.i += 1
            self.skip_ws()

    def parse_subtree(self) -> Node:
        self.skip_ws()
        if self.peek() == "(":
            # internal node
            self.consume("(")
            children: List[Node] = []
            while True:
                child = self.parse_subtree()
                children.append(child)
                self.skip_ws()
                if self.peek() == ",":
                    self.consume(",")
                    continue
                elif self.peek() == ")":
                    self.consume(")")
                    break
                else:
                    raise ValueError(f"Unexpected char '{self.peek()}' in internal node at pos {self.i}")
            # optional internal label
            label = self.parse_name()
            self.skip_branchlen()
            # ignore support numbers etc; keep only N\d+ labels if present, else None
            label = label if re.fullmatch(r"N\d+", label) else (label if label else None)
            return Node(name=label, children=children)
        else:
            # leaf
            label = self.parse_name()
            if not label:
                raise ValueError(f"Missing leaf label at pos {self.i}")
            self.skip_branchlen()
            return Node(name=label, children=[])

    def parse(self) -> Node:
        node = self.parse_subtree()
        self.skip_ws()
        if self.peek() == ";":
            self.consume(";")
        self.skip_ws()
        if self.i != self.n:
            # tolerate trailing whitespace
            rest = self.s[self.i:].strip()
            if rest:
                raise ValueError(f"Unexpected trailing content after Newick: {rest[:50]}")
        return node


def parse_newick_file(path: Path) -> Node:
    s = path.read_text().strip()
    return NewickParser(s).parse()


def build_node_to_tips(root: Node) -> Dict[str, Set[str]]:
    """
    Map internal node label (e.g., N9) -> set of canonical tip taxa under that node.
    Only includes nodes whose name matches N\\d+.
    """
    node_to_tips: Dict[str, Set[str]] = {}

    def dfs(n: Node) -> Set[str]:
        if n.is_leaf():
            return {canon_taxon(n.name or "")}
        tips: Set[str] = set()
        for c in n.children:
            tips |= dfs(c)
        if n.name and re.fullmatch(r"N\d+", n.name):
            node_to_tips[n.name] = set(tips)
        return tips

    dfs(root)
    return node_to_tips


def build_species_to_ancestors(root: Node) -> Dict[str, List[str]]:
    """
    Map canonical tip taxon -> list of ancestor node labels (N\\d+) from youngest -> oldest.
    """
    species_to_anc: Dict[str, List[str]] = {}

    def walk(n: Node, anc_stack: List[str]) -> None:
        # push this node label if it's N\d+
        next_stack = anc_stack
        if n.name and re.fullmatch(r"N\d+", n.name):
            next_stack = anc_stack + [n.name]
        if n.is_leaf():
            sp = canon_taxon(n.name or "")
            # anc_stack currently oldest->youngest? We appended as we go down, so it is root->... order.
            # We want youngest->oldest, so reverse.
            nlist = [x for x in next_stack if re.fullmatch(r"N\d+", x)]
            species_to_anc[sp] = list(reversed(nlist))
            return
        for c in n.children:
            walk(c, next_stack)

    walk(root, [])
    return species_to_anc


def classify_node_tblastn_tree_based(
    pres: Dict[str, int],
    focal_species: str,
    node_to_tips: Dict[str, Set[str]],
    ancestors_young_to_old: List[str],
) -> str:
    """Classify Node_tblastn from the explicit N0-N21 node definitions.

    Rule used:
      - positives = canonical species with value > 0
      - if no positives: NA
      - if only focal positive: <focal>-specific
      - otherwise choose the youngest explicit node on the focal lineage whose
        clade contains all positives

    The explicit user-supplied node memberships are treated as the source of truth.
    The parsed Newick is kept only as a fallback if a focal species somehow lacks
    an explicit lineage.
    """
    focal = canon_taxon(focal_species)
    positives: Set[str] = {sp for sp, v in pres.items() if int(v) > 0 and sp in CANON_SET}

    if not positives:
        return "NA"

    if positives == {focal}:
        if focal == "Arabidopsis_thaliana":
            return "Athaliana-specific"
        return f"{focal}-specific"

    explicit_anc = explicit_ancestors_young_to_old(focal)
    if explicit_anc:
        for node_label in explicit_anc:
            clade = NODE_CLADE_CANON[node_label]
            if positives.issubset(clade):
                return node_label
        # This should not normally happen because N0 contains all species, but
        # keep a defensive fallback.
        return explicit_anc[-1]

    for node_label in ancestors_young_to_old:
        clade = node_to_tips.get(node_label, set())
        if positives.issubset(clade):
            return node_label

    return ancestors_young_to_old[-1] if ancestors_young_to_old else "NA"


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    ap = argparse.ArgumentParser()

    ap.add_argument("query_proteins_faa", help="Protein multi-FASTA query")
    ap.add_argument("subject_folder", help="Folder containing subject DNA FASTA files")
    ap.add_argument("output_tsv", help="Output TSV of per-HSP hits")
    ap.add_argument("summary_output_tsv", help="Output TSV summary by geneID")
    ap.add_argument("ancseq_output_tsv", help="Output TSV summary of ancestral-sequence ORF presence by node")

    ap.add_argument("--tblastn", default="tblastn", help="Path to tblastn executable")
    ap.add_argument("--makeblastdb", default="makeblastdb", help="Path to makeblastdb executable")
    ap.add_argument("--first-match", action="store_true", help="If multiple subject files match geneID, use the first match.")
    ap.add_argument("--evalue", type=float, default=1.0)
    ap.add_argument("--word_size", type=int, default=2)
    ap.add_argument("--matrix", default="BLOSUM80")
    ap.add_argument("--seg", choices=["yes", "no"], default="no")
    ap.add_argument("--comp_based_stats", type=int, default=0)
    ap.add_argument("--max_hsps", type=int, default=10)
    ap.add_argument("--max_target_seqs", type=int, default=50)
    ap.add_argument(
        "--continuous",
        action="store_true",
        help="Encourage a single gapped HSP per subject: sets max_hsps=1 and adds -xdrop_gap 50 -xdrop_ungap 50.",
    )

    ap.add_argument("--sanity-check", action="store_true", help="Append longest_NoPTC_segment column to output #1.")

    # summary filters
    ap.add_argument("--frac_ORFss", type=float, default=80.0, help="Minimum frac_ORFss threshold (default 80).")
    ap.add_argument("--pident", type=float, default=40.0, help="Minimum pident threshold (default 40).")

    # tree-based node classification
    ap.add_argument("--tree", type=str, default=None, help="Newick tree file with node labels (N#).")
    ap.add_argument("--focal-species", type=str, default=None, help="Focal species (canonical or abbrev).")
    ap.add_argument("--debug-nodes", action="store_true", help="Print node-debug info per geneID to stderr.")
    ap.add_argument("--verbose", action="store_true", help="Print progress stats (queries, matched subjects, HSPs) to stderr.")

    args = ap.parse_args()

    query_path = Path(args.query_proteins_faa)
    subject_folder = Path(args.subject_folder)
    out_path = Path(args.output_tsv)
    summary_path = Path(args.summary_output_tsv)
    ancseq_path = Path(args.ancseq_output_tsv)

    tblastn_exe = args.tblastn
    makeblastdb_exe = args.makeblastdb

    # Parse tree if provided
    node_to_tips: Dict[str, Set[str]] = {}
    species_to_ancestors: Dict[str, List[str]] = {}
    focal_canon: Optional[str] = None
    if args.tree and args.focal_species:
        root = parse_newick_file(Path(args.tree))
        node_to_tips = build_node_to_tips(root)
        species_to_ancestors = build_species_to_ancestors(root)
        try:
            focal_canon = resolve_species_name(args.focal_species, SPECIES_COLUMNS)
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(2)
        if focal_canon not in species_to_ancestors and not explicit_ancestors_young_to_old(focal_canon):
            print(
                f"ERROR: focal species '{args.focal_species}' resolved to '{focal_canon}' but was not found in tree tips or explicit nodes.",
                file=sys.stderr,
            )
            sys.exit(2)
        if args.verbose:
            print(f"[FOCAL] input={args.focal_species} resolved={focal_canon}", file=sys.stderr)

    out_header = list(OUT_HEADER_BASE)
    if args.sanity_check:
        # insert sanity fields right before EXTRA_FIELDS? user asked at end; keep at end for simplicity
        out_header = ["geneID", "subject_file"] + BLAST_FIELDS + DERIVED_FIELDS + (SANITY_FIELDS) + EXTRA_FIELDS

    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    ancseq_path.parent.mkdir(parents=True, exist_ok=True)

    # Index subjects
    subject_files = build_subject_file_list(subject_folder)

    outfmt = "6 " + " ".join(BLAST_FIELDS)
    max_hsps_effective = 1 if args.continuous else args.max_hsps

    # Collect for summary: geneID -> species -> count
    n_queries = 0
    n_matched = 0
    n_hsps = 0
    n_hsps_passing = 0
    summary_counts: Dict[str, Dict[str, int]] = {}
    ancseq_counts: Dict[str, Dict[str, int]] = {}
    all_query_gene_ids: List[str] = []

    # Write output #1
    with out_path.open("w") as out:
        out.write("\t".join(out_header) + "\n")

        for qhdr, qseq in iter_fasta(query_path):
            n_queries += 1
            gene_id = geneid_from_query_header(qhdr)
            all_query_gene_ids.append(gene_id)
            matches = find_subjects_for_gene(gene_id, subject_files)
            if args.verbose and not matches:
                print(f"[NO MATCH] gene_id={gene_id} variants={geneid_variants(gene_id)}", file=sys.stderr)
            if args.verbose and matches:
                print(f"[MATCH] gene_id={gene_id} matched={len(matches)} first={matches[0].name}", file=sys.stderr)

            if not matches:
                na_tail_len = len(out_header) - 2
                out.write("\t".join([gene_id, "NA"] + ["NA"] * na_tail_len) + "\n")
                continue

            if len(matches) > 1 and not args.first_match:
                print(f"ERROR: multiple subject files match {gene_id}. Use --first-match.", file=sys.stderr)
                for m in matches:
                    print(f"  - {m}", file=sys.stderr)
                sys.exit(1)

            subject_file = matches[0]
            n_matched += 1

            with tempfile.TemporaryDirectory(prefix=f"tblastn_{gene_id}_") as tmpdir:
                tmpdir_p = Path(tmpdir)
                qfile = tmpdir_p / f"{gene_id}.faa"
                db_prefix = tmpdir_p / "subject_db"

                write_single_fasta(qfile, qhdr.split()[0], qseq)

                cmd_db = [makeblastdb_exe, "-in", str(subject_file), "-dbtype", "nucl", "-out", str(db_prefix)]
                p_db = run(cmd_db)
                if p_db.returncode != 0:
                    na_tail_len = len(out_header) - 2
                    out.write("\t".join([gene_id, subject_file.name] + ["NA"] * na_tail_len) + "\n")
                    continue

                cmd_blast = [
                    tblastn_exe,
                    "-query", str(qfile),
                    "-db", str(db_prefix),
                    "-evalue", str(args.evalue),
                    "-word_size", str(args.word_size),
                    "-matrix", str(args.matrix),
                    "-seg", str(args.seg),
                    "-comp_based_stats", str(args.comp_based_stats),
                    "-max_hsps", str(max_hsps_effective),
                    "-max_target_seqs", str(args.max_target_seqs),
                    "-outfmt", outfmt,
                ]
                if args.continuous:
                    cmd_blast += ["-xdrop_gap", "50", "-xdrop_ungap", "50"]

                p_bl = run(cmd_blast)
                if p_bl.returncode != 0:
                    na_tail_len = len(out_header) - 2
                    out.write("\t".join([gene_id, subject_file.name] + ["NA"] * na_tail_len) + "\n")
                    continue

                txt = p_bl.stdout.strip()
                if not txt:
                    na_tail_len = len(out_header) - 2
                    out.write("\t".join([gene_id, subject_file.name] + ["NA"] * na_tail_len) + "\n")
                    continue

                for line in txt.splitlines():
                    n_hsps += 1
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) != len(BLAST_FIELDS):
                        continue

                    row = dict(zip(BLAST_FIELDS, parts))

                    try:
                        qlen = int(float(row["qlen"]))
                    except Exception:
                        qlen = 0

                    # Replace '*' -> 'X' in sseq
                    row["sseq"] = replace_stops_with_X(row["sseq"])

                    qseq_ng = nogaps(row["qseq"])
                    sseq_ng = nogaps(row["sseq"])

                    ptc_pos = first_ptc_pos_1based(sseq_ng)
                    frac = frac_cds_seq(qlen, sseq_ng)
                    longest_seg = longest_no_ptc_segment_len(sseq_ng)
                    orf_start_stop = longest_orf_start_stop_len(sseq_ng)
                    frac_orf = frac_orfss(qlen, orf_start_stop)

                    # Species per-HSP from sseqid prefix
                    species_raw = species_from_sseqid(row["sseqid"])
                    species = canon_taxon(species_raw)  # map abbrev to canonical if needed
                    geneid_species = f"{gene_id}_{species_raw}"  # keep raw label in this composite for traceability

                    derived_vals = [qseq_ng, sseq_ng, str(orf_start_stop), f"{frac_orf:.6f}", str(ptc_pos), f"{frac:.6f}"]
                    sanity_vals: List[str] = [str(longest_seg)] if args.sanity_check else []
                    extra_vals = [species_raw, geneid_species]

                    out.write(
                        "\t".join(
                            [gene_id, subject_file.name]
                            + [row[f] for f in BLAST_FIELDS]
                            + derived_vals
                            + sanity_vals
                            + extra_vals
                        )
                        + "\n"
                    )

                    # Summary counting: only if species maps to a canonical tip
                    # Apply thresholds on this HSP (frac_CDS_seq and pident)
                    try:
                        pident_val = float(row["pident"])
                    except Exception:
                        pident_val = float("nan")

                    if (not math.isnan(pident_val)) and (frac_orf >= args.frac_ORFss) and (pident_val >= args.pident):
                        if species in CANON_SET:
                            summary_counts.setdefault(gene_id, {}).setdefault(species, 0)
                            summary_counts[gene_id][species] += 1
                            n_hsps_passing += 1
                        if re.fullmatch(r"N\d+", species_raw) and orf_start_stop > 0:
                            ancseq_counts.setdefault(gene_id, {}).setdefault(species_raw, 0)
                            ancseq_counts[gene_id][species_raw] = 1

    if n_queries == 0:
        print("WARNING: No query sequences were read from the query FASTA. Output files will contain only headers.", file=sys.stderr)
    if args.verbose:
        n_summary_rows = len(dict.fromkeys(all_query_gene_ids))
        print(f"[STATS] queries={n_queries} matched_subject_files={n_matched} hsps_parsed={n_hsps} hsps_passing_summary_filters={n_hsps_passing} summary_rows={n_summary_rows}", file=sys.stderr)

    # Write summary output #2
    with summary_path.open("w") as sf:
        sf.write("\t".join(["geneID", "Node_tblastn"] + SPECIES_COLUMNS) + "\n")

        for gene_id in sorted(dict.fromkeys(all_query_gene_ids).keys()):
            counts = summary_counts.get(gene_id, {})
            pres = {sp: (1 if counts.get(sp, 0) > 0 else 0) for sp in SPECIES_COLUMNS}

            node = "NA"
            if focal_canon and node_to_tips and species_to_ancestors:
                node = classify_node_tblastn_tree_based(
                    pres=pres,
                    focal_species=focal_canon,
                    node_to_tips=node_to_tips,
                    ancestors_young_to_old=species_to_ancestors.get(focal_canon, []),
                )

                if args.debug_nodes:
                    positives = sorted([sp for sp, v in pres.items() if v == 1])
                    anc = explicit_ancestors_young_to_old(focal_canon) or species_to_ancestors.get(focal_canon, [])
                    print(f"[DEBUG] geneID={gene_id} focal={focal_canon}", file=sys.stderr)
                    print(f"[DEBUG] positives={','.join(positives) if positives else 'NONE'}", file=sys.stderr)
                    print(f"[DEBUG] ancestors_young_to_old={'|'.join(anc) if anc else 'NONE'}", file=sys.stderr)
                    for nl in anc:
                        cl = NODE_CLADE_CANON.get(nl, node_to_tips.get(nl, set()))
                        outside = sorted(set(positives) - set(cl))
                        print(f"[DEBUG]   node={nl} outside={','.join(outside) if outside else 'NONE'}", file=sys.stderr)
                    print(f"[DEBUG] chosen={node}", file=sys.stderr)

            row = [gene_id, node] + [str(counts.get(sp, 0)) for sp in SPECIES_COLUMNS]
            sf.write("\t".join(row) + "\n")

    # Write AncSeq output #3
    ancseq_nodes = [f"N{i}" for i in range(21, -1, -1)]
    with ancseq_path.open("w") as af:
        af.write("\t".join(["geneID", "Node_AncSeq_ORF"] + ancseq_nodes) + "\n")

        for gene_id in sorted(dict.fromkeys(all_query_gene_ids).keys()):
            counts = ancseq_counts.get(gene_id, {})
            pres = {node: (1 if counts.get(node, 0) > 0 else 0) for node in ancseq_nodes}

            focal_specific_label = f"{focal_canon}-specific" if focal_canon else "NA"
            focal_lineage_nodes = explicit_ancestors_young_to_old(focal_canon) if focal_canon else []
            if not focal_lineage_nodes and focal_canon:
                focal_lineage_nodes = species_to_ancestors.get(focal_canon, [])

            node_anc = focal_specific_label
            for node in reversed(focal_lineage_nodes):  # oldest -> youngest, focal lineage only
                if pres.get(node, 0) > 0:
                    node_anc = node
                    break

            row = [gene_id, node_anc] + [str(pres[node]) for node in ancseq_nodes]
            af.write("\t".join(row) + "\n")

    print("Done.")
    print(f"Output TSV        : {out_path}")
    print(f"Summary output TSV: {summary_path}")
    print(f"AncSeq output TSV : {ancseq_path}")


if __name__ == "__main__":
    main()
