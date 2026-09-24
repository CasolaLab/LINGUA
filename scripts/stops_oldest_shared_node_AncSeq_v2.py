#!/usr/bin/env python3

"""
Usage:
python3 stops_oldest_shared_node_AncSeq_v2.py stops_dir [outdir]

Description:
    Batch-process all .csv and .tsv files in an input directory.

    For each input file, the script:

      1. Builds stops_by_species_<FocalSpecies>.tsv
         - one row per (gff_transcriptID, StopCodonStart_RefPos)
         - counts duplicate occurrences within a species if they occur on
           different chromosomes/scaffolds
         - ancestral N0...N21 rows are excluded from this output

      2. Builds nodes_stops_summary_<FocalSpecies>.tsv
         - one row per (gff_transcriptID, StopCodonStart_RefPos)
         - counts UNIQUE extant species only
         - summarizes species presence into focal-lineage node bins
         - reports the oldest shared stop node based on extant species support

      3. Builds stops_by_AncSeq_<FocalSpecies>.tsv
         - one row per (gff_transcriptID, StopCodonStart_RefPos)
         - counts ancestral sequence rows N0...N21
         - columns are N21, N20, ..., N0

      4. Builds all_stops_nodes_<FocalSpecies>.tsv
         - one row per geneID
         - Nodes_stops contains all ancestral nodes on the focal species lineage
           where at least one STOP codon is found for that gene
         - nodes are separated by vertical bars, e.g. N1|N2|N3|N4|N6
         - if no focal-lineage ancestral node has a STOP codon, prints na

Notes:
    - Species names are parsed from the 'Species' column as the text before
      the first dot.
    - Extant species rows and ancestral N0...N21 rows are separated after
      parsing Species.
    - The focal species is inferred from the file_name column.
    - The script can read both .csv and .tsv input files.
    - Missing node calls are printed as lowercase 'na'.
"""

import os
import re
import sys
import glob
import pandas as pd

ALIASES = {
    "Athaliana": "Arabidopsis_thaliana",
    "Ahalleri": "Arabidopsis_halleri",
    "Alyrata": "Arabidopsis_lyrata",
    "Bstricta": "Boechera_stricta",
    "Cgrandiflora": "Capsella_grandiflora",
    "Crubella": "Capsella_rubella",
    "Cneglecta": "Camelina_neglecta",
    "Csativa": "Camelina_sativa",
    "Camaraamara": "Cardamine_amara_amara",
    "Cardamine_amara": "Cardamine_amara_amara",
    "Aalpina": "Arabis_alpina",
    "Anemorensis": "Arabis_nemorensis",
    "Boleracea": "Brassica_oleracea",
    "BrapaFPsc": "Brassica_rapa_FPsc",
    "Brassica_rapa": "Brassica_rapa_FPsc",
    "Esalsugineum": "Eutrema_salsugineum",
    "Iindigotica": "Isatis_indigotica",
    "Merraticum": "Microthlaspi_erraticum",
    "Rsativus": "Raphanus_sativus",
    "Sparvula": "Schrenkiella_parvula",
    "Salba": "Sinapis_alba",
    "TarvensevarMN106": "Thlaspi_arvense_MN106",
    "Thlaspi_arvense": "Thlaspi_arvense_MN106",
    "Aarabicum": "Aethionema_arabicum",
    "Thassleriana": "Tarenaya_hassleriana",
    "Tcacao": "Theobroma_cacao",
}

def canon_taxon(x: str) -> str:
    return ALIASES.get(x, x)

SPECIES_ORDER = [
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

ANCSEQ_NODE_ORDER = [f"N{i}" for i in range(21, -1, -1)]

NODE_CLADE_ABBREV = {
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

NODE_CLADE_CANON = {
    node: {canon_taxon(x) for x in members}
    for node, members in NODE_CLADE_ABBREV.items()
}

def parse_species_field(s: str) -> str:
    if pd.isna(s):
        return s
    return canon_taxon(str(s).split(".", 1)[0])

def is_ancestral_name(s: str) -> bool:
    if pd.isna(s):
        return False
    return bool(re.fullmatch(r"N\d+", str(s)))

def detect_focal_species(file_name: str) -> str:
    matches = [sp for sp in SPECIES_ORDER if sp in str(file_name)]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        return sorted(matches, key=len, reverse=True)[0]
    return None

def read_table_auto(infile: str) -> pd.DataFrame:
    if infile.lower().endswith(".csv"):
        return pd.read_csv(infile, sep=",", dtype=str)
    elif infile.lower().endswith(".tsv"):
        return pd.read_csv(infile, sep="\t", dtype=str)
    else:
        return pd.read_csv(infile, sep=None, engine="python", dtype=str)

def focal_lineage_nodes(focal: str):
    nodes = [node for node, members in NODE_CLADE_CANON.items() if focal in members]
    return sorted(nodes, key=lambda x: int(x[1:]))

def build_lineage_bins(lineage_nodes):
    bins = {}
    for i, node in enumerate(lineage_nodes):
        if i == len(lineage_nodes) - 1:
            bins[node] = set(NODE_CLADE_CANON[node])
        else:
            younger = lineage_nodes[i + 1]
            bins[node] = set(NODE_CLADE_CANON[node]) - set(NODE_CLADE_CANON[younger])
    return bins

def assign_species_to_lineage_bin(species: str, lineage_nodes, lineage_bins):
    for node in reversed(lineage_nodes):
        if species in lineage_bins[node]:
            return node
    return None

def oldest_shared_node_from_counts(node_counts, lineage_nodes):
    present = [node for node in lineage_nodes if node_counts.get(node, 0) > 0]
    if len(present) < 2:
        return "na"
    return present[0]

def safe_positive_count(value) -> bool:
    try:
        return int(value) > 0
    except (TypeError, ValueError):
        return False

def oldest_and_younger_from_ancseq_row(row, lineage_nodes):
    present_nodes = []
    for node in lineage_nodes:
        if node in row.index and safe_positive_count(row[node]):
            present_nodes.append(node)

    if len(present_nodes) < 2:
        return "na", "na"

    return present_nodes[0], present_nodes[1]


def collect_gene_lineage_stop_nodes(ancseq_counts: pd.DataFrame, lineage_nodes):
    """
    Build rows for all_stops_nodes_<FocalSpecies>.tsv.

    New rules:
        1. Exclude N0 and N1 before evaluating nodes.
        2. Report every STOP position of a gene that is present in at least two
           remaining focal-lineage ancestral nodes.
        3. Add StopCodonStart_RefPos as the second output column.

    Output columns:
        geneID  StopCodonStart_RefPos  Nodes_stops

    Nodes are printed oldest -> youngest according to lineage_nodes, after
    excluding N0 and N1.
    """
    rows = []

    lineage_nodes_no_N0_N1 = [
        node for node in lineage_nodes
        if node not in {"N0", "N1"}
    ]

    for _, row in ancseq_counts.iterrows():
        geneID = row["gff_transcriptID"]
        stop_start = row["StopCodonStart_RefPos"]

        present_nodes = []

        for node in lineage_nodes_no_N0_N1:
            if node not in row.index:
                continue

            value = pd.to_numeric(row[node], errors="coerce")
            if pd.notna(value) and value > 0:
                present_nodes.append(node)

        if len(present_nodes) < 2:
            continue

        rows.append({
            "geneID": geneID,
            "StopCodonStart_RefPos": stop_start,
            "Nodes_stops": "|".join(present_nodes),
        })

    if not rows:
        return pd.DataFrame(columns=["geneID", "StopCodonStart_RefPos", "Nodes_stops"])

    return pd.DataFrame(rows, columns=["geneID", "StopCodonStart_RefPos", "Nodes_stops"])

def process_one_file(infile: str, outdir: str):
    df = read_table_auto(infile)

    required = {"file_name", "gff_transcriptID", "Species", "StopCodonStart_RefPos"}
    missing = required - set(df.columns)
    if missing:
        print(f"WARNING: skipping {os.path.basename(infile)}; missing columns: {', '.join(sorted(missing))}")
        return

    df = df.copy()
    df["StopCodonStart_RefPos"] = df["StopCodonStart_RefPos"].astype(str)
    df["Species_clean"] = df["Species"].apply(parse_species_field)

    df_species = df[~df["Species_clean"].apply(is_ancestral_name)].copy()
    df_ancseq = df[
        df["Species_clean"].apply(is_ancestral_name) &
        df["Species_clean"].isin(ANCSEQ_NODE_ORDER)
    ].copy()

    if df_species.empty and df_ancseq.empty:
        print(f"WARNING: skipping {os.path.basename(infile)}; no usable species or ancestral rows found")
        return

    focal = detect_focal_species(df["file_name"].iloc[0])
    if focal is None:
        print(f"WARNING: skipping {os.path.basename(infile)}; could not detect focal species from file_name")
        return

    lineage_nodes = focal_lineage_nodes(focal)
    lineage_bins = build_lineage_bins(lineage_nodes)

    if not df_species.empty:
        species_counts = (
            df_species.groupby(["gff_transcriptID", "StopCodonStart_RefPos", "Species_clean"])
            .size()
            .unstack(fill_value=0)
        )

        for sp in SPECIES_ORDER:
            if sp not in species_counts.columns:
                species_counts[sp] = 0

        species_counts = species_counts[SPECIES_ORDER].reset_index()
    else:
        species_counts = df[["gff_transcriptID", "StopCodonStart_RefPos"]].drop_duplicates().copy()
        for sp in SPECIES_ORDER:
            species_counts[sp] = 0

    out1 = os.path.join(outdir, f"stops_by_species_{focal}.tsv")
    species_counts.to_csv(out1, sep="\t", index=False)

    if not df_species.empty:
        unique_presence = (
            df_species[["gff_transcriptID", "StopCodonStart_RefPos", "Species_clean"]]
            .drop_duplicates()
            .copy()
        )

        unique_presence["AssignedNode"] = unique_presence["Species_clean"].apply(
            lambda sp: assign_species_to_lineage_bin(sp, lineage_nodes, lineage_bins)
        )

        node_summary = (
            unique_presence
            .dropna(subset=["AssignedNode"])
            .groupby(["gff_transcriptID", "StopCodonStart_RefPos", "AssignedNode"])
            .size()
            .unstack(fill_value=0)
        )

        for node in lineage_nodes:
            if node not in node_summary.columns:
                node_summary[node] = 0

        node_summary = node_summary[lineage_nodes].reset_index()
    else:
        node_summary = df[["gff_transcriptID", "StopCodonStart_RefPos"]].drop_duplicates().copy()
        for node in lineage_nodes:
            node_summary[node] = 0

    oldest_nodes_per_stop = []
    for _, row in node_summary.iterrows():
        node_counts = {node: int(row[node]) for node in lineage_nodes}
        oldest_nodes_per_stop.append(oldest_shared_node_from_counts(node_counts, lineage_nodes))

    output_node_cols = list(reversed(lineage_nodes))
    rename_map = {node: f"{node}_#species_stop" for node in output_node_cols}

    node_summary_out = node_summary[["gff_transcriptID", "StopCodonStart_RefPos"] + output_node_cols].copy()
    node_summary_out.rename(columns=rename_map, inplace=True)
    node_summary_out["Oldest_shared_stop_node"] = oldest_nodes_per_stop

    out2 = os.path.join(outdir, f"nodes_stops_summary_{focal}.tsv")
    node_summary_out.to_csv(out2, sep="\t", index=False)

    if not df_ancseq.empty:
        ancseq_counts = (
            df_ancseq.groupby(["gff_transcriptID", "StopCodonStart_RefPos", "Species_clean"])
            .size()
            .unstack(fill_value=0)
        )

        for node in ANCSEQ_NODE_ORDER:
            if node not in ancseq_counts.columns:
                ancseq_counts[node] = 0

        ancseq_counts = ancseq_counts[ANCSEQ_NODE_ORDER].reset_index()
    else:
        ancseq_counts = df[["gff_transcriptID", "StopCodonStart_RefPos"]].drop_duplicates().copy()
        for node in ANCSEQ_NODE_ORDER:
            ancseq_counts[node] = 0

    out4 = os.path.join(outdir, f"stops_by_AncSeq_{focal}.tsv")
    ancseq_counts.to_csv(out4, sep="\t", index=False)

    # -----------------------------
    # Output 4: all focal-lineage ancestral nodes with STOP codons per gene
    # -----------------------------
    all_nodes_df = collect_gene_lineage_stop_nodes(ancseq_counts, lineage_nodes)

    out3 = os.path.join(outdir, f"all_stops_nodes_{focal}.tsv")
    all_nodes_df.to_csv(out3, sep="\t", index=False)

    print(f"Done: {os.path.basename(infile)}")
    print(f"  Focal species: {focal}")
    print(f"  Wrote: {os.path.basename(out1)}")
    print(f"  Wrote: {os.path.basename(out2)}")
    print(f"  Wrote: {os.path.basename(out4)}")
    print(f"  Wrote: {os.path.basename(out3)}")

def main():
    if len(sys.argv) not in (2, 3):
        print(__doc__)
        sys.exit(1)

    indir = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) == 3 else indir

    if not os.path.isdir(indir):
        sys.exit(f"ERROR: input directory not found: {indir}")

    os.makedirs(outdir, exist_ok=True)

    files = sorted(
        f for f in (
            glob.glob(os.path.join(indir, "*.csv")) +
            glob.glob(os.path.join(indir, "*.tsv"))
        )
        if not any(x in os.path.basename(f) for x in [
            "stops_by_species",
            "nodes_stops_summary",
            "stops_by_AncSeq",
            "oldest_stop_node",
            "all_stops_nodes"
        ])
    )

    if not files:
        sys.exit("ERROR: no valid input files found in input directory")

    for infile in files:
        process_one_file(infile, outdir)

if __name__ == "__main__":
    main()
