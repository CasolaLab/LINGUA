#!/usr/bin/env python3

"""
python nodes_by_coordinates.py \
    coordinates_alignments_all_species \
    nodes_summary_by_species \
    ../cdngs_tblastn_results_Apr-2026 \
    parsed_coordinates

nodes_by_coordinates.py

Parse coordinate alignment files, merge with node summaries and tblastn presence/absence,
and output filtered coordinate tables per focal species.

Input folders:
1. coordinates_alignments_all_species
   Files like:
   Athaliana_DNG_45k_cds_all_species_summary.tsv

2. nodes_summary_by_species
   Files like:
   athaliana_nodes_summary_filtered.tsv

3. cdngs_tblastn_results_Apr-2026
   Files like:
   athaliana_cdngs_tblastn_1078_summary.tsv

Output:
For each focal species, writes:
<species_prefix>_coords_kept.tsv

Usage:
python nodes_by_coordinates.py \
    coordinates_alignments_all_species \
    nodes_summary_by_species \
    ../cdngs_tblastn_results_Apr-2026 \
    parsed_coordinates

Notes:
- The numeric value in *_cdngs_tblastn_<NUMBER>_summary.tsv can vary.
- Some species use different prefixes in different input folders.
- Rows without ORF_Node are removed from final output.
- Empty Enablers_Nodes values are converted to NA.
"""

import os
import sys
import re
import glob
import pandas as pd


# -----------------------------
# Species aliases
# -----------------------------
def canonicalize_species(s: str) -> str:
    if s is None:
        return None

    aliases = {
        # original / title-case short prefixes
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
        "Tarvense_MN106": "Thlaspi_arvense_MN106",
        "TarvensevarMN106": "Thlaspi_arvense_MN106",
        "Aalpina": "Arabis_alpina",
        "Anemorensis": "Arabis_nemorensis",
        "Boleracea": "Brassica_oleracea",
        "Brapa": "Brassica_rapa_FPsc",
        "BrapaFPsc": "Brassica_rapa_FPsc",
        "Brapa_FPsc": "Brassica_rapa_FPsc",
        "Esalsugineum": "Eutrema_salsugineum",
        "Iindigotica": "Isatis_indigotica",
        "Merraticum": "Microthlaspi_erraticum",
        "Rsativus": "Raphanus_sativus",
        "Rsativus_GCF": "Raphanus_sativus",
        "Sparvula": "Schrenkiella_parvula",
        "Salba": "Sinapis_alba",
        "Aarabicum": "Aethionema_arabicum",
        "Thassleriana": "Tarenaya_hassleriana",
        "Tcacao": "Theobroma_cacao",

        # lowercase / filename prefixes
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
        "tarvense_mn106": "Thlaspi_arvense_MN106",
        "tarvensevarmn106": "Thlaspi_arvense_MN106",
        "aalpina": "Arabis_alpina",
        "anemorensis": "Arabis_nemorensis",
        "boleracea": "Brassica_oleracea",
        "brapa": "Brassica_rapa_FPsc",
        "brapafpsc": "Brassica_rapa_FPsc",
        "brapa_fpsc": "Brassica_rapa_FPsc",
        "esalsugineum": "Eutrema_salsugineum",
        "iindigotica": "Isatis_indigotica",
        "merraticum": "Microthlaspi_erraticum",
        "rsativus": "Raphanus_sativus",
        "rsativus_gcf": "Raphanus_sativus",
        "sparvula": "Schrenkiella_parvula",
        "salba": "Sinapis_alba",
        "aarabicum": "Aethionema_arabicum",
        "thassleriana": "Tarenaya_hassleriana",
        "tcacao": "Theobroma_cacao",
    }

    return aliases.get(s, aliases.get(str(s).lower(), s))


# -----------------------------
# Standard output/file prefix
# -----------------------------
def species_to_file_prefix(species):
    """
    Default short species prefix.
    Used for output unless folder-specific prefixes are needed.
    """

    prefix_map = {
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
        "Thlaspi_arvense_MN106": "tarvense",
        "Aethionema_arabicum": "aarabicum",
        "Tarenaya_hassleriana": "thassleriana",
        "Theobroma_cacao": "tcacao",
    }

    return prefix_map.get(species, str(species).lower())


# -----------------------------
# Folder-specific prefixes
# -----------------------------
def species_to_nodes_prefix(species):
    """
    Prefix used in nodes_summary_by_species.

    Special cases:
    - Thlaspi_arvense_MN106 nodes file is tarvensevarmn106_nodes_summary_filtered.tsv
    - Brassica_rapa_FPsc nodes file is brapa_nodes_summary_filtered.tsv
    - Raphanus_sativus nodes file is rsativus_nodes_summary_filtered.tsv
    """

    prefix_map = {
        "Thlaspi_arvense_MN106": "tarvensevarmn106",
        "Brassica_rapa_FPsc": "brapa",
        "Raphanus_sativus": "rsativus",
    }

    return prefix_map.get(species, species_to_file_prefix(species))


def species_to_tblastn_prefix(species):
    """
    Prefix used in cdngs_tblastn_results_Apr-2026.

    Special cases:
    - Thlaspi_arvense_MN106 tblastn file is tarvense_cdngs_tblastn_*_summary.tsv
    - Brassica_rapa_FPsc tblastn file is brapa_cdngs_tblastn_*_summary.tsv
    - Raphanus_sativus tblastn file is rsativus_cdngs_tblastn_*_summary.tsv
    """

    prefix_map = {
        "Thlaspi_arvense_MN106": "tarvense",
        "Brassica_rapa_FPsc": "brapa",
        "Raphanus_sativus": "rsativus",
    }

    return prefix_map.get(species, species_to_file_prefix(species))


# -----------------------------
# Extract geneID from coordinate file_name
# -----------------------------
def extract_geneID(filename):
    """
    Extract geneID from file_name values like:

    Arabidopsis_thaliana_AT1G01335.1_chr1_130735_130858.aln
    Tarenaya_hassleriana_XP_010545484.1_chr1_123_456.aln
    Brassica_rapa_FPsc_Brara.A01406.1_chr1_123_456.aln

    This function first removes the trailing chromosome/start/end fields and then
    removes the known focal species prefix. This allows geneIDs with underscores,
    such as XP_010545484.1, to be preserved.
    """

    base = os.path.basename(str(filename))
    base = re.sub(r"\.(aln|fa|fasta|fna)$", "", base)

    # Remove trailing _chromosome_start_end
    # chromosome can contain dots or other characters but not underscores
    m = re.match(r"(.+)_([^_]+)_([0-9]+)_([0-9]+)$", base)
    if not m:
        return None

    left = m.group(1)

    known_species = [
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

    for sp in sorted(known_species, key=len, reverse=True):
        if left.startswith(sp + "_"):
            gid = left[len(sp) + 1:]
            # Strip chromosome first-part that leaks when chromosome contains underscores
            # e.g. scaffold_1 -> "_scaffold" leaks; Lachesis_group3 -> "_Lachesis" leaks
            gid = re.sub(r"_[A-Za-z][A-Za-z0-9]*$", "", gid)
            return gid

    # Fallback: last token before chromosome/start/end
    return left.split("_")[-1]


# -----------------------------
# Load nodes summary
# -----------------------------
def load_nodes_summary(nodes_dir, species):
    prefix = species_to_nodes_prefix(species)

    patterns = [
        os.path.join(nodes_dir, f"{prefix}_nodes_summary_filtered.tsv"),
        os.path.join(nodes_dir, f"*{prefix}*nodes_summary_filtered.tsv"),
    ]

    matches = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            break

    if len(matches) == 0:
        print("WARNING: nodes file not found. Tried:")
        for p in patterns:
            print(f"  {p}")
        return None

    if len(matches) > 1:
        print("WARNING: multiple nodes files found, using first:")
        for m in matches:
            print(f"  {m}")

    fname = matches[0]
    print(f"Using nodes file: {fname}")

    df = pd.read_csv(fname, sep="\t")

    if "geneID" not in df.columns:
        print(f"WARNING: geneID column not found in nodes file: {fname}")
        return None

    df["geneID"] = df["geneID"].str.replace(r"\.p$", "", regex=True)
    return df.set_index("geneID")


# -----------------------------
# Load tblastn summary
# -----------------------------
def load_tblastn(tblastn_dir, species):
    prefix = species_to_tblastn_prefix(species)

    pattern = os.path.join(
        tblastn_dir,
        f"{prefix}_cdngs_tblastn_*_summary.tsv"
    )

    matches = sorted(glob.glob(pattern))

    if len(matches) == 0:
        print(f"WARNING: no tblastn file found with pattern: {pattern}")
        return None

    if len(matches) > 1:
        print("WARNING: multiple tblastn files found, using first:")
        for m in matches:
            print(f"  {m}")

    fname = matches[0]
    print(f"Using tblastn file: {fname}")

    df = pd.read_csv(fname, sep="\t")

    if "geneID" not in df.columns:
        print(f"WARNING: geneID column not found in tblastn file: {fname}")
        return None

    df["geneID"] = df["geneID"].str.replace(r"\.p$", "", regex=True)
    return df.set_index("geneID")


# -----------------------------
# Main
# -----------------------------
def main(coord_dir, nodes_dir, tblastn_dir, outdir):

    coord_dir = os.path.abspath(coord_dir)
    nodes_dir = os.path.abspath(nodes_dir)
    tblastn_dir = os.path.abspath(tblastn_dir)
    outdir = os.path.abspath(outdir)

    print(f"Coordinates dir: {coord_dir}")
    print(f"Nodes dir: {nodes_dir}")
    print(f"TBLASTN dir: {tblastn_dir}")
    print(f"Output dir: {outdir}")

    os.makedirs(outdir, exist_ok=True)

    for file in sorted(os.listdir(coord_dir)):

        if not file.endswith("_all_species_summary.tsv"):
            continue

        infile = os.path.join(coord_dir, file)
        print(f"\nProcessing: {infile}")

        # Detect focal species from coordinate input filename.
        # Example:
        # Anemorensis_DNG_45k_cds_all_species_summary.tsv -> Arabis_nemorensis
        short_prefix = file.replace("_DNG_45k_cds_all_species_summary.tsv", "")
        focal_species = canonicalize_species(short_prefix)

        print(f"Detected focal species: {focal_species}")

        df = pd.read_csv(infile, sep="\t")

        if "file_name" not in df.columns:
            print(f"WARNING: skipping file without file_name column: {infile}")
            continue

        if "species_name" not in df.columns:
            print(f"WARNING: skipping file without species_name column: {infile}")
            continue

        # Add geneID and focal_species columns
        df.insert(0, "geneID", df["file_name"].apply(extract_geneID))
        df.insert(1, "focal_species", focal_species)

        bad_geneids = df["geneID"].isna().sum()
        if bad_geneids:
            print(f"WARNING: {bad_geneids} rows had geneID=None")

        nodes_df = load_nodes_summary(nodes_dir, focal_species)
        tblastn_df = load_tblastn(tblastn_dir, focal_species)

        if nodes_df is None or tblastn_df is None:
            print(f"WARNING: missing nodes or tblastn file for {focal_species}")
            continue

        # Append selected node columns
        for col in ["ORF_Node", "Synteny_Node", "Enablers_Nodes", "Gene_type"]:
            if col not in nodes_df.columns:
                print(f"WARNING: {col} not found in nodes summary for {focal_species}")
                df[col] = "NA"
            else:
                df[col] = df["geneID"].map(nodes_df[col])

        # -----------------------------
        # Clean Enablers_Nodes column
        # -----------------------------
        df["Enablers_Nodes"] = (
            df["Enablers_Nodes"]
            .replace(["", "na", "Na", "NA", "NaN", "nan", "None", None], pd.NA)
            .fillna("NA")
        )

        # -----------------------------
        # Remove rows without ORF_Node
        # -----------------------------
        before_orf_filter = len(df)

        df = df[df["ORF_Node"].notna()]
        df = df[df["ORF_Node"] != "NA"]
        df = df[df["ORF_Node"] != ""]

        after_orf_filter = len(df)

        print(f"Removed {before_orf_filter - after_orf_filter} rows without ORF_Node")

        # Filter by tblastn presence/absence
        keep = []

        for _, row in df.iterrows():
            gene = row["geneID"]
            sp = canonicalize_species(row["species_name"])

            if gene not in tblastn_df.index:
                keep.append(False)
                continue

            if sp not in tblastn_df.columns:
                keep.append(False)
                continue

            val = tblastn_df.loc[gene, sp]

            try:
                keep.append(float(val) > 0)
            except Exception:
                keep.append(False)

        kept_df = df[keep].copy()

        out_prefix = species_to_file_prefix(focal_species)
        outfile = os.path.join(outdir, f"{out_prefix}_coords_kept.tsv")

        kept_df.to_csv(outfile, sep="\t", index=False)

        print(f"Rows after ORF_Node filter: {len(df)}")
        print(f"Rows kept after tblastn filter: {len(kept_df)}")
        print(f"Written: {outfile}")


if __name__ == "__main__":

    if len(sys.argv) != 5:
        print("Usage:")
        print("python nodes_by_coordinates.py <coord_dir> <nodes_dir> <tblastn_dir> <outdir>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
