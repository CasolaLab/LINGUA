# Stage 4. Enabler detection and node assignment

This folder holds the last stage of LINGUA (see the top-level `README.md`). It detects
enabler mutations in whole-genome alignments and classifies candidate de novo genes as DNGs
or PDNGs across the Brassicaceae.

For each candidate de novo gene (cDNG) the pipeline finds the oldest node of the species
tree at which the open reading frame (ORF) is intact, the oldest node at which the locus is
syntenic, and the nodes that carry the frameshifting indels and premature stop codons
(disablers) whose removal created the ORF. A locus is a DNG when an enabler can be placed on
the branch immediately ancestral to the ORF, and a PDNG otherwise.

## Inputs from the earlier stages

This stage starts from the stage 3 results, packed as one archive per species
(`../3_alignment_extraction/`), and the candidate proteins from stage 1
(`../1_gene_stratigraphy/`) screened in stage 2 (`../2_homolog_exclusion/`). Step 7b removes
the loci in `REMOVE_LIST` (default `data/lsgs_to_remove.txt`, the loci removed by the original
CDD, InterProScan and transposable-element screens of the published run). With stage 2 of
LINGUA v1.1 or later these loci never reach stage 3, so set `REMOVE_LIST=` (empty) to skip it.
The tree settings and method options that steps 6b-7e read are `data/tree_config.json`
(written by `scripts/make_tree_config.py`) and `data/method_options.json`; the variables
`LINGUA_CONFIG` and `LINGUA_OPTIONS` point to other files.

## Requirements

- Python 3.10 or later with pandas
- BLAST+ (`tblastn`, `makeblastdb`) on the PATH, or set `TBLASTN` and `MAKEBLASTDB`

A conda environment is provided.

```
conda env create -f ../environment.yml
conda activate lingua
```

## Input layout

```
WORKDIR/
  _all_species_data/                   Stage 3 result archives, one .zip per species
  clsgs_45k_1line/                     cDNG proteins, one single-line FASTA per species
  coordinates_alignments_all_species/  MAF block coordinates per species (step 9 only)
```

Species names, file prefixes and the tree labels they map to are listed in
`data/species.tsv`. The species tree with internal nodes N0 to N21 is
`data/SpeciesTree_rooted_node_labels_oct_16.txt`.

## Running

```
bash 4_enablers/run_pipeline.sh WORKDIR            # all steps (from the repository root)
bash 4_enablers/run_pipeline.sh WORKDIR 4 4        # only step 4
ONLY=athaliana bash 4_enablers/run_pipeline.sh WORKDIR   # restrict step 4 to one species
```

Set `PYTHON`, `TBLASTN` and `MAKEBLASTDB` to use specific executables.

## Steps

| Step | Script | What it does | Main output |
|---|---|---|---|
| 1 | `unzip_and_move_LINGUA_results.py` | Unpacks the Stage 3 archives and sorts alignments, indel tables and stop-codon tables into separate folders | `brassicaceae_23species_LINGUA_results/`, `indels_dir/`, `stops_dir/` |
| 2 | `remove_gaps_batch_recursive.py` | Removes alignment gaps from every per-species and per-node sequence of each exon-stitched locus | `brassicaceae_23species_LINGUA_nogaps/` |
| 3 | `rename_file_getIDs.py` | Renames the candidate protein files and tabulates gene IDs by species | `45k_renamed/`, `brassicaceae_geneIDs_species_names.tsv` |
| 4 | `blast_prts-vs-nogaps_AncSeq_v2.py` | tBLASTn of each candidate protein against the gap-free sequence of every species and ancestral node in its alignment, scoring ORF presence and the oldest node with an intact ORF | `cdngs_tblastn_results/*_cdngs_tblastn_<n>{,_summary,_AncSeq}.tsv` |
| 5 | `summary_synteny_tree_AncSeq.py` | Scores synteny conservation per species and node and reports the oldest syntenic node | `synteny_40perc/*_nogaps-synteny-40{,_summary,_AncSeq}.tsv` |
| 6 | `indel_oldest_shared_node_AncSeq_v2.py`, `stops_oldest_shared_node_AncSeq_v2.py` | Assigns indels and premature stop codons shared by descendant sequences to nodes | `indels_dir/all_indels_nodes_<species>.tsv`, `stops_dir/all_stops_nodes_<species>.tsv` |
| 6b | `step6b_orf_bearing_nodes_filter.py --drop-whole-event` | Drops every disabler event carried by a node whose ancestral sequence still bears the ORF (the pre-filter tables are kept as `*_unfiltered.tsv`) | filtered `all_indels_nodes_*`, `all_stops_nodes_*` |
| 7 | `node_summary_tree_v2.py`, then `filter_nodes_summary_by_geneids.py` | Combines ORF, synteny, indel and stop nodes and classifies each locus under the options in `data/method_options.json` (ORF-clade consistency, recorded enabler events), renames each summary to the short species prefix that step 9 expects, then (step 7b) removes the loci in `REMOVE_LIST` | `stage4_rules/nodes_summary_by_species/` |
| 7c | `annotate_evidence.py` | Annotates the evidence behind every enabler node, including support by any valid disabler at that node and support by extant side-branch species | `<prefix>_evidence.tsv` next to each summary |
| 7d | `apply_option_b_reclassification.py` | Option B: a DNG keeps its class only when both enabler nodes are corroborated by an extant genome carrying the same event | `stage4_optionB/nodes_summary_by_species/` |
| 7e | `apply_consecutive_enablers.py` | A DNG's two enabler nodes must be the parent and the grandparent of its ORF node | `stage4_final/nodes_summary_by_species/` (the final classification) |
| 8 | `summarize_nodes_by_species_split.py` | Family-wide counts of DNGs by node and of DNGs and PDNGs by species, and one list of all classified loci | `brassicaceae_DNGs_by_age.tsv`, `brassicaceae_species_totals.tsv`, `brassicaceae_all_gene_types.tsv` |
| 9 | `nodes_by_coordinates.py` | Places DNGs and PDNGs at internal nodes by intersecting alignment coordinates across species | `parsed_coordinates/` |

## Parameters and classification rules

| Setting | Value |
|---|---|
| tBLASTn | e-value 1, word size 2, BLOSUM80, no SEG filtering, no composition-based statistics, at most 10 HSPs and 50 target sequences |
| Intact ORF | the matched sequence, from its first methionine to the first stop codon or the end of the match, spans at least 80% of the query length (`--frac_ORFss 80`), at an identity of at least 40% (`--pident 40`); a stop codon is not required |
| Synteny | at least 40% of the focal coding sequence aligned (`--synteny 40`) |
| Shared disablers | an indel or stop is assigned only when present in at least two nodes older than the ORF node and not older than the synteny node |

A locus is a **DNG** when its ORF node is younger than its synteny node and the two nearest
disabler-carrying nodes are both older than the ORF node and within one node of it, so that
an enabler can be placed on the branch immediately ancestral to the ORF
(`ORF-Enablers_Node_dist` of 1). A locus meeting the first condition whose nearest disabler
nodes lie more than one node away is a **PDNG**. Loci whose ORF node is N0, N1, N2 or N3 are
excluded from the summaries in step 8.

The published catalog (NEE revision, 1,811 DNGs) adds four rules to this definition.
First, a disabler event carried by any node whose reconstructed sequence still bears the ORF
is not a disabler of that ORF and is dropped (step 6b). Second, an event used as an enabler
must be absent from the ORF clade, namely from every node inside the clade of the ORF node and
from every extant species in it other than the focal one; otherwise the ORF node moves one
node younger and the events are chosen again (step 7, `orf_clade_consistency`). Third, both
enabler nodes must be corroborated by at least one extant species on their side branch that
carries the same event, at the same position and, for indels, with the same length (option B,
step 7d). Fourth, the two enabler nodes must be the parent and the grandparent of the ORF node
(step 7e). Indels inside exons that stage 3 marks as missing data are not disablers (stage 3
step 14).

## Notes on reproducing the published analysis

The defaults reproduce the NEE catalog of 1,811 DNGs. The earlier classification of the
first submission (2,611 DNGs, without steps 6b and 7c-7e, the stage 3 step 8c and the
missing-data mask) is the stage 4 code at tag `paper-published` (commit `f2e5a02`).

- Disablers enter step 6 exactly as tabulated in Stage 3, with no filter on their position
  along the alignment.
- The default synteny threshold of `summary_synteny_tree_AncSeq.py` is 40%, the value used
  in the paper, and `run_pipeline.sh` also passes it explicitly.
- `data/lsgs_to_remove.txt` (18,431 gene IDs) is the list of loci removed by the upstream
  CDD, InterProScan and transposable-element screens. It is applied after classification,
  in step 7, as in the published run. Each locus is classified independently, so the
  timing of this removal does not change any per-locus result.
- `summary_synteny_tree_AncSeq.py` now restricts the ancestral-sequence synteny node to
  ancestors of the focal species. The published run, commit `f2e5a02`, took the oldest
  syntenic ancestor anywhere in the alignment, which placed 41 of 45,017 loci on a node
  off their own lineage (none in *A. thaliana*). To reproduce the published tables
  exactly, check out that commit.
- The BLAST+ version used for the published run was not recorded. A rerun with BLAST+ 2.16.0
  reproduced it exactly (see Validation).

## Validation

Final catalog (2026-11-09). Steps 6-9 were rerun with this code on the inputs of the
published run (tBLASTn, synteny and alignments of all 23 species; indel tables from stage 3
with step 8c and the missing-data mask, the *A. thaliana* table regenerated by
`../3_alignment_extraction/run_stage3.sh`). All 233 output files (classification summaries,
evidence tables, option B and consecutive-enabler results, step 8 tables and step 9
coordinates) hold the same rows as the published run, giving 1,811 DNGs and 20,565 PDNGs in the
20 ingroup species. The two node-list columns of the evidence tables are now written in
lineage order; the published run wrote them in an arbitrary order, which differed between
runs, with the same nodes in 258 cells. Two reruns with different Python hash seeds were
byte-identical.

Earlier classification (first submission):
The complete pipeline was rerun on *Arabidopsis thaliana* (1,078 candidates) from the Stage 3
archive. Every intermediate and final table was byte-identical to the published run, the
only exception being the per-HSP tBLASTn table, which holds the same rows in a different
order. The final classification, 281 DNGs and 627 PDNGs among 908 loci, matches
Supplementary Table 1 of the paper on every column for every locus.

## License

MIT, see `../LICENSE`.
