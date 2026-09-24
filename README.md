# LINGUA: LINeage-specific Gene Uniform Annotator

LINGUA detects enabler mutations in whole-genome alignments and classifies candidate
de novo genes as DNGs or PDNGs across the Brassicaceae, as used in

> Owoyemi A., Sanders S., Chaison C. A., Marano M. A. & Casola C. Phylogenomic
> reconstruction across the Brassicaceae reveals the tempo, genomic determinants, and
> functional emergence of de novo genes. (manuscript; citation to be added on publication)

For each candidate de novo gene (cDNG) the pipeline finds the oldest node of the species
tree at which the open reading frame (ORF) is intact, the oldest node at which the locus is
syntenic, and the nodes that carry the frameshifting indels and premature stop codons
(disablers) whose removal created the ORF. A locus is a DNG when an enabler can be placed on
the branch immediately ancestral to the ORF, and a PDNG otherwise.

## Where this repository sits

This repository holds the enabler and node-assignment stage of LINGUA, the last stage of
the analysis. It starts from outputs produced upstream, which are not part of this
repository.

1. Candidate identification. Proteomes of 23 species were clustered with OrthoFinder
   v3.0.1b1 (default settings, DIAMOND search), and Brassicaceae-restricted genes were
   filtered by homology searches against broad protein databases, an organelle screen,
   RPS-BLAST against CDD, InterProScan and a BLASTp screen against transposable-element
   proteins. The loci removed by the domain and TE screens are listed in
   `data/lsgs_to_remove.txt`, which step 7 applies.
2. Whole-genome alignment. A reference-free Progressive Cactus alignment of 21 Brassicaceae
   and two outgroups, deposited at Zenodo (https://doi.org/10.5281/zenodo.22137464).
3. LINGUA alignment extraction. MAF blocks for every candidate locus were extracted from the
   alignment, filtered for extraction artifacts and stitched across exons, and indels and stop
   codons were tabulated per species and per reconstructed ancestral node. This stage is
   maintained separately and writes one result archive per species; these archives are the
   input to step 1 here.

## Requirements

- Python 3.10 or later with pandas
- BLAST+ (`tblastn`, `makeblastdb`) on the PATH, or set `TBLASTN` and `MAKEBLASTDB`

A conda environment is provided.

```
conda env create -f environment.yml
conda activate lingua
```

## Input layout

```
WORKDIR/
  _all_species_data/                   LINGUA result archives, one .zip per species
  clsgs_45k_1line/                     cDNG proteins, one single-line FASTA per species
  coordinates_alignments_all_species/  MAF block coordinates per species (step 9 only)
```

Species names, file prefixes and the tree labels they map to are listed in
`data/species.tsv`. The species tree with internal nodes N0 to N21 is
`data/SpeciesTree_rooted_node_labels_oct_16.txt`.

## Running

```
bash run_pipeline.sh WORKDIR            # all steps
bash run_pipeline.sh WORKDIR 4 4        # only step 4
ONLY=athaliana bash run_pipeline.sh WORKDIR   # restrict step 4 to one species
```

Set `PYTHON`, `TBLASTN` and `MAKEBLASTDB` to use specific executables.

## Steps

| Step | Script | What it does | Main output |
|---|---|---|---|
| 1 | `unzip_and_move_LINGUA_results.py` | Unpacks the LINGUA archives and sorts alignments, indel tables and stop-codon tables into separate folders | `brassicaceae_23species_LINGUA_results/`, `indels_dir/`, `stops_dir/` |
| 2 | `remove_gaps_batch_recursive.py` | Removes alignment gaps from every per-species and per-node sequence of each exon-stitched locus | `brassicaceae_23species_LINGUA_nogaps/` |
| 3 | `rename_file_getIDs.py` | Renames the candidate protein files and tabulates gene IDs by species | `45k_renamed/`, `brassicaceae_geneIDs_species_names.tsv` |
| 4 | `blast_prts-vs-nogaps_AncSeq_v2.py` | tBLASTn of each candidate protein against the gap-free sequence of every species and ancestral node in its alignment, scoring ORF presence and the oldest node with an intact ORF | `cdngs_tblastn_results/*_cdngs_tblastn_<n>{,_summary,_AncSeq}.tsv` |
| 5 | `summary_synteny_tree_AncSeq.py` | Scores synteny conservation per species and node and reports the oldest syntenic node | `synteny_40perc/*_nogaps-synteny-40{,_summary,_AncSeq}.tsv` |
| 6 | `indel_oldest_shared_node_AncSeq_v2.py`, `stops_oldest_shared_node_AncSeq_v2.py` | Assigns indels and premature stop codons shared by descendant sequences to nodes | `indels_dir/all_indels_nodes_<species>.tsv`, `stops_dir/all_stops_nodes_<species>.tsv` |
| 7 | `node_summary_tree_v2.py`, then `filter_nodes_summary_by_geneids.py` | Combines ORF, synteny, indel and stop nodes, classifies each locus, renames each summary to the short species prefix that step 9 expects, then removes the loci in `data/lsgs_to_remove.txt` | `nodes_summary_by_species/<species>_nodes_summary_filtered.tsv` |
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

## Notes on reproducing the published analysis

- Disablers enter step 6 exactly as tabulated by LINGUA, with no filter on their position
  along the alignment.
- The default synteny threshold of `summary_synteny_tree_AncSeq.py` is 40%, the value used
  in the paper, and `run_pipeline.sh` also passes it explicitly.
- `data/lsgs_to_remove.txt` (18,431 gene IDs) is the list of loci removed by the upstream
  CDD, InterProScan and transposable-element screens. It is applied after classification,
  in step 7, as in the published run. Each locus is classified independently, so the
  timing of this removal does not change any per-locus result.
- The BLAST+ version used for the published run was not recorded. A rerun with BLAST+ 2.16.0
  reproduced it exactly (see Validation).

## Validation

The complete pipeline was rerun on *Arabidopsis thaliana* (1,078 candidates) from the LINGUA
archive. Every intermediate and final table was byte-identical to the published run, the
only exception being the per-HSP tBLASTn table, which holds the same rows in a different
order. The final classification, 281 DNGs and 627 PDNGs among 908 loci, matches
Supplementary Table 1 of the paper on every column for every locus.

## Data

- Progressive Cactus alignment, https://doi.org/10.5281/zenodo.22137464
- Proteomes and OrthoFinder results, https://doi.org/10.5281/zenodo.22926478
- De novo gene catalog and gap-free alignments, https://doi.org/10.5281/zenodo.22926717

## License

MIT, see `LICENSE`.

## Contact

Casola Lab, Department of Ecology and Conservation Biology, Texas A&M University
(ccasola@tamu.edu).
