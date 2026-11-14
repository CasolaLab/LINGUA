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
The tree settings are `data/tree_config.json` (written by `scripts/make_tree_config.py`), and
the method options file decides which classification rules apply; the variables
`LINGUA_CONFIG` and `LINGUA_OPTIONS` point to other files.

## Rule sets

| Options file | Rules |
|---|---|
| `data/method_options.json` (default) | the rules of the T2T mammal analysis (2026-11-13), adopted as the LINGUA default for new runs |
| `presets/method_options_brassicaceae_published.json` | the rules of the published Brassicaceae catalog (1,811 DNGs) |

```
bash 4_enablers/run_pipeline.sh WORKDIR                                   # default rules
LINGUA_OPTIONS=4_enablers/presets/method_options_brassicaceae_published.json \
    bash 4_enablers/run_pipeline.sh WORKDIR                               # published Brassicaceae rules
```

The options file of every run is copied to `WORKDIR/logs/method_options_used.json`. Each key is
described under Method options below.

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
REUSE_TBLASTN=1 bash 4_enablers/run_pipeline.sh WORKDIR 4 4   # rerun only step 4b on the existing tBLASTn tables
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
| 4b | `orf_length_node.py` (`orf_test` = `length`) | Replaces the tBLASTn ORF calls by a length test on every ancestral sequence (the step 4 table is kept as `*_AncSeq.tblastn.tsv`) | rewritten `*_AncSeq.tsv`, `*_AncSeq.orf_length_details.tsv` |
| 5b | `add_start_codon_events.py` (`start_codon_events`) | Adds start-codon gains (a non-ATG codon at the focal ATG) to the stop tables as events at position 1 | stop tables (originals kept as `*.pre_start`) |
| 6 | `indel_oldest_shared_node_AncSeq_v2.py`, `stops_oldest_shared_node_AncSeq_v2.py` | Assigns indels and premature stop codons shared by descendant sequences to nodes | `indels_dir/all_indels_nodes_<species>.tsv`, `stops_dir/all_stops_nodes_<species>.tsv` |
| 6b | `filter_enabler_events.py` (`drop_terminal_indels`, `homopolymer_needs_extant`) | Removes indel events that are not evidence of a mutation (gap runs reaching either end of the focal CDS; homopolymer indels without an extant carrier) | filtered `all_indels_nodes_*`, `indels_dir/removed_enabler_events_<species>.tsv` |
| 6c | `step6b_orf_bearing_nodes_filter.py --drop-whole-event` (`orf_bearing_nodes_filter`) | Drops every disabler event carried by a node whose ancestral sequence still bears the ORF (inputs kept as `*.pre_orf_bearing.tsv`) | filtered `all_indels_nodes_*`, `all_stops_nodes_*` |
| 7 | `node_summary_tree_v2.py`, then `filter_nodes_summary_by_geneids.py` | Combines ORF, synteny, indel and stop nodes and classifies each locus, renames each summary to the short species prefix that step 9 expects, then (step 7b) removes the loci in `REMOVE_LIST` | `stage4_rules/nodes_summary_by_species/` |
| 7c | `remove_conserved_orfs.py` (`conserved_orf_filter`, on unless set to false) | Removes candidates, DNG or PDNG, whose ORF is present in an ancestor outside the ingroup (`CONSERVED_NODES`, by default the `excluded_orf_nodes` of the tree config) | `<prefix>_conserved_orf_removed.tsv` |
| 7d | `annotate_evidence.py` (`annotate_evidence`) | Annotates the evidence behind every enabler node, including support by any valid disabler at that node and by extant side-branch species | `<prefix>_evidence.tsv` next to each summary |
| 7e | `remove_literature_etes.py` (`remove_literature_etes`, off) | Removes known exapted transposable-element genes by BLASTp against a protein FASTA given as `ete_proteins`; written for the T2T mammals, off by default | `<prefix>_ete_removed.tsv` |
| 7f | `apply_option_b_reclassification.py` (`option_b_extant_support`) | Option B: a DNG keeps its class only when each enabler node is corroborated by an extant genome carrying the same event | `stage4_optionB/nodes_summary_by_species/` |
| 7g | `apply_consecutive_enablers.py` (`consecutive_enablers`) | A DNG's two enabler nodes must be the parent and the grandparent of its ORF node | `stage4_final/nodes_summary_by_species/` (the final classification) |
| 8 | `summarize_nodes_by_species_split.py` | Family-wide counts of DNGs by node and of DNGs and PDNGs by species, and one list of all classified loci | `brassicaceae_DNGs_by_age.tsv`, `brassicaceae_species_totals.tsv`, `brassicaceae_all_gene_types.tsv` |
| 9 | `nodes_by_coordinates.py` | Places DNGs and PDNGs at internal nodes by intersecting alignment coordinates across species | `parsed_coordinates/` |

Steps 4b, 5b, 6b, 6c, 7c, 7d, 7e, 7f and 7g run only when their option is on. In the published
catalog the steps now numbered 6c, 7d, 7f and 7g were numbered 6b, 7c, 7d and 7e. `ONLY` restricts
steps 4, 4b, 7c, 7d and 7e.

## Method options

| Key | Step | Effect | Default | Published |
|---|---|---|---|---|
| `orf_test` = `length`, `orf_min_len`, `orf_min_cov` | 4b | ancestral ORF presence by length instead of tBLASTn, namely the longest ORF in any forward frame starting at ATG or at the focal gene's own first codon, at least `orf_min_len` % of the focal CDS length and spanning at least `orf_min_cov` % of it, with no identity cutoff | length, 80, 80 | tBLASTn |
| `orf_count_stop` | 4b | an ancestral ORF ended by a stop codon includes it, as the focal CDS length does | on | off |
| `orf_min_in_frame` | 4b | an ancestor has the ORF only when, besides the length test, at least this % of the focal CDS is an unbroken focal-frame stretch in it (no stop in focal codon columns, no frameshifting gap run) | 50 | off |
| `orf_node_contiguous` | 4, 4b | the ORF node is the oldest ancestor of the unbroken run of ORF-bearing focal-lineage ancestors above the focal species | on | off |
| `start_codon_events` | 5b | a start-codon gain is an enabler, entered as a stop-table event at position 1 so that every stop rule applies to it | on | off |
| `indel_events_by_length` | 6 | an indel event is a position and a length, so different indels starting at one base are not merged | on | off |
| `contiguous_event_nodes` | 6 | an event is inherited only along contiguous lineage ancestors; each run of at least two carriers is a separate event | on | off |
| `event_gap_living_inferred` | 6 | a single lineage ancestor lacking an event between two carriers counts as a carrier when an extant species on its side branch carries the same event | on | off |
| `drop_terminal_indels` | 6b | indel events whose gap run reaches the first or last focal CDS base are removed (unaligned ends, not deletions) | on | off |
| `homopolymer_needs_extant`, `homopolymer_support_any_outgroup` | 6b | homopolymer indel events count only when an extant species carries the same indel, on a carrier's side branch or, with the second key, anywhere outside the clade that lacks the event | on, on | off |
| `orf_bearing_nodes_filter` | 6c | a disabler event touching any ORF-bearing ancestor is dropped whole | on | on |
| `record_enabler_events` | 7 | records the CDS positions of the indel and stop events used as enablers | on | on |
| `skip_events_in_younger_nodes` | 7 | an event is not used when a lineage ancestor younger than the ORF node carries it | on | off |
| `orf_clade_consistency` | 7 | the ORF node moves younger while the chosen disabler is present inside the ORF clade | on | on |
| `conserved_orf_filter` | 7c | removes candidates with an ORF in an ancestor outside the ingroup | on | off |
| `annotate_evidence`, `evidence_any_valid_disabler` | 7d | writes the evidence tables; extant support counts any valid disabler at an enabler node | on, on | on, on |
| `remove_literature_etes`, `ete_proteins` | 7e | removes known exapted TE genes (mammals; needs an ETE protein FASTA) | off | off |
| `option_b_extant_support` | 7f | a DNG stays DNG only when every enabler node has an extant carrier; otherwise PDNG | on | on |
| `option_b_outer_side_branch` | 7d, 7f | when an enabler node's own side branch has no aligned sequence near the event, an extant carrier on an older side branch may stand in, each node through a different branch | on | off |
| `consecutive_enablers` | 7g | a DNG's enabler nodes must be the parent and grandparent of the ORF node; otherwise PDNG | on | on |

A key that is absent is off, except `conserved_orf_filter`. With the default rules the T2T
mammal analysis gave 3,034 DNGs and 2,059 PDNGs in 11 species, against 2,194 and 2,401 without the
last four rules adopted (`orf_min_in_frame`, `start_codon_events`, `event_gap_living_inferred`,
`option_b_outer_side_branch`; LINGUA_t2t, `stage4/README.md`). On the Brassicaceae (same stage 3
inputs, tBLASTn and synteny as the published run; 2026-11-13) the default rules give 2,064 DNGs and
20,145 PDNGs in the 20 ingroup species, against 1,811 and 20,565 with the published rules. Of the
published DNGs, 1,557 stay DNGs, 167 become PDNGs, 49 are removed by step 7c and 38 drop out at step
7 because the length test places their ORF at N0-N3; 493 published PDNGs become DNGs, and 14 loci
that the published rules dropped at step 7 are now DNGs.

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
is not a disabler of that ORF and is dropped (step 6c). Second, an event used as an enabler
must be absent from the ORF clade, namely from every node inside the clade of the ORF node and
from every extant species in it other than the focal one; otherwise the ORF node moves one
node younger and the events are chosen again (step 7, `orf_clade_consistency`). Third, both
enabler nodes must be corroborated by at least one extant species on their side branch that
carries the same event, at the same position and, for indels, with the same length (option B,
step 7f). Fourth, the two enabler nodes must be the parent and the grandparent of the ORF node
(step 7g). Indels inside exons that stage 3 marks as missing data are not disablers (stage 3
step 14).

## Notes on reproducing the published analysis

`presets/method_options_brassicaceae_published.json` reproduces the NEE catalog of 1,811 DNGs
(the defaults before 2026-11-13; the code of that time is tag `nee-catalog-stages3-4`). The
earlier classification of the first submission (2,611 DNGs, without steps 6c and 7d-7g, the
stage 3 step 8c and the missing-data mask) is the stage 4 code at tag `paper-published`
(commit `f2e5a02`).

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

T2T scripts merged in (2026-11-13). Steps 6-9 were rerun on Dayhoff with
`presets/method_options_brassicaceae_published.json` on the same inputs as the 2026-11-09 run
below. All 233 output files hold the same rows as the 2026-11-09 run of tag
`nee-catalog-stages3-4`, giving 1,811 DNGs and 20,565 PDNGs in the 20 ingroup species. The new
steps were first used on the T2T mammals (LINGUA_t2t); steps 4-9 under the default rules were then run
on the Brassicaceae (all 23 species; step 7e off), with every locus of the published run accounted for
(kept, reclassified, removed by step 7c, or dropped at step 7 with its new ORF node at N0-N3). Loci
without a stage 3 alignment (8 of the 45,491 loci of the step 4 tables, all without an ORF in any ancestor) keep their step 4 row in
step 4b.

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
