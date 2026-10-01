# LINGUA: LINeage-specific Gene Uniform Annotator

LINGUA identifies candidate de novo genes across a clade, extracts their syntenic regions
from a reference-free whole-genome alignment, and classifies each locus as a de novo gene
with an enabler mutation placed on the phylogeny (DNG) or a putative de novo gene (PDNG).
It was developed for, and used in,

> Owoyemi A. O., Sanders S. C., Chaison C. A., Marano M. A. & Casola C. Phylogenomic
> reconstruction across the Brassicaceae reveals the tempo, genomic determinants, and
> functional emergence of de novo genes. (manuscript; citation to be added on publication)

## Stages

LINGUA runs in four stages. Each folder is self-contained, with its own README, scripts and
example or reference data, and the output of each stage is the input of the next.

| Stage | Folder | What it does | Main output |
|---|---|---|---|
| 1 | [`1_protein_preprocessing/`](1_protein_preprocessing/) | Protein quality control and longest-isoform selection from per-species GFF and protein FASTA files | One standardized, non-redundant proteome per species |
| 2 | [`2_gene_stratigraphy/`](2_gene_stratigraphy/) | Places every gene on the species tree from OrthoFinder hierarchical orthogroups and extracts candidate lineage-specific genes | Candidate de novo gene (cDNG) sequences per species |
| 3 | [`3_alignment_extraction/`](3_alignment_extraction/) | Extracts each candidate locus from the Progressive Cactus HAL alignment with mafExtractor, removes alignment and extraction artifacts, stitches exons and tabulates start, stop and frameshift changes per species | Per-species, exon-stitched alignments with indel and stop-codon tables |
| 4 | [`4_enablers/`](4_enablers/) | Scores intact ORFs, synteny and disablers on every species and reconstructed ancestral node, places enablers and classifies DNGs and PDNGs | Per-species node summaries and family-wide DNG and PDNG catalogs |

Between stages 2 and 3, candidates were filtered by homology searches against broad protein
databases, an organelle screen, RPS-BLAST against CDD, InterProScan and a BLASTp screen
against transposable-element proteins, as described in the paper. The loci removed by the
domain and transposable-element screens are listed in `4_enablers/data/lsgs_to_remove.txt`.

## Reproducing the published analysis

The tag `paper-published` marks the stage 4 code exactly as used for the paper (commit
`f2e5a02`). `4_enablers/README.md` describes how that run was validated and the one later
change to synteny-node assignment.

## Data

- Progressive Cactus alignment of 21 Brassicaceae and two outgroups,
  https://doi.org/10.5281/zenodo.22137464
- Proteomes and OrthoFinder results, https://doi.org/10.5281/zenodo.22926478
- De novo gene catalog and gap-free alignments, https://doi.org/10.5281/zenodo.22926717

## Authors

Stages 1 and 2, Adekola O. Owoyemi. Stage 3, Sierra C. Sanders. Stage 4, the Casola Lab.
Stages 1 to 3 were developed in separate repositories
([protein-preprocessing-isoform-pipeline](https://github.com/Ludtson/protein-preprocessing-isoform-pipeline),
[gene-stratigraphy-pipeline](https://github.com/Ludtson/gene-stratigraphy-pipeline),
[progressive-cactus-gene-analysis-pipeline](https://github.com/sierras64/progressive-cactus-gene-analysis-pipeline))
and were merged here with their full commit history.

## License

MIT, see `LICENSE`. This applies to every stage.

## Contact

Casola Lab, Department of Ecology and Conservation Biology, Texas A&M University
(ccasola@tamu.edu).
