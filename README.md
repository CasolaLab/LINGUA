# LINGUA: LINeage-specific Gene Uniform Annotator

LINGUA identifies candidate de novo genes across a clade, extracts their syntenic regions
from a reference-free whole-genome alignment, and classifies each locus as a de novo gene
with an enabler mutation placed on the phylogeny (DNG) or a putative de novo gene (PDNG).
It was developed for, and used in,

> Owoyemi A. O., Sanders S. C., Chaison C. A., Marano M. A. & Casola C. Phylogenomic
> reconstruction across the Brassicaceae reveals the tempo, genomic determinants, and
> functional emergence of de novo genes. (manuscript; citation to be added on publication)

## Stages

LINGUA runs in five stages, numbered 0 to 4. Stage 0 is a preparatory step that standardizes
the input proteomes. Each folder is self-contained, with its own README, scripts and example or
reference data, and the output of each stage is the input of the next.

| Stage | Folder | What it does | Main output |
|---|---|---|---|
| 0 | [`0_protein_preprocessing/`](0_protein_preprocessing/) | Protein quality control and longest-isoform selection from per-species GFF and protein FASTA files | One standardized, non-redundant proteome per species |
| 1 | [`1_gene_stratigraphy/`](1_gene_stratigraphy/) | Places every gene on the species tree from OrthoFinder hierarchical orthogroups and extracts candidate lineage-specific genes | Candidate de novo gene (cDNG) sequences per species |
| 2 | [`2_homolog_exclusion/`](2_homolog_exclusion/) | Removes candidates with a known protein domain (InterProScan, RPS-BLAST against CDD), a homolog outside the clade (BLAST+, jackhmmer), an organelle-derived origin, or similarity to transposable-element proteins over at least half of their length (BLASTp against a TE protein library) | Candidates with no evidence of homology outside the clade |
| 3 | [`3_alignment_extraction/`](3_alignment_extraction/) | Extracts each candidate locus, in every species and in the ancestral genomes reconstructed by Progressive Cactus, from the HAL alignment with mafExtractor, removes alignment and extraction artifacts, stitches exons and tabulates start, stop and frameshift changes per species | Per-species, exon-stitched alignments with indel and stop-codon tables |
| 4 | [`4_enablers/`](4_enablers/) | Scores intact ORFs, synteny and disablers on every species and on the Progressive Cactus ancestral sequences of every internal node, places enablers and classifies DNGs and PDNGs | Per-species node summaries and family-wide DNG and PDNG catalogs |

In the published run, the loci removed by the domain and transposable-element screens are listed
in `4_enablers/data/lsgs_to_remove.txt`.

Each stage has its own software requirements, listed in its README. `environment.yml` at the top level is the
conda environment for stage 4.

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

Stages 0 to 2, Adekola O. Owoyemi. Stage 3, Sierra C. Sanders. Stage 4, the Casola Lab.
Stages 0 to 3 were developed in separate repositories
([protein-preprocessing-isoform-pipeline](https://github.com/Ludtson/protein-preprocessing-isoform-pipeline),
[gene-stratigraphy-pipeline](https://github.com/Ludtson/gene-stratigraphy-pipeline),
[homolog-exclusion-pipeline](https://github.com/Ludtson/homolog-exclusion-pipeline),
[progressive-cactus-gene-analysis-pipeline](https://github.com/sierras64/progressive-cactus-gene-analysis-pipeline))
and were merged here with their full commit history.

## License

MIT, see `LICENSE`. This applies to every stage.

## Contact

Casola Lab, Department of Ecology and Conservation Biology, Texas A&M University
(ccasola@tamu.edu).
