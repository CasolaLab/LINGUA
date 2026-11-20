# 07 Running LINGUA on your own clade

LINGUA was built on 21 Brassicaceae and two outgroups and has since been run on bryophytes and on T2T mammal genomes.
Most of the code is clade-independent, but a few inputs encode the tree and the species, and they must be rebuilt for
any other group. This page lists them in the order the stages need them.

## Before you start, choose the species

- **Ingroup and outgroups.** De novo genes are sought inside a densely sampled ingroup. Outgroups outside it are needed
  to polarize mutations, because an enabler can only be placed when relatives on both sides of the branch keep the
  ancestral state. More than one outgroup, at increasing distance, gives more of these comparisons.
- **Several genomes per lineage.** In the Brassicaceae, 78.5% of DNGs (1,864 of 2,375) rely on a single extant genome on
  at least one side of the tree, so losing one genome loses DNGs. Genomes spread over the whole tree help more than many
  genomes from one genus.
- **Annotation quality.** Stage 0 starts from GFF3 and protein files, and the summary stage uses the GFF3 again for the
  overlap step. Check completeness (for example with BUSCO) before stage 0.

## Stage by stage

| Stage | What to change |
|---|---|
| 0 | nothing; give one GFF3 and one protein FASTA per species, and use `--format-headers --add-species` so species names travel with the IDs |
| 1 | run OrthoFinder with your rooted species tree (`-s`), and choose the HOG file of the node that holds the ingroup (`--hog-filename`, see [02](02_stage1_stratigraphy.md)) |
| 2 | build your own databases with `make-db`, leaving out the analysis species and their close relatives (`--exclude-inputs`); choose a TE protein library for your clade; give organelle genomes of your clade, or leave the organelle screen out with a reason (`SKIP_SCREENS=organelle`) |
| Cactus | align all genomes, ingroup and outgroups, with Progressive Cactus and your species tree; the tree given to stage 3 and to `make_tree_config.py` must use the ancestor names found in the alignment |
| 3 | one run per reference species, with that species' BED and CDS and the same tree |
| 4 | write a tree configuration and a species table (below) |
| summary | a branch-length table from a time-calibrated tree (`branch`, `length` in million years), a GFF table, and your HOG table |

## The stage 4 tree configuration

Every stage 4 and summary-stage script reads the tree from one file, `tree_config.json`, written by
`4_enablers/scripts/make_tree_config.py` from a rooted Newick tree with named tips and named internal nodes.

```bash
python3 4_enablers/scripts/make_tree_config.py my_tree.nwk --ingroup-crown <crown node> --outdir my_config
```

- The internal nodes are relabeled N0, N1, ... from the root down, so a parent always has a lower number than its
  descendants. `node_name_map.tsv` in the output links each new label to the original name and its species.
- `--ingroup-crown` names the crown node of the ingroup. By default no ORF may be placed at the crown, its ancestors or
  outside it, and no enabler at the ancestors of the crown except its parent.
- `--allow-orf-node` lets named nodes carry an ORF anyway, when outgroups polarize them. The Brassicaceae default allows
  N2 and N3 (`--ingroup-crown N3 --allow-orf-node N2 --allow-orf-node N3`).
- `--conditional-enabler-node NODE=REQ` lets an otherwise excluded node count as an enabler node only for a disabler that
  the required node also carries (used for the root of the T2T mammal tree).

Then write a `species.tsv` for your species with the same four columns as `4_enablers/data/species.tsv` (short query
prefix, gap-free folder name, output prefix, full species name as in the alignment), and point stage 4 at the three new
files. Steps 4, 5, 7 and 8 also read the Newick tree itself, and its default is the Brassicaceae tree, so all three
settings are needed.

```bash
LINGUA_CONFIG=my_config/tree_config.json TREE=my_config/species_tree_N.nwk SPECIES_TSV=my_species.tsv \
  bash 4_enablers/run_pipeline.sh WORKDIR 2 9
```

## Lessons from the bryophyte and mammal runs

- **Gene IDs shared between species.** Some annotations reuse the same gene ID in different species. The summary stage
  keys every gene by species and ID, and it stops when one ID maps to two genes, rather than merging them silently.
- **ID suffixes.** Phytozome protein IDs often end in `.p` while the GFF and BED IDs do not. Stage 3 and the summary
  stage match IDs with and without the suffix, and stop when a join leaves IDs unmatched. Check that the IDs of stage 0,
  the BED and the GFF agree before a long run.
- **Several ingroups in one alignment.** When one alignment holds several densely sampled clades (primates, rodents and
  bovids in the T2T run), run stage 4 and the summary stage once per clade, each with its own tree configuration.
- **Organelle and TE screens.** Give the organelle genomes (mitochondrial only, in animals) and a TE protein library that
  fit your clade. Leaving a screen out is recorded in the stage 2 `STATUS`.

Next, [08 Troubleshooting](08_troubleshooting.md).
