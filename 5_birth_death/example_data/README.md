# Summary-stage example data

Inputs of the tutorial page [`docs/tutorial/06_stage5_birth_death.md`](../../docs/tutorial/06_stage5_birth_death.md),
which reads the stage 4 example output from `tutorial_runs/stage4/`.

- `hogs_subset.tsv` holds the header and the three rows of the paper's OrthoFinder HOG table (`N0.tsv` of 2025-10-16)
  whose *A. thaliana* column holds a tutorial gene.
- `Arabidopsis_thaliana_subset.gff3` holds the genes of the Araport11 annotation (Phytozome release 447, the one the
  paper's overlap step used) that lie within 30 kb of the seven loci, with all their features; `gff_table.tsv` points
  to it.
- `stage5_example.conf` is the run configuration; `expected/` holds its output.
