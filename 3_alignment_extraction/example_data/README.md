# Stage 3 example data

Seven *Arabidopsis thaliana* candidates of the paper's final run (LINGUA v1.2.0, 32,210 lineage-specific genes), used by
the tutorial page [`docs/tutorial/04_stage3_alignment_extraction.md`](../../docs/tutorial/04_stage3_alignment_extraction.md).

- `Arabidopsis_thaliana_tutorial.bed`, `Arabidopsis_thaliana_tutorial_cds.fna`, `species_tree.nwk` are the BED rows,
  CDS records and tree of the *A. thaliana* stage 3 run of 2026-11-13 that the final stage 4 read, restricted to the
  seven transcripts (12 exons).
- `step1/` holds that run's step 1 MAF extractions of the 12 exons (copied from the run, not regenerated).
- `expected/` holds that run's own step 8c, 9 and 12 files of the seven loci and their `step17_all_species.tsv` rows,
  and the step 13 and 14 files written by the current scripts (commit 4fd67bc, 2026-11-15), which move the stage 4
  step 5a and 5c filters into stage 3. Stage 4 run on these tables or on the paper's gives identical final tables for
  the seven loci.
