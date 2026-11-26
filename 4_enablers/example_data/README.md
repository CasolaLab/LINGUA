# Stage 4 example data

`workdir/` is the stage 4 input of the seven tutorial genes (see [`docs/tutorial/05_stage4_enablers.md`](../../docs/tutorial/05_stage4_enablers.md)),
in the layout the paper's final run assembled from its stage 3 results, namely the step 9 alignments, the step 13 stop and
step 14 indel tables, the first five columns of step 17 and the candidate proteins, restricted to the seven loci. `scripts/stage3_to_stage4.py` writes the same layout from the stage 3 tutorial run, with its own folder and file names
and the stop and indel tables of the current stage 3 scripts (shorter, see the stage 3 tutorial page); stage 4 gives
identical final tables from either.
`expected/` holds the output of `ONLY=athaliana REMOVE_LIST= bash 4_enablers/run_pipeline.sh <workdir> 2 9`, whose
`stage4_final` rows equal those of the paper's final run for these genes.
