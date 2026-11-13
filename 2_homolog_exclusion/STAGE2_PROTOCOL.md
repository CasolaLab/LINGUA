# Stage 2 protocol of the Brassicaceae de novo gene study (re-screen of November 2026)

This file records the stage-2 rules used for the 20 Brassicaceae genomes plus *Aethionema arabicum*, *Tarenaya hassleriana* and
*Theobroma cacao*, how to run them with this module, and how the study's own numbers were produced. The settings are in
`presets/brassicaceae_2026.conf` (homology and domain searches, HOG rule) and `presets/brassicaceae_2026_te.conf` (TE screen).

## Input

192,193 candidate lineage-specific genes (the longest protein isoform of each) from 23 genomes, namely genes in
orthogroups without members outside the Brassicaceae plus genes OrthoFinder left unclustered (stage 1). Every screen searches
every candidate, so no gene is left unsearched because another screen removed it.

## Databases

| Database | Proteins | Taxa |
|---|---|---|
| prokaryotes | 137,634 | 47 |
| protists | 260,343 | 15 |
| metazoans | 284,185 | 8 |
| fungi | 369,049 | 53 |
| green plants (incl. *Tarenaya hassleriana*, whose hits are ignored) | 1,962,820 | 73 |

The same five protein FASTA files serve BLASTp and jackhmmer. Organelle screen: 19 mitochondrial genomes of Brassicaceae and 22
chloroplast genomes (21 Brassicaceae and *T. cacao*). TE screen: 23,230 TE proteins (1,459 from TREP and 21,771 retrotransposon
and transposase entries of UniProtKB restricted to Brassicaceae). PANTHER family membership: PANTHER 19.0 sequence
classification files (143 genomes), built into a table with `panther-families`.

## What counts as evidence

| Screen | A gene is removed when | Preset keys |
|---|---|---|
| BLASTp | hits with E ≤ 1e-5 from at least two genera (search E 1e-3, `-seg yes`, ten targets per database) | `filter_evalue`, `min_genera`, `max_target_seqs` |
| jackhmmer | one iteration (`-N 1 -E 1e-4 --incE 1e-5`); full-sequence hits with E ≤ 1e-5 from at least two genera, with no composition-bias filter (`max_bias_frac` blank since 0.3.3) | `iterations`, `filter_evalue`, `min_genera` |
| RPS-BLAST vs CDD | a hit with E ≤ 0.01 to any domain except four Brassicaceae-restricted DUFs (pfam02713 DUF220, pfam05212 DUF707, pfam06683 DUF1184, pfam06880 DUF1262) | `ignore_domains`, `min_domain_cov=0` |
| InterProScan | a hit with E ≤ 1e-3 (PROSITE profiles and HAMAP at their own thresholds), ignoring Coils, MobiDB-lite, PROSITE patterns, AntiFam, the signatures PS51257, PR01217 and PR01228, and the four DUFs | `max_evalue`, `require_evalue`, `exclude_analyses`, `ignore_domains` |
| PANTHER | a family hit (no E-value cutoff) whose family includes a species other than *A. thaliana*, *Brassica napus*, *B. rapa* and *T. cacao* | `max_evalue_exempt`, `panther_family_species`, `panther_exclude_species` |
| TE | one BLASTp alignment (word size 2, PAM30, E ≤ 1e-3) to a TE protein covering at least 50% of the candidate | `brassicaceae_2026_te.conf` |
| Organelle | any tBLASTn hit (E ≤ 1e-3, ungapped) to an organelle genome | `organelle_tblastn.conf` |
| HOG rule | after the above, the candidate's BLASTp alignment to a removed member of its HOG overlaps that member's hit region by at least 20 aa (members removed by BLASTp, jackhmmer, CDD, InterProScan or PANTHER) | `[hog-rule]` |
| Paralog rule | after the HOG rule, one BLASTp alignment (SEG off, E ≤ 1e-5, up to 1,000 targets) to a gene of the candidate's own proteome that is not a candidate has at least 50% identity over at least 50% of the candidate's length | `[paralog-rule]` |

Hits to the OrthoFinder outgroup genomes do not count (only *T. hassleriana* is in the databases).

## Commands

```bash
P=presets/brassicaceae_2026.conf; S=scripts/homolog_exclusion_pipeline.sh
bash $S panther-families PTHR19.0_* > panther_family_species.tsv
bash $S homology-blast     -i cand/ -o out/blast    --label blast    --config $P --db db/prokaryotes --db db/protists --db db/animals --db db/fungi --db db/plants
bash $S homology-jackhmmer -i cand/ -o out/jack     --label jackhmmer --config $P --db db/all_prt.fa
bash $S domain-cdd         -i cand/ -o out/cdd      --label cdd      --config $P --db Cdd
bash $S domain-interpro    -i cand/ -o out/ipr      --label interpro --config $P
bash $S homology-blast     -i cand/ -o out/te       --label te       --config presets/brassicaceae_2026_te.conf --db db/te_proteins
bash $S homology-blast     -i cand/ -o out/org      --label organelle --config presets/organelle_tblastn.conf --db db/organelle
bash $S combine --runs out/blast out/jack out/cdd out/ipr out/te out/org -o out/combine
bash $S hog-rule -i cand/ -o out/hog --combine out/combine --hogs N0.tsv --config $P
bash $S paralog-rule -i cand/ -o out/paralog --passed out/hog --proteomes proteomes/ --config $P
```

From 0.3.4, `scripts/run_stage2.sh` runs the same sequence from one config (`CONFIG=presets/brassicaceae_2026.conf`,
`TE_CONFIG=presets/brassicaceae_2026_te.conf`, `ORGANELLE_CONFIG=presets/organelle_tblastn.conf`); it was tested on a small
hand-built example, not rerun on the study's 192,193 candidates. The `panther-families` table is built beforehand and given as
`PANTHER_FAMILY_SPECIES`.

The study kept ten BLASTp targets per database and counted genera across all five databases together. From 0.3.2 the module
does both: with several `--db`, each database is searched as its own run (`blast-prokaryotes`, `blast-plants` and so on, ten
targets each), and one more run labelled `blast` pools their hits and applies `min_genera` across all five databases
(`min_genera_scope=all_databases`, the default). `combine` reads the pooled `blast` run, which the `blast*` pattern of the
preset matches, and skips the per-database runs. Before 0.3.2 genera were counted within each database. In the study every
database hit had a species label, so the module's rule that hits of unknown species count together as one genus never applied.

## Counts in the study

| Step | Genes |
|---|---|
| candidates | 192,193 |
| flagged by BLASTp / jackhmmer / CDD / InterProScan / PANTHER / TE / organelle | 114,241 / 109,688 / 115,966 / 113,982 / 113,816 / 46,318 / 8,118 |
| removed by their own hits (any screen) | 148,575 (15,817 by one screen only) |
| candidates with no counted hit of their own | 43,618 |
| in a HOG with a removed member | 16,182 |
| removed by the HOG rule | 10,422 |
| passing the HOG rule | 33,196 |
| removed by the paralog rule | 986 |
| lineage-specific genes | 32,210 |

Until 0.3.2 jackhmmer hits also needed composition bias below half the bit score; with that rule the same data gave 108,116
jackhmmer-flagged genes, 44,306 candidates with no counted hit, 34,020 passing the HOG rule and 33,014 lineage-specific genes.

## How the study's numbers were produced

The study did not run this module end to end. The re-screen combined new searches of every candidate (BLAST+ 2.17.0, HMMER 3.4,
RPS-BLAST against CDD 3.12, InterProScan 5.77-108.0 for PANTHER, PANTHER 19.0 classification files) with earlier searches of the
same proteins (BLASTp with BLAST+ 2.2.31, earlier jackhmmer runs, an InterProScan run of March 2026 whose version was not recorded,
and an RPS-BLAST run of March 2026), and project scripts merged them under the rules above. The rules were then added to this
module. What has been checked:

- with the 0.3.2 rules (jackhmmer bias rule on), the `hog-rule` decision, applied to the study's within-HOG BLASTp output and hit
  regions, reproduces the study's call for all 16,279 genes in HOGs with a removed member (10,286 removed; best overlap identical for every gene), and
  `hog-rule` running BLAST+ itself reproduces the verdicts and overlaps of four of them;
- `paralog-rule`, run with BLAST+ 2.17.0 on the 34,020 genes passing the HOG rule under the 0.3.2 rules, removes 1,006 genes with
  1,000 targets; with 50 targets it removed 954, the same 954 that the project's analysis of a separate BLASTp run removed. Its
  verdict depends only on a gene and its own proteome, so the 0.3.3 count reuses those verdicts for the 33,196 genes (all among the
  34,020), which removes 986;
- each new setting has a unit test with hand-worked answers (`tests/test_stage2_rules.py`).

Not checked: a full rerun of the 192,193 candidates through the module. Known differences from the project scripts are small
and listed here so they can be tested. The TE rule here measures an alignment's query span over the protein length, whereas the
study divided the alignment length by the protein length minus one; InterProScan matches must have status `T` here; and `min_genera`
is applied per database run unless the databases are merged (see Commands).
