# Stage 5. Births, deaths and rates

Stage 4 classifies each candidate locus of each species as a DNG or a PDNG and gives the node at which its ORF
arose. The same ancestral locus is often classified separately in several descendant species, and duplication
makes one origin look like several. This stage turns the per-species classifications into origins and losses on
the branches of the species tree, and divides them by branch length to give rates.

It reads the tree from the stage 4 configuration (`4_enablers/data/tree_config.json` or the one you made with
`4_enablers/scripts/make_tree_config.py`), so it runs on any tree that stage 4 ran on. When stage 4 was run
separately for several clades of one tree (each with its own crown and ORF window), run stage 5 once per clade with
that clade's `tree_config.json`; a single run with one clade's configuration would leave out the origins of the
others.

## Steps

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `step1_coordinates.py` | Joins the stage 3 step 17 coordinates of each locus in every genome with the stage 4 classification and the stage 4 tBLASTn table, keeping a genome only when tBLASTn finds the ORF there (stage 4 step 9 without its species lists) | `coords.tsv` |
| 1a | `step1a_orf_presence.py` (`PRESENCE=orf`) | Replaces the tBLASTn presence by presence of the ORF in every extant genome under the stage 4 ORF rule (the thresholds of the stage 4 method options: length, coverage, in-frame share), and lists the genomes with no usable call (row missing, or bases at under 50% of the focal CDS) | `presence/*_cdngs_tblastn_*_summary.tsv`, `presence/unknown_presence.tsv`, `presence/presence_comparison.tsv` |
| 2 | `step2_loci.py` | Joins the genes of different species that occupy the same aligned position into one locus (union-find over coordinate overlap within the clade of the gene's ORF node, on either strand); a locus takes the oldest ORF node of its members, and is a DNG when any member is | `loci.tsv`, `genes.tsv` |
| 3 | `step3_overlap.py` | Classifies each locus against the protein-coding genes of the species' GFF3 (same-strand CDS overlap, same-strand overlap without CDS overlap, antisense overlap, or the nearest neighbor); only loci whose coding exons overlap coding exons of a same-strand gene are left out in step 5. Overlaps between two candidate loci do not count. Required: the run stops when a species has no GFF, unless it is named in `NO_GFF_SPECIES` (`--no-gff`), and step 5 needs `--overlap` (`--no-overlap` is for tests and comparisons only) | `overlap.tsv` |
| 4 | `step4_unclustered_families.py` + BLASTp | For loci in no OrthoFinder orthogroup, BLASTp all against all within each species (E ≤ 10⁻⁵); connected loci form a family and only the longest protein is kept, whatever its type or age (as in the paper) | `unclustered/unclustered_families.tsv` |
| 5 | `step5_rates.py` | Collapses loci sharing a HOG into one DNG origin and one PDNG origin per lineage, drops PDNGs whose every species also carries a DNG of that HOG, places each origin on the branch of its oldest locus, counts deaths by Dollo parsimony on the presence of step 1 or 1a, and divides by branch length | `rates_by_branch.tsv`, `origins.tsv`, `gene_status.tsv` |

## Rules

- **One origin per orthogroup.** Loci sharing a hierarchical orthogroup (HOG) of the OrthoFinder table are taken as
  one origin followed by duplication. A locus takes the HOG most frequent among its genes. Within a HOG, DNG loci
  form one origin and PDNG loci another.
- **One origin per lineage** (`--hog-span split`, the default). OrthoFinder groups proteins by similarity, not by
  position, so a HOG can join loci that no single origin explains, for instance a locus born at N5 and a
  species-specific locus in a species outside clade N5. Loci are taken oldest first; a locus joins the origin whose
  clade contains its branch, otherwise it starts a new origin (labeled `HOG#2`, `HOG#3`...) with its own birth and
  deaths. `--hog-span merge` counts one origin per HOG, as in the paper.
- **A PDNG paralog of a DNG is a duplicate.** A PDNG locus is dropped when every species in which it occurs also has
  a DNG of the same HOG (the DNG is the parent). Otherwise it stays in the PDNG origin of the HOG.
- **Birth branch.** The branch of the oldest locus of the origin: the internal node with the fewest ancestors, or the
  species for a species-specific locus. Loci whose ORF node is excluded in `tree_config.json` are not counted.
- **Deaths.** For an origin at internal node N, a species of clade(N) carries it when tBLASTn finds the ORF of any of
  its genes in that genome, an internal node when any descendant does, and N always. A death is placed on every
  branch whose parent carries the origin and which does not. Origins on a species branch have no deaths.
- **Presence.** With `PRESENCE=tblastn` (default; the published analysis) a genome carries a gene when stage 4
  tBLASTn finds its ORF there. With `PRESENCE=orf` (step 1a, for stage 4 runs with the default rules, whose ORF test
  is by length) it carries the gene when its aligned row holds an ORF passing the stage 4 ORF rule, the definition
  stage 4 uses for ancestors. Step 1a refuses a stage 4 options file without `orf_test` = `length`.
- **Unknown genomes** (`UNKNOWN_PRESENCE`, step 5 `--unknown`, a TSV `species`, `geneID`, `unknown_in`; with
  `PRESENCE=orf` the table of step 1a by default). A genome with no usable call for a gene takes whichever state gives
  the fewest deaths for the origin, so missing data never adds a death; ties put the death on the higher branch, as
  plain Dollo does. Without unknown genomes the result is plain Dollo, which step 5 checks on every origin. Both rules
  come from the T2T mammal analysis (2026-11-13).
- **Rates** are counts divided by branch length, written for the ingroup crown, every node below it and every
  species (outgroups included). Rates are NA where the length is missing or zero.
- **Gene IDs** are compared with the HOG table without isoform or `.p` suffixes (a trailing `.p`, then one trailing
  `.<number>`), within each species (`--hog-id-match normalized`). The Brassicaceae paper compared raw IDs (`exact`);
  see Validation. Step 4 always uses the normalized form to decide which loci are in no orthogroup, also when step 5
  runs with `exact`.
- **ID checks.** The run stops, rather than counting wrongly, when the IDs of two inputs differ in form: fewer than
  95% of a species' classified loci named in its coordinate table (step 1), a species none of whose genes is in the
  HOG table or some of whose genes match a HOG-table ID only after case, prefixes such as `ref|`, a species tag,
  `_` vs `.` or version suffixes are ignored (steps 4 and 5), two different IDs of a species that become one after
  the suffixes are removed, as when the trailing number is part of the gene name (`evm.model.Chr1.12`; steps 4 and
  5), fewer than 95% of a species' unclustered loci found in its protein file or several sequences for one of them
  (step 4), or a gene listed in two HOGs (step 5).
- **Overlap filter.** Loci are compared with the protein-coding genes of the GFF (a gene whose stated biotype is not
  `protein_coding` is ignored), and every overlapping gene is found, however far upstream it starts.
  - Antisense: the locus span overlaps any part of a gene on the opposite strand. Counted.
  - Same strand, coding exons do not overlap: the locus lies in an intron or a UTR of the other gene. Counted.
  - Same strand, a coding exon of the locus overlaps a coding exon of the other gene (`same_strand_cds_overlap`).
    Left out, as the locus may be a fragment or an alternative frame of that gene.
  - No overlap (intergenic, or the nearest neighbor class within 1 kb). Counted.

  - Not checkable, counted: `unknown` (no GFF gene at the locus, or its chromosome is not in the GFF) and
    `own_gene_not_found` (GFF genes overlap the locus but none carries its geneID).

  The own gene of a locus is the overlapping GFF gene that carries the locus geneID as a transcript or protein ID
  (NP_..., XP_..., ENSMUST...), as a transcript ID followed by a separator (Thhalv10029120m.v1.0), or as a whole word
  of its gene ID. It is never guessed from coordinates, because a long host gene covers the locus as well as its own
  gene; step 3 stops when more than 5% of a species' loci have no own gene (`--max-unmatched`), when more than 5% are
  `unknown` (`--max-unknown`; the GFF belongs to another assembly or names chromosomes differently), or when a `.fai`
  is given for a GFF without `##sequence-region` lines. The coding exons of
  the locus are those of the matched transcript within the locus span (those of all its matching transcripts when the
  locus isoform is not in the GFF or the gene was matched by gene ID; the log counts these loci). When a `.fai` is given to rename GFF
  chromosomes, only lengths found once in each file are renamed; genes on other chromosomes are dropped, so their
  loci become `unknown` instead of being compared with genes of another scaffold. `--exclude-overlap any` leaves out
  every overlap class, as in the paper; `unknown` and `own_gene_not_found` loci are counted under either setting.

## Inputs

| Input | From |
|---|---|
| `tree_config.json` | stage 4 |
| `*_all_species_summary.tsv`, one per focal species | stage 3 step 17 |
| `*_nodes_summary_filtered.tsv` | stage 4 final classification (`stage4_final/nodes_summary_by_species/`) |
| `*_cdngs_tblastn_*_summary.tsv` | stage 4 step 4 |
| stage 3 alignment folders (`ALIGNMENTS_DIR`) and the stage 4 method options of the run (`STAGE4_OPTIONS`), with `PRESENCE=orf` only | stage 3 / stage 4 (`WORKDIR/brassicaceae_23species_LINGUA_results/`, `WORKDIR/logs/method_options_used.json`) |
| HOG table (`N0.tsv`, or the node table that holds every species) | OrthoFinder `Phylogenetic_Hierarchical_Orthogroups/`; its species columns must be species names or aliases in `tree_config.json` |
| candidate proteins, one FASTA per species | the stage 4 input |
| branch lengths: TSV `branch`, `length` (or `node`, `age`) | a time-calibrated tree; `data/branch_lengths_brassicaceae.tsv` is the paper's (million years) |
| GFF3 per species (required): TSV `species`, `gff`, optional `fai` (`GFF_TABLE`); species with no annotation are named in `NO_GFF_SPECIES` and counted without the overlap check | the annotations; give the `.fai` only when the GFF names chromosomes differently from the genome used in stage 3 |

File-name patterns of step 1 can be changed with `--coords-glob`, `--nodes-glob` and `--tblastn-glob`. The species of
each file is read from the start of its name and must be a species or alias in `tree_config.json`.

## Running

Requirements: Python 3.10 or later (standard library only) and, for step 4, BLAST+ (`blastp`, `makeblastdb`).

```
cp 5_birth_death/config/stage5.conf.example my_stage5.conf      # edit the paths
bash 5_birth_death/run_stage5.sh my_stage5.conf                   # all steps
bash 5_birth_death/run_stage5.sh my_stage5.conf 5 5               # only step 5
python 5_birth_death/tests/test_toy.py                            # hand-checked example + the ID and overlap checks
```

`OUTDIR/STATUS` reads RUNNING, SUCCEEDED or FAILED at step N; on a failure the driver also prints the end of that
step's log. Each step logs to `OUTDIR/stepN.log`. A step refuses to overwrite its output, so a rerun goes to a new
folder (a refused rerun leaves `STATUS` as it was). `OUTDIR` must not contain spaces when step 4 runs, because
BLAST+ cannot use such paths.

## Output

- `rates_by_branch.tsv`: Branch, Length, DNG, PDNG and TDNG (DNG + PDNG) births and deaths, and their rates.
- `origins.tsv`: one row per counted origin, with its HOG, type, birth branch, death branches and number of loci.
- `gene_status.tsv`: every gene of `loci.tsv`, its locus, its HOG unit and why it is counted or not (CDS overlaps a
  same-strand gene, duplicate of an unclustered locus, PDNG paralog of a DNG, ORF node outside the window, gene type
  other than DNG or PDNG).

## Validation

Unknown genomes and step 1a (2026-11-13). With the default settings (`PRESENCE=tblastn`, no `UNKNOWN_PRESENCE`) a
full driver run on the Brassicaceae inputs gives output files byte-identical to the run of commit 710a07a. Step 1a
and `--unknown` are tested on the toy example (`test_orf_presence`, `test_unknown_genomes`); step 1a has not yet been
run on the Brassicaceae data.

Run on the Brassicaceae inputs of the published analysis (2026-11-12), each step against the output of the
script it replaces:

- Step 1 equals stage 4 step 9 (`nodes_by_coordinates.py`) on the same inputs, 30,918 of 30,918 rows, except that
  an empty Synteny_Node is written as NA.
- Step 2 equals `anc_nodes_estimates.py`, 25,821 of 25,821 loci with the same members, node and type.
- Step 3 as first written equaled `orientation_all_species.py`, 26,124 of 26,124 rows with the same class. The
  current step 3 differs on purpose (2026-11-12): the original missed loci inside long genes and, on tied overlaps,
  took a neighboring candidate for the locus's own gene, so it counted overlaps between candidates that it meant to
  ignore. Comparing whole spans, 227 genes are newly excluded (220 of them *Theobroma cacao* loci inside gene
  models of 620 and 835 kb) and 129 no longer are (126 overlap only other candidates, 3 *Camelina sativa* loci
  overlap only non-coding genes): 2,004 excluded instead of 1,906 when every overlap is left out. Under the current
  rule only the 51 same-strand CDS overlaps are excluded.
- Step 4 equals `unclustered_paralog_blast.py` (protein files byte-identical in all 23 species) and
  `build_unclustered_families.py` (17,822 of 17,822 rows; 243 duplicates).
- Step 5 with `--hog-id-match exact --hog-span merge --exclude-overlap any` and the paper's overlap file equals the
  paper's rate table
  (`dng_pdng_rates_by_branch_no_overlap_unc_dedup_noDNGparalog.tsv`) in all 504 cells of its 42 branches:
  1,427 DNG and 17,783 PDNG births, 474 DNG and 1,470 PDNG deaths.
- Births and deaths with `normalized` matching and the current overlap rule (51 genes excluded), 2026-11-12:

  | `--hog-span` | DNG births | PDNG births | DNG deaths | PDNG deaths |
  |---|---|---|---|---|
  | merge (the paper) | 1,632 | 18,758 | 546 | 1,529 |
  | split (default) | 1,700 | 19,429 | 556 | 1,563 |

  These counts use the loci of the paper (26,124 genes). A full driver run from the local stage 4 tBLASTn tables
  (`cdngs_tblastn_results_Apr-2026`) gives 26,078 genes in 25,776 loci, and 1,700 DNG and 19,386 PDNG births, 563 DNG
  and 1,573 PDNG deaths (2026-11-13); the paper's step 9 read the tBLASTn tables of the strict run on the server.
  With the split, 56 of 1,632 DNG HOG units become 124 origins and 465 of 18,758 PDNG units become 1,136.
  N0.HOG0019995 (four genes at four alignment positions) becomes three PDNG origins, at N5 (*Microthlaspi
  erraticum* and *Brassica oleracea*), N9 (*Arabidopsis lyrata*) and *Camelina neglecta*.
- Same-strand overlaps (2026-11-12). In the Brassicaceae, 705 loci overlap a same-strand gene and 51 of them overlap
  its coding exons (28 *Theobroma cacao*, 11 *Arabis alpina*, 6 *Microthlaspi erraticum*, 6 in three other species).
  The first version renamed *Raphanus sativus* scaffolds by length although 1,250 GFF scaffolds share a length; its 3
  *R. sativus* exclusions were overlaps with genes of a neighboring scaffold, and those loci are now `unknown`.
- Independent check (`nee_revision_stage4/scripts/audit_stage5_overlap_20261112.py` of the Brassicaceae analysis:
  own transcript by ID, a scan of every gene, the candidate transcript's CDS): the class agrees for 26,044 of 26,044
  Brassicaceae loci it could resolve (34 could not be resolved, 29 of them on *R. sativus* scaffolds of shared
  length) and 4,568 of 4,568 mammal loci; every locus span equals the CDS span of its transcript.
  In the T2T mammals (stage 4 of 2026-11-12, 4,568 loci), 417 loci overlap a same-strand gene and 29 overlap its
  coding exons, 23 of them in *Mus musculus*, whose annotation combines Ensembl genes lifted over to the T2T genome
  with a second gene prediction (17 of the 23 are predicted genes overlapping Ensembl genes). The goat and sheep
  annotations have almost no overlapping genes (0 and 14 genome-wide), so no locus there overlaps a gene.
- With the default `normalized` matching, the loci of 633 genes are assigned to a HOG instead of counting as origins
  of their own: *Capsella rubella* 202, *Brassica rapa* 143, *Boechera stricta* 128, *Thlaspi arvense* 81, *Capsella
  grandiflora* 77, and 2 *Arabidopsis lyrata* genes that share a locus with such genes. The paper missed these
  because the IDs of those five species carry `.p` in the HOG table only. With the paper's overlap file,
  `--exclude-overlap any` and `--hog-span merge`, births become 1,389 DNG and 17,391 PDNG (paper 1,427 and 17,783).
- Step 2 joins two genes when their intervals overlap or touch (one ending where the other starts), as the paper's
  script did; requiring a shared base gives identical loci on the Brassicaceae inputs, so the original test was kept.

## License

MIT, see `../LICENSE`.
