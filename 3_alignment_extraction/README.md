# Genomes2Genes

Extends Progressive Cactus's HAL alignments into gene- and exon-level comparative
analysis — extraction, multi-exon stitching, frame validation, and codon-level
detection of start/stop loss and frameshift indels across species. Built and
validated on 21 Brassicaceae species plus 2 outgroups; developed for de novo gene
detection but applicable broadly across comparative genomics.

Progressive Cactus's primary output is a Hierarchical Alignment (HAL) file — a
graphical representation of the multi-genome alignment. Existing tools (the UCSC
Genome Browser, HalTools) let you visualize a HAL file, but beyond phylogenomic
or pangenome studies, there has not been an efficient way to use it to answer
questions about genes, synteny, and downstream biology. This pipeline works in
conjunction with a Cactus HAL output (and the existing HAL/MAF toolkits) to
convert and analyze alignment results in a way that generalizes to any research
question.

## What the pipeline does

Starting from a Progressive Cactus HAL file and a reference species' exon/
transcript coordinates (BED), the pipeline:

1. Converts the HAL alignment to MAF and per-species FASTA (reference-free bias
   is preserved — the reference species is only used to pull from the overall
   alignment, not to constrain it).
2. Extracts localized MAF blocks at the coordinates of interest, and stops if an
   extraction fails or comes back empty instead of passing empty files on.
3. Audits species presence and missing base pairs between alignment blocks
   (a byproduct of MAF's coordinate system pulling non-contiguous sequence from
   the graphical HAL alignment).
4. Handles species/sequences showing atypical genomic jumps (>50kb on the same
   chromosome, a jump between chromosomes, or a change of strand) — likely alignment
   artifacts, though possibly genome assembly or annotation errors. By default each
   species keeps its dominant collinear chain of blocks in each exon and loses only
   the other blocks; the strict alternative removes the species from the exon.
5. Converts each exon to an aligned FASTA by concatenating the MAF blocks and
   reverse-complements sequences per the reference strand (MAF tools assume the
   positive strand by default). Bases a species skips between blocks are unaligned
   sequence and are not added, so no genome FASTAs are needed.
6. Stitches individual exon blocks into full-length transcripts using the
   reference species' transcript IDs. Each species and ancestral node then keeps
   one row, built from its dominant collinear chain of exons (one contig and
   strand, in reference exon order, within 50 kb); its other exons become gaps,
   and every stretch where a species has no aligned sequence is recorded in a
   missing-data mask.
7. Removes the one extra reference base per exon that the extraction adds, at its
   position from the BED, so every species is cut at the same correct column, and
   validates sequence length and identity against the reference nucleotide FASTA
   (and the BED exon lengths against that FASTA).
8. Analyzes the validated, per-species alignment to locate start codons (ATG) and
   stop codons (TAA/TAG/TGA), read in each species' own frame and reported in
   reference coordinates, identify frameshift-length indels (in practice deletions
   relative to the reference, since the trimmed reference has no gaps; gaps inside
   the missing-data mask are not reported, and deletions overlapping it are trimmed
   to their unmasked part), and extract the first/last three bases of each sequence.
9. Aggregates every species' results, relative to the reference, into a master
   summary table.

Every step produces its own CSV/TSV summary; the final step compiles a master
summary across all species. Rows are written in a fixed (sorted) order, so a rerun
on the same inputs gives identical files (except the per-exon step 2 gap tables, which no later step reads). Full usage and arguments for each step are
documented individually in [`docs/`](docs/); the step-by-step technical detail
lives in [`docs/Genomes2Genes_README.txt`](docs/Genomes2Genes_README.txt).

| Step | Script | Purpose |
|------|--------|---------|
| all | `run_stage3.sh` | Run steps 1-17 for one reference species, with checks, a status file and resume |
| 1 | `step1_mafextractor.sh` | Extract localized MAF blocks at target coordinates; stop on a failed or empty extraction |
| 2 | `step2_missingbp_summaries.sh` | Summarize species presence / missing bp between blocks |
| 3 | `step3_dominant_chain.py` | Keep each species' dominant collinear chain of blocks per exon (default) |
| 3 (alt.) | `step3_filtering_species.py` | Remove species with atypical genomic jumps from the exon (`STEP3=strict`) |
| 4 | `step4_reverse_complement_add_missingbp.sh` | MAF → aligned FASTA, strand-correct (no gap filling) |
| 5 | `step5_validate_reverse_complement.py` | Validate strand orientation via local alignment |
| 6 | `step6_update_summary_files.py` | Sync FASTA-extraction summary with the BED file |
| 7 | `step7_first_stage_validation.py` | Nucleotide-count integrity check across raw/filtered/extracted files |
| 8 | `step8_stitching_exons.py` | Stitch exon blocks into full transcripts (one row per species.chromosome) |
| 8c | `step8c_one_collinear_row.py` | One row per species from its dominant collinear exon chain; writes the missing-data mask (from LINGUA) |
| 9 | `step9_trimming_duplicate_bases.py` | Remove the extra base per exon at its BED position; keep only files whose reference equals the CDS |
| 10 | `step10_second_stage_validation.py` | Confirm final sequences match reference length/content |
| 11 | `step11_update_stitched_exons_summary.py` | Compare BED, FNA and alignment lengths (`expected_bed`, `expected_fna`, `actual_nt`, `bed_fna_check`) |
| 12 | `step12_identify_start_codons.sh` | Locate start codons, reported in reference coordinates |
| 13 | `step13_identify_stop_codons.sh` | Locate stop codons, reported in reference coordinates |
| 14 | `step14_identify_indels.sh` | Identify frameshift-length indels relative to the reference, skipping or trimming gaps in the step 8c mask (from LINGUA) |
| 15 | `step15_identify_first_and_last_three_bases.sh` | Extract first/last 3bp per species |
| 16 | `step16_extract_start_and_stops.py` | Directory-wide start/stop motif extraction and validation |
| 17 | `step17_combine_all_species_info.py` | Per-species genomic span of the exons kept (step 8c mask applied), typed living or ancestor from the tree |

`clean_nucleotide_fasta_headers.py` is a required pre-processing step for the
reference nucleotide FASTA before running step 9 (see that step's README);
`run_stage3.sh` runs it automatically.

## Requirements

- HAL tools, as bundled with Progressive Cactus (`hal2maf`), to convert the HAL alignment to MAF
- mafTools (`mafExtractor`), unless existing step 1 files are reused (`STEP1_SRC`)
- GNU parallel (step 4)
- Python 3 with Biopython and pandas
- samtools only for `scripts/build_cds_from_genome.py` (building the reference CDS FASTA from a genome)

## Running

`run_stage3.sh` runs every step in order for one reference species:

```
bash run_stage3.sh <workdir> <reference_species> <bed> <cds_fasta> <maf> <tree_newick>
```

- `<reference_species>` is the species prefix used in the MAF sequence names (e.g. `Arabidopsis_thaliana`).
- `<bed>`: tab-separated chrom, start (0-based), end, name, score, strand, transcript ID; a header line is optional.
- `<cds_fasta>`: the reference CDS, named by the BED transcript IDs.
- Existing step 1 extractions can be reused with `STEP1_SRC=<folder>` (mafExtractor is then not needed and <maf> can be replaced with -);
  `ALLOW_MISSING=1` continues when some exons have no step 1 file.
- `<bed>` must be 0-based half-open: the step 8c mask is in BED coordinates, and the run stops before step 14
  if step 9 detects a 1-based BED. The workdir must not hold any other `.bed` file.
- `STEP3=dominant|strict` chooses the step 3 behavior (default dominant); `MAXSPAN` (default 50000) is the
  largest span of the exon chain step 8c keeps per species. `ONE_ROW` is no longer used.
- Steps 12-17 see only deletions relative to the reference; insertions (including 3-bp insertions forming a
  stop codon) are not detected. In the Capra hircus test run (2026-11-10) the reference rows had no gaps already
  in the step 1 MAFs (13,071 reference rows, 0 gap characters), so insertions never reached the alignments;
  step 9 would also remove any column where the reference has a gap.
- The run checks its inputs and tools before starting, writes `STATUS` (`RUNNING`, then `SUCCEEDED` or
  `FAILED at <step>`) and `PROVENANCE.txt` in the workdir, and can be resumed from any step with
  `START_STEP=<n>` after moving that step's outputs aside.

Full options and outputs: [`docs/run_stage3_README.txt`](docs/run_stage3_README.txt).

Each step can also be run on its own, in order; arguments for each are given in its README in `docs/`. For example,
step 1 is

```
bash scripts/step1_mafextractor.sh <species_name> <reference_species_coordinates_file> <maf_file> <output_dir>
```


## Credits

Developed by Sierra C. Sanders in the Casola Lab, Texas A&M University. Step 8c (`step8c_one_collinear_row.py`)
and the masked step 14 come from LINGUA stage 3 ([CasolaLab/LINGUA](https://github.com/CasolaLab/LINGUA),
`3_alignment_extraction/`, commit 81849bc). MIT license (`../LICENSE`).
