# Progressive Cactus Gene Analysis Pipeline

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
2. Extracts localized MAF blocks at the coordinates of interest.
3. Audits species presence and missing base pairs between alignment blocks
   (a byproduct of MAF's coordinate system pulling non-contiguous sequence from
   the graphical HAL alignment).
4. Filters out species/sequences showing atypical genomic jumps (>50kb on the
   same chromosome, or a jump between chromosomes) — likely alignment artifacts,
   though possibly genome assembly or annotation errors.
5. Fills missing bases from each species' genome FASTA and reverse-complements
   sequences per the reference strand (MAF tools assume the positive strand by
   default).
6. Stitches individual exon blocks into full-length transcripts using the
   reference species' transcript IDs.
7. Trims padding and intron bleed-over introduced by the whole-genome aligner,
   validating sequence length, identity, and reading-frame periodicity against
   the reference nucleotide FASTA.
8. Analyzes the validated, per-species alignment to locate in-frame start
   codons (ATG) and stop codons (TAA/TAG/TGA), identify frameshift-causing
   indels — including specific 3bp insertions that create a premature stop
   codon — and extract the first/last three bases of each sequence.
9. Aggregates every species' results, relative to the reference, into a master
   summary table.

Every step produces its own CSV/TSV summary; the final step compiles a master
summary across all species. Full usage and arguments for each step are
documented individually in [`docs/`](docs/); the step-by-step technical detail
lives in [`docs/progressive_cactus_gene_analysis_pipeline_README.txt`](docs/progressive_cactus_gene_analysis_pipeline_README.txt).

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `step1_mafextractor.sh` | Extract localized MAF blocks at target coordinates |
| 2 | `step2_missingbp_summaries.sh` | Summarize species presence / missing bp between blocks |
| 3 | `step3_filtering_species.py` | Drop species with atypical genomic jumps |
| 4 | `step4_reverse_complement_add_missingbp.sh` | MAF → aligned FASTA, gap-fill, strand-correct |
| 5 | `step5_validate_reverse_complement.py` | Validate strand orientation via local alignment |
| 6 | `step6_update_summary_files.py` | Sync FASTA-extraction summary with the BED file |
| 7 | `step7_first_stage_validation.py` | Nucleotide-count integrity check across raw/filtered/extracted files |
| 8 | `step8_stitching_exons.sh` | Stitch exon blocks into full transcripts |
| 9 | `step9_trimming_duplicate_bases.py` | Trim padding/intron bleed-over, enforce reading frame |
| 10 | `step10_second_stage_validation.py` | Confirm final sequences match reference length/content |
| 11 | `step11_update_stitched_exons_summary.py` | Compare expected vs. actual stitched-exon lengths |
| 12 | `step12_identify_start_codons.sh` | Locate in-frame start codons |
| 13 | `step13_identify_stop_codons.sh` | Locate in-frame stop codons |
| 14 | `step14_identify_indels.sh` | Identify frameshift indels, incl. premature-stop-causing 3bp insertions |
| 15 | `step15_identify_first_and_last_three_bases.sh` | Extract first/last 3bp per species |
| 16 | `step16_extract_start_and_stops.py` | Directory-wide start/stop motif extraction and validation |
| 17 | `step17_combine_all_species_info.py` | Final aggregation across all species into a master TSV |

`clean_nucleotide_fasta_headers.py` is a required pre-processing step for the
reference nucleotide FASTA before running step 9 (see that step's README).

## Development and test case

The pipeline was built and validated using 21 species from the Brassicaceae
plant family plus 2 outgroups, chosen specifically to stress-test the pipeline
against known genomic hurdles in plants (fast divergence, poor synteny) as it
was assembled.

## Status

Actively maintained, and expected to keep evolving as new tools and genomic
data become available.
