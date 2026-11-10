================================================================================
SCRIPT: step1_mafextractor.sh
================================================================================
Description: A wrapper for mafExtractor that reads a list of target coordinates and extracts localized MAF blocks from a master MAF file, one file per exon. Each exon is written to a temporary file and kept as <species>.<chrom>_<start>_<stop>.txt only if mafExtractor exits successfully and the output contains at least one reference line (s <species>.). Failed exons are listed on stderr and in <output_dir>/step1_failed.tsv, and the script exits with status 1; it also stops at once if mafExtractor is not on PATH or the MAF or coordinate file cannot be read. The coordinate file is read, not modified (Windows line endings are removed while reading).

Usage:
  ./step1_mafextractor.sh <species_name> <reference_species_coordinates_file> <maf_file> <output_dir>

Inputs:
  - <species_name>: The reference species prefix as used in the MAF sequence names (e.g., Arabidopsis_thaliana).
  - <reference_species_coordinates_file>: A tab-separated file with a header line, columns sequence_name (chromosome ID), start, end (BED coordinates).
  - <maf_file>: The master reference-anchored MAF file.
  - <output_dir>: Destination folder for the extracted localized MAF files.

Notes:
  - mafExtractor treats --stop as inclusive, so each file holds end - start + 1 reference bases: one base past the BED exon end. Step 9 removes it.
  - run_stage3.sh can copy existing step 1 files instead (STEP1_SRC) and then checks that every exon has a non-empty file with reference lines.
