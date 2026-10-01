================================================================================
SCRIPT: step1_mafextractor.sh
================================================================================
Description: A wrapper for mafExtractor that reads a list of target coordinates and extracts localized MAF blocks from a master MAF file.

Usage:
./step1_mafextractor.sh <species_name> <bed_coordinates_file> <maf_file> <output_dir>


Inputs:
  - <species_name>: The reference species prefix (e.g., Arabidopsis_thaliana).
  - <reference_species_coordinates_file>: A tab-separated file with headers, where columns are sequence_name (chromosome ID), start coordinate, end coordinate.
  - <maf_file>: The master raw MAF file.
  - <output_dir>: Destination folder for the extracted localized MAF files.
