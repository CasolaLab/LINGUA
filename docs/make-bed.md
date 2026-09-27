# make-bed

Given a GFF3 annotation and a list of gene, transcript or protein IDs, writes a BED of their CDS coordinates and a
summary table (exon and CDS segment counts, gene/mRNA/CDS lengths). This is general purpose: the ID list can be
`passed_all_ids.tsv` from `combine`, a full candidate set before any filtering, or any list of IDs from any source.
It needs only a GFF3 file and an ID list; no FASTA, no search tools, no network.

## Before you start

You need a GFF3 file for each species and a list of the IDs you want coordinates for. Nothing here is shipped or
downloaded; both are yours.

## Run it

```bash
# one species
bash scripts/homolog_exclusion_pipeline.sh make-bed --ids passed_all_ids.tsv --gff annotation.gff3 -o out/

# several species: one ID-list file and one GFF file per species, in two folders
bash scripts/homolog_exclusion_pipeline.sh make-bed --ids-dir ids/ --gff-dir gff/ -o out/ --species-map species.tsv
```

`--ids` accepts a plain one-ID-per-line file, or a table with a header (tab- or comma-separated); the ID column is
guessed (the first column whose name contains "gene", "protein", "transcript" or "id") unless you give `--id-column`.

In batch mode the species name is the ID file's own name, with underscores read as spaces (`Arabidopsis_thaliana.tsv`
→ `Arabidopsis thaliana`), matching the `Species` column Stage 2 already writes. Pairing an ID file to the right GFF
file works two ways:

- **`--species-map FILE`** (recommended): `species name<TAB>GFF file basename`, the same two-column format used
  elsewhere in this tool. Exact, no guessing.
- **Without a map**: a GFF is guessed by prefix from the ID file's name (the full name, its first word, and, for a
  two-word "Genus_species" name, the "Gspecies" abbreviation some annotation sources use). If more than one GFF file
  fits, or two ID files would both claim the same GFF, the run refuses rather than guessing which one you meant.

## Matching

An ID is first looked up against each `mRNA` line's own `ID`, `Name`, parent gene ID and `protein_id` attributes.
If nothing matches that way, each `CDS` line is matched directly by its own `protein_id`, `locus_tag`, `Name` or `ID`
attribute — some GFFs, especially from smaller or non-model genomes, don't carry a clean `Name` on the mRNA line that
lines up with your ID list. This fallback exists because real annotation is inconsistent between sources; it is
tested with its own mock data, since the real files this was checked against never needed it (see below).

An ID that matches nothing is written to `not_found.tsv`, never silently dropped. A run with any such ID, or any
species whose ID list or GFF could not even be read, counts as incomplete: exit status 2 (see `docs/interpreting_hits.md`).

## What you get

| File | Content |
|---|---|
| `<species>.bed` | one row per CDS segment: chrom, start, end (0-based half-open), name, score (always 0), strand, matched ID |
| `master.tsv` | one row per matched ID: species, gene ID (from the GFF's parent gene, when known), matched ID, exon count, CDS segment count, gene/mRNA/CDS length |
| `not_found.tsv` | species and ID, for every requested ID that matched nothing |
| `stats.tsv` | per species: IDs requested, IDs found, CDS rows written, status |

No table has a blank cell: a length or gene ID that could not be resolved is `NA`, never 0 (a real 0-length feature
and "not known" must not look the same). A CDS segment that appears twice in the GFF with identical coordinates (a
real annotation-export quirk) is written once in the BED and counted once in the segment count and length, not twice.

## What has been tested

- **Logic** (`tests/test_make_cds_bed.py`, 33 checks): two hand-built GFF3 files in the exact layout of the two real
  annotation sources this was checked against (below), covering the primary match path, the CDS-only fallback path
  (never exercised by the real data), a duplicated CDS line, an unmatched ID, batch-mode file pairing by guess and by
  `--species-map`, the collision refusals, and every refusal (mixed single/batch flags, an unknown `--id-column`, an
  existing output folder without `--force`). The tests were checked by deliberately breaking the code in eleven ways;
  each break made them fail.
- **Real GFF3 data, network cut off**: all 23 Brassicaceae species of a real Stage-3 project, about 2.5 GB of GFF3
  (Phytozome-style with a clean per-isoform `Name`, and MAKER-style with one `Name` shared by every isoform), each
  matched against its real `C_LSG_Iso.tsv` gene list. **79,631 of 79,631 requested IDs found, 219,111 CDS rows
  written, 0 not found, 0 warnings**, both with `--species-map` and with the file-name guess. Confirmed exact parity
  (same found-count and same CDS-row-count per species) against the original script this was rewritten from.
- **Real NCBI RefSeq and Ensembl GFF3** (*Saccharomyces cerevisiae* R64, downloaded 2026-09-27 from
  `ftp.ncbi.nlm.nih.gov` and `ftp.ensembl.org`): ID lists built from each file's own real `protein_id` accessions
  found 150 of 150 in both, through the CDS-only fallback path — neither source's `mRNA` line carries a `Name` or
  `protein_id` that lines up with protein accessions, so this is a real (not mocked) test of that path. A separate
  list of 150 real NCBI `locus_tag` values found 134; the 16 misses were checked one by one and every one is a
  non-coding gene (`gene_biotype=tRNA`, `ncRNA`, `rRNA` or `snoRNA`) with no `CDS` feature to find — correctly
  reported in `not_found.tsv`, not guessed at.
