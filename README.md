# Homolog Exclusion Pipeline

Stage 3 of the LINGUA comparative genomics framework. It takes candidate lineage-specific genes (from Stage 2) and removes those
that show evidence of not being specific to your lineage: a known protein domain, or a homolog in species outside your clade. The
genes that remain have no such evidence in the searches you ran.

It was built for the identification and evolutionary analysis of de novo genes. It is dataset-agnostic: it does not assume a
clade, and the same approach applies to plants, mammals, insects or anything else.

**Everything runs on your machine.** No module sends sequences anywhere or needs the internet (downloading reference data is a
separate, optional step). The repository ships no programs and no databases; you install and build those (`docs/setup.md`).

## What it does

Four search modules, each runnable on its own, plus tools to build databases and merge results:

| Command | What it does | Uses |
|---|---|---|
| `domain-interpro` | removes genes with a known protein domain or family | InterProScan |
| `domain-cdd` | removes genes that match a domain of NCBI's Conserved Domain Database | RPS-BLAST |
| `homology-blast` | removes genes with a homolog in a database you built (or, with `tblastn`, in organelle genomes) | BLAST+ |
| `homology-jackhmmer` | the same, more sensitive (remote homologs) | HMMER |
| `make-db` | builds the species-tagged FASTA that the homology modules search | Python only |
| `combine` | merges the runs into one per-gene table and a verdict (`PASS`, `REMOVED`, `INCOMPLETE`) | Python only |
| `db-report` (optional) | shows how well your database covers the tree of life outside your clade | NCBI taxonomy |
| `make-bed` (optional) | looks up CDS coordinates for a list of gene/protein IDs in a GFF3, for handoff to downstream steps | Python only |
| `check` | reports which tools and databases are installed (changes nothing) | |

Run them all on the same genes and merge, or chain them so each searches only what the previous one passed (`--sequential` does
this across several databases in one command). See `docs/workflows.md`.

## Quick start

```bash
conda env create -f setup/environment.yml && conda activate hep     # BLAST+, HMMER, Java, Python 3.12
bash scripts/homolog_exclusion_pipeline.sh check                     # what is installed
bash scripts/homolog_exclusion_pipeline.sh homology-blast -i candidates/ -o out/blast --db db/mosses.faa
bash scripts/homolog_exclusion_pipeline.sh combine --runs out/blast -o out/combine
```

`candidates/` is a folder with one protein FASTA per species. Full setup, including how to get InterProScan and CDD, is in
`docs/setup.md`.

## Principles

- **A search that did not finish is never "no hit".** A gene whose search failed, was cut short, or was never run is `NOT_RUN`
  (or `NOT_SEARCHED`) and does not pass. Completion is checked from each tool's own end-of-output marker. A run with any such
  gene exits with status 2 (results still written), so scripts stop instead of passing partial results on.
- **No blank cells, and no 0 that hides "not searched".** Every table uses `NA` for a missing value. A gene an earlier run removed is
  `NOT_SEARCHED` in the later columns, not `0`.
- **Species matter.** Hits to your own analysis species are reported but not counted; subspecies and strains count as the same
  species. Every run warns if an analysis species is found inside a database.
- **Standard defaults, all changeable.** Every cut-off is documented (`docs/thresholds.md`), can be changed by a flag or a config
  file, and is recorded with the tool versions in each run's `run.json`.
- **Nothing is overwritten by accident.** Output folders are protected; `--resume` continues an interrupted run.

## Documentation

| | |
|---|---|
| `docs/setup.md`, `setup/data_sources.txt` | install everything; where the data comes from |
| `docs/workflows.md` | recipes: all at once, iterative, several databases, organelle |
| `docs/databases.md` | build a database, and design one that suits your clade |
| `docs/thresholds.md` | every default, where it comes from, how far it was checked |
| `docs/interpreting_hits.md` | what the statuses mean and what a pass does not prove |
| `docs/domain-interpro.md`, `domain-cdd.md`, `homology-blast.md`, `homology-jackhmmer.md`, `combine.md`, `db-report.md`, `make-bed.md` | one page per module |
| `docs/citing.md` | what to cite |
| `REVIEW.md` | a ten-minute guide to checking the tool without reading the code |

## Testing

```bash
for t in tests/test_*.py; do python3 "$t" | tail -1; done      # every file must say ALL PASSED
bash tests/run_offline.sh                                       # the same, with the network cut off (Linux/WSL)
```

The tests need no external program: they use example output in each tool's exact format and stand-in programs. Each module was
also run for real (InterProScan 5.78, CDD 3.21, BLAST+ 2.17.0, HMMER 3.4) with the network cut off; `docs/setup.md` says what has
and has not been tested, and `REVIEW.md` lists what is still open.

## Credits and license

Developed by Adekola Owoyemi in the Protein Evolution Lab (Casola Lab), Texas A&M University, as part of dissertation work on the
identification and evolutionary analysis of de novo genes. MIT license (`LICENSE`). To cite: `CITATION.cff` and `docs/citing.md`.

> DRAFT: this wording copies Stages 1 and 2 and awaits the authors' confirmation.

Stage 1: `protein-preprocessing-isoform-pipeline`. Stage 2: `gene-stratigraphy-pipeline`.
