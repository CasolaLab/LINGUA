# db-report (optional)

Checks how well a homology database covers the tree of life **outside** your analysis clade. It is a helper for building
databases (see "Designing a database for your clade" in `docs/databases.md`); nothing else in the tool needs it, and it changes and
downloads nothing. It is guidance, not a verdict: the tiers and thresholds are heuristics, and public databases are uneven.

## What it does

1. Reads the species of your database: the `__Species_Name` tags of a FASTA (`make-db` writes them), or a plain list of names.
2. Places each species on the NCBI taxonomy.
3. Sorts each species into a **tier** by how far it lies outside your clade. The ladder is your clade's own ancestors: for a
   family, its order, then its class, phylum, kingdom, domain, and finally "outside all of these". Nearest first.
4. Reports species and sequences per tier, tiers that are `EMPTY` or `THIN`, species that hold a very large share of the database
   (`DOMINANT`), species with few sequences (`FEW_SEQUENCES`), species that could not be placed (`UNMATCHED`), and species that
   sit **inside** your analysis clade (`INSIDE_CLADE`: these do not belong in a database that tests for outside homologs).
5. Optionally writes a Newick tree of your database species and your clade.

## Get the taxonomy once (the only step that needs the internet)

```bash
bash setup/get_taxonomy.sh taxonomy/        # 80 MB download, about 550 MB unpacked
```

The dump is public data from NCBI and is not part of this repository. It changes daily; `db-report` copies the archive date and
checksum into its report. No taxonomy is needed if you write the tiers yourself (below).

## Run it

```bash
bash scripts/homolog_exclusion_pipeline.sh db-report --db db/mosses.faa db/fungi.faa \
    --clade Brassicaceae --taxonomy taxonomy/ -o report/ --tree report/tree.nwk
```

Your clade can be given four ways (exactly one): `--clade NAME`; `--clade-taxid ID` (if the name is ambiguous); `--clade-species
FILE` (a list of your analysis species, one per line; the clade is their common ancestor); or `--clade-inputs DIR` (your analysis
FASTA folder, species from the file names, with `--species-map` to give full names). If your analysis species include an
outgroup, the clade is the ancestor of all of them, and the outgroup is then inside the clade.

If your database has no `__` tags, use `--species-list FILE` (one species per line, optionally a tab and a sequence count). A
name the taxonomy cannot match is reported `UNMATCHED` (never guessed): fix it with `--taxid-map FILE` (`species<TAB>taxid`).

### Without the taxonomy: your own tiers

```bash
bash scripts/homolog_exclusion_pipeline.sh db-report --species-list names.txt --tiers my_tiers.tsv \
    --expected-tiers close,plants,fungi,animals,protists,prokaryotes -o report/
```

`my_tiers.tsv` is `species<TAB>tier`. Tiers you list in `--expected-tiers` but leave empty are reported.

## Settings that decide the warnings (all adjustable)

| Setting | Default | Meaning |
|---|---|---|
| `--ranks` | family,order,class,phylum,kingdom,domain | which ancestors of your clade form the ladder; `all` uses every ancestor |
| `--thin` | 3 | a tier with fewer species is `THIN` |
| `--max-share` | 0.25 | a species with a larger share of all sequences is `DOMINANT` |
| `--min-seqs` | 500 | a species with fewer sequences is `FEW_SEQUENCES` (a fragmentary proteome is likely) |

These numbers are proposals that have not been tested against a range of projects. Use them to look, not to decide. A tier can be
`EMPTY` simply because no species of that group has a proteome; that is a fact about your sources, not an error.

## Species matching

A database species is matched to the taxonomy by its whole name, then by its first two words (so a strain or subspecies is placed
at its species), then by its first word as a genus (reported, so you can check it). A name that matches more than one taxon is
`UNMATCHED` and never guessed. Only scientific names, synonyms and equivalent names are used: **common names ("human", "thale
cress") and authority strings are not**, because they are not names of a species.
In IDs, text after `__` counts as a species tag only if it looks like one (`Genus_species...`, a capital first letter, then
letters); an accidental `__` inside a random-looking ID is not read as a species, and is counted as untagged.

## What you get

| File | Content |
|---|---|
| `species_tiers.tsv` | per species: sequences, share of the database, taxid, how it was matched, tier, flags |
| `tier_summary.tsv` | per tier: species, sequences, share, status (`ok`, `THIN`, `EMPTY`, `INSIDE_CLADE`, `UNMATCHED`) |
| `report.json` | the clade, the taxonomy copy used (archive date and md5), the settings, and every warning |
| the `--tree` file | Newick, taxonomy-induced (no branch lengths), with the analysis clade marked `ANALYSIS_CLADE_<name>` |

No table has a blank cell: a missing value is `NA`.

## What has been tested

- **Logic** (`tests/test_db_report.py`, 38 checks): a small taxonomy in the exact NCBI file layout with answers worked out by hand:
  tiers, a strain matched through its species, ambiguous names, common names refused, the four ways to give a clade, the ladder
  order, `THIN`/`EMPTY`/`DOMINANT`/`FEW_SEQUENCES`, `INSIDE_CLADE`, the tree, your own tiers table, and each refusal. The tests were
  checked by deliberately breaking the code in twelve ways; eleven made them fail and one (a safety net that no input can reach)
  cannot be told apart from the original.
- **Real taxonomy** (NCBI dump of 2026-09-26): the 8 real bryophyte proteomes of the test database against the clade Funariaceae
  (the analysis species *Physcomitrium patens* was correctly found inside it, and the tree came out as expected), and 66 species of a
  real plants database against Brassicaceae (none inside it; the nearest outside tier is a single species). About 23 seconds
  each, most of it reading the dump.
- **Not yet tested**: a clade other than a family; a database with tens of thousands of species.
