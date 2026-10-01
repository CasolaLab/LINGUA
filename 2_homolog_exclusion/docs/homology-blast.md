# homology-blast

Searches each candidate protein with BLAST+ against a database you built, and removes candidates that have a homolog
there. It covers two of the screens:

- **Against a protein database of other species** (`blastp`): a young gene should not match genes from other species.
- **Against organelle genomes** (`tblastn`, a nucleotide database): a hit means a leaked organellar gene posing as a
  lineage-specific gene. Organelle genes are highly conserved and inherited from one parent, so any hit counts.

## Everything runs on your machine

BLAST+ runs locally. The module never uses the network: it does not run BLAST+ with `-remote`, and it fetches nothing.

## Before you start

You install everything and build every database; the repository ships no programs or data (see `docs/setup.md`).

- **BLAST+** (`blastp`, `tblastn`, `blastdbcmd`, `makeblastdb`). Tested with 2.17.0.
- **A database: a species-tagged protein FASTA** (`make-db` builds one; see `docs/databases.md`). Give the FASTA to `--db`
  and the module builds a BLAST index from it the first time (`makeblastdb -parse_seqids`, local), keeps it in the run's
  output folder or in `--index-dir`, and reuses it until the FASTA changes. The index is a disposable cache. A database you
  already formatted can be given by its prefix instead; it must have been built with `-parse_seqids` (the module refuses one
  built without it, because BLAST+ then reports the sequences as `gnl|BL_ORD_ID|n` instead of their IDs):

```bash
makeblastdb -in proteins.faa -dbtype prot -out db/proteins -parse_seqids
makeblastdb -in organelle_genomes.fna -dbtype nucl -out db/organelle -parse_seqids
```

  (For `tblastn`, give the organelle genomes as a nucleotide FASTA and the index is built as nucleotide.)

### The database must be built with care: this is the most important rule

The tool counts **any** hit. So the database must **not contain the species you are analysing, or any close relative that
could belong in the tree.** If it does, every gene "has a homolog" (itself, or its sister species), and nothing passes.

As a safeguard, and to catch mistakes, the module checks each hit's species:

- **Every database sequence ID should carry its species after `__`**, as Stage 0 writes with `--add-species`:
  `g13875.t1__Marchantia_paleacea`. Bare IDs like `g13875.t1` repeat across species, so without the tag there is no way
  to know where a hit came from. (For databases like UniProt, use a taxon map instead: `species_source=map` and
  `taxon_map=file`, two columns, accession then species.)
- **Each run reads the database's IDs** (locally, with `blastdbcmd`) and records how many carry each species in `run.json`.
  If an analysis species is found INSIDE the database, it is named in a warning (its hits are exempt, or counted if
  `ignore_in_phylogeny=no`); `--skip-db-inventory` skips this on very large databases.
- **The analysis species** are taken from your input file names (`Physcomitrium_patens_final.faa` is
  `Physcomitrium_patens`; `--species-map` overrides this with the proper name). Hits to an analysis species, or to a species
  in `ignore_taxa`, are reported but **not counted**: the gene is `IN_PHYLOGENY_ONLY`, and it passes.
- **A hit whose species is unknown is counted**, never exempted, and says so. If none of the first 1,000 IDs in a database
  carries a species tag, the module warns that every such hit will be counted.
- `match_level` sets how species names are compared: `species` (default: the first two words, so subspecies, varieties and
  strains count as the same species: `Marchantia polymorpha subsp. ruderalis` matches `Marchantia_polymorpha`), `exact`
  (the whole name), or `genus` (also exempts hits to congeners, so use with care).

## Run it

```bash
# protein database of other species
bash scripts/homolog_exclusion_pipeline.sh homology-blast -i proteins/ -o out/ --db db/proteins.faa \
    --species-map species.tsv --threads 8

# organelle genomes: the preset turns off the analysis-species exemption and the coverage requirement
bash scripts/homolog_exclusion_pipeline.sh homology-blast -i proteins/ -o out_org/ --db db/organelle \
    --config presets/organelle_tblastn.conf --label organelle
```

`-i` is one FASTA file or a folder of per-species FASTA files. Sequences are split into chunks (`chunk_size`, longest
first), `parallel_jobs` chunks at a time.

To re-filter results you already have, skip the search: `--parse-only earlier_out/raw` (a folder of
`<species>.blast.tsv` files made with `-outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send
evalue bitscore qlen slen"`; another column layout is refused, not guessed at).

## Which hits count

All alignments of one candidate to one database sequence are merged into a single hit. The hit counts when **all** hold:

| Rule | Setting | Default | Why |
|---|---|---|---|
| E-value at most this (per alignment) | `evalue` | 1e-3 | a common filtering cutoff; BLAST's own default (10) is a search default |
| The alignments cover at least this percent of the **query** | `min_qcov` | 50 | a short local match should not count as homology; the query is your candidate |
| ... of the **target** (`blastp` only) | `min_tcov` | 0 (off) | would wrongly discard a short candidate that matches one domain of a long protein |
| Identity at least this percent | `min_pident` | empty (off) | identity depends on how diverged the species are |
| The target's species is not an analysis species or in `ignore_taxa` | `ignore_in_phylogeny`, `ignore_taxa` | yes, empty | see above |

Coverage counts every position once: alignments covering 60 and 51 residues of a 200-residue query, without overlapping,
cover 111 (55.5%). For `tblastn` the target coordinates are nucleotides, so target coverage is not reported and `min_tcov`
is refused.

Other settings, passed straight to BLAST+: `program` (`blastp` or `tblastn`), `seg` (mask low-complexity regions in the
query; default yes), `max_hsps`, `ungapped` (needs `comp_based_stats=0`; BLAST+ refuses it otherwise, and the module says
so before running), `comp_based_stats`, `dbsize` (make E-values comparable across databases of different sizes),
`max_target_seqs` (5000: BLAST+ keeps the *first* N hits it finds, not the best N, so a small value can hide hits).

Change any setting with a flag (`--min-qcov 30`), `--set key=value`, or a config file. Each run records what it used in
`run.json`, and a warning names every setting that differs from the default.

## Organelle genomes

`presets/organelle_tblastn.conf` is an example (`program=tblastn`, ungapped, `comp_based_stats=0`, `max_hsps=5`,
`min_qcov=0`, `ignore_in_phylogeny=no`): **any** hit to an organelle genome counts, including a hit to the analysis
species' own organelle genome, because that is the signal. Give the run its own label (`--label organelle`) so `combine`
shows it as a separate column.

## Several databases, and dropping genes as you go

- **Chaining (works today, for any module and any database).** Each run writes its survivors to `passed/`; give that folder
  to the next run as `-i`. A gene that hit the first database is never searched against the next.

```bash
bash scripts/homolog_exclusion_pipeline.sh homology-blast -i proteins/ -o out/1_mosses --db db/mosses --label mosses
bash scripts/homolog_exclusion_pipeline.sh homology-blast -i out/1_mosses/passed -o out/2_fungi --db db/fungi --label fungi
```

- **`--db` more than once** searches every database against **all** the genes, each as its own labelled run in its own
  subfolder (`out/<label>/`). A database that fails does not stop the others; the command exits with an error at the end
  if any failed. This is the "all at once" mode, not iterative.
- **`--sequential` (iterative filtering across databases in one command).** With several `--db`, add `--sequential`: the
  databases are searched in the order you give them, and each one searches **only the genes that passed the previous one**.

```bash
bash scripts/homolog_exclusion_pipeline.sh homology-blast -i proteins/ -o out/ --sequential \
    --db db/mosses --db db/fungi --db db/animals --db db/plants
```

  Each database gets its own labelled folder (`out/mosses/`, `out/fungi/`, ...) with the usual files, and `summary.tsv`
  shows how many genes each one received. `run.json` records the position in the chain and where the input came from.
  The genes that passed every database are in the last folder's `passed/`. **If a database fails, the chain stops** and
  the command says which databases were not searched, because genes cannot pass a database that was never searched. (In the
  all-at-once mode the databases are independent, so a failure does not stop the others.) Put the databases most likely to
  hit first: every gene removed early is a gene the later, slower searches never see.

## What you get

| File | Content |
|---|---|
| `gene_status.tsv` | one row per candidate: `HIT`, `NO_HIT`, `IN_PHYLOGENY_ONLY`, or `NOT_RUN` |
| `hits.tsv` | one row per candidate and database sequence: target, target species, E-value, bit score, identity, query and target coverage, whether it counts, and why not |
| `summary.tsv` | one row per species and a total, with the count of each status |
| `passed/<species>.faa`, `passed_ids.tsv` | the candidates that pass (`NO_HIT` and `IN_PHYLOGENY_ONLY`), with original headers |
| `run.json`, `citations.txt` | settings, BLAST+ version, the database (type, size, how many IDs carry a species tag), the analysis species; the end-of-run note |
| `raw/<species>.blast.tsv` | the raw BLAST+ results, so a run can be re-filtered with different settings |

**`NOT_RUN` is never "no hit".** A gene is `NOT_RUN` if its BLAST job failed, if the output for it never reached BLAST+'s
completion line (a killed or truncated run), or if it is missing from the results. BLAST+ writes `# BLAST processed N
queries` at the end of a finished run; a group of queries is trusted only if that line follows it and N equals the number of
queries actually seen. Output folders are protected: the tool will not overwrite an earlier run (`--resume` continues an
interrupted run, `--force` overwrites the same run).

## What has been tested

- **Logic, on mock and real results** (`tests/test_homology_blast.py`): every rule above on 14 hand-made genes (one rule
  each), the species levels (subspecies and strains), the taxon map, the organelle settings, unfinished and malformed
  output, and real BLAST+ output kept in `example_data/homology-blast/real_sample*`. The tests were checked by deliberately
  breaking the code in eight ways; each break made them fail.
- **The run path with stand-in `blastp`, `tblastn` and `blastdbcmd`** (Linux and WSL): the command line, a failed chunk giving
  `NOT_RUN`, output without the completion line not being trusted, a database of the wrong type or built without
  `-parse_seqids`, a database with untagged IDs, and several databases with one broken.
- **Real BLAST+ 2.17.0, network cut off** (2026-09-26, WSL): 200 random proteins from a real moss proteome, 5 shuffled
  controls and 40 plastid-encoded proteins, searched against a database of 8 real bryophyte proteomes (249,000 sequences)
  that deliberately included the query species. 245 genes, none `NOT_RUN`, 97 seconds. The module's hits were identical to a
  direct BLAST+ run (42,209 gene-to-target pairs, none missing). With the species check, 187 of the 200 genes hit and 13
  (whose only hits were to their own species) passed; **without it, all 200 hit and none passed**, which is the mistake the
  check exists to catch. The shuffled controls had no hits.
- **Real organelle search:** the same inputs against the real *Marchantia paleacea* chloroplast genome (NCBI NC_001319.1)
  with the preset: 37 of 40 plastid-encoded proteins hit, and 3 of 200 random nuclear proteins (1.5%). Without the preset
  (default 50% coverage, own species exempt), 1 of 200.
- **Not yet tested:** databases of the size the full pipeline uses; a search that fails part-way on real data; nucleotide
  databases of many genomes.
