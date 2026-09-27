# homology-jackhmmer

Searches each candidate protein with jackhmmer (HMMER 3) against a database you built, and removes candidates that have a
homolog there. jackhmmer builds a profile from each candidate and searches again (iterations), so it finds remote relatives
that `homology-blast` misses. It is slower; use it as the sensitive screen, usually after the fast one has removed the easy
cases.

## Everything runs on your machine

HMMER runs locally. The module never uses the network and fetches nothing.

## Before you start

You install everything and build every database; the repository ships no programs or data (see `docs/setup.md`).

- **HMMER 3** (`jackhmmer`). Tested with 3.4.
- **A database: a species-tagged protein FASTA.** `make-db` builds one (see `docs/databases.md`). jackhmmer reads the FASTA
  directly, so there is no index to make. A BLAST database prefix cannot be used, and the module says so.

The database rules are the same as for `homology-blast` (see `docs/homology-blast.md`): it must not contain the species you
analyse or close relatives; every ID should end `__Species_Name`; hits from the analysis species and from `ignore_taxa` are
reported but not counted (`IN_PHYLOGENY_ONLY`); a hit of unknown species is counted; `match_level` sets how names compare
(`species` by default, so subspecies and strains are the same species). Each run reads the database's IDs, records how many
carry each species in `run.json`, and warns by name when an analysis species is inside the database.

## Run it

```bash
bash scripts/homolog_exclusion_pipeline.sh homology-jackhmmer -i proteins/ -o out/ --db db/mosses.faa --parallel-jobs 8
```

`-i` is one FASTA file or a folder of per-species FASTA files. Several databases work as in `homology-blast`: `--db` more than
once gives one labelled run each; add `--sequential` and each database searches only the genes that passed the previous one
(the chain stops if a database fails).

To re-filter results you already have, skip the search: `--parse-only earlier_out/raw` (a folder of
`<species>.jackhmmer.domtbl`, made with `--domtblout`). A file you made yourself must end with HMMER's `# [ok]` line.

## Which hits count

A hit is one candidate against one database sequence (jackhmmer reports one row per domain; the rows are merged). It counts
when **all** hold:

| Rule | Setting | Default | Why |
|---|---|---|---|
| The whole-sequence E-value is at most this | `evalue` | 1e-3 | a common filtering cutoff; HMMER's own default (10) is a search default |
| The domains whose own E-value also passes cover at least this percent of the **query** | `min_qcov` | 50 | a short local match should not count as homology |
| ... of the **target** | `min_tcov` | 0 (off) | would wrongly discard a short candidate that matches one domain of a long protein |
| The target's species is not an analysis species or in `ignore_taxa` | `ignore_in_phylogeny`, `ignore_taxa` | yes, empty | see `homology-blast` |

Coverage counts every position once (overlapping domains are merged, not added), and a weak domain of a strong hit does not
add coverage. Other settings: `iterations` (`-N`, default 3), `incE` (the inclusion E-value that decides which hits build the
profile, default 1e-3), `dbsize` (`-Z`: fixes the database size so E-values compare across databases; blank = actual size).
Change any setting with a flag, `--set key=value`, or a config file; `run.json` records what was used.

The tables record the **final** iteration only. jackhmmer's own text output says in which round each hit first appeared, but it
is very large (4,682 lines for 3 queries even without alignments), so the module discards it and does not report the round.

## Speed

jackhmmer's cost per query depends on how many hits it finds and how large its profile grows, not on the protein's length. So
the default is one query per job (`chunk_size=1`), which lets the slowest single query, not the slowest group of them, set the
run time. In one test with chunks of 10, 19 of 20 chunks finished in 18 minutes and the last took 37 minutes, with 7 of 8 CPUs
idle. One query per job costs about 7% more CPU (the database is read again for each query); the wall-clock gain was **not
measured**, because the test machine's clock jumped during that run. `parallel_jobs` is the number of jobs at once and `--threads` the CPUs
each job uses (as in the other modules; the default is 1, so a run uses `parallel_jobs` CPUs).

## What you get

| File | Content |
|---|---|
| `gene_status.tsv` | one row per candidate: `HIT`, `NO_HIT`, `IN_PHYLOGENY_ONLY`, or `NOT_RUN` |
| `hits.tsv` | one row per candidate and database sequence: target, target species, whole-sequence E-value, bit score, query and target coverage, whether it counts, and why not (identity is `NA`: jackhmmer does not report it) |
| `summary.tsv` | one row per species and a total, with the count of each status |
| `passed/<species>.faa`, `passed_ids.tsv` | the candidates that pass, with original headers |
| `run.json`, `citations.txt` | settings, HMMER version, the database and its species inventory, the analysis species; the end-of-run note |
| `raw/<species>.jackhmmer.domtbl`, `raw/<species>.jackhmmer.searched` | the raw results of the jobs that finished, and the list of queries they covered, so a run can be re-read |

**`NOT_RUN` is never "no hit".** jackhmmer writes nothing at all for a query without hits, so the only proof that such a query
was searched is that HMMER finished. HMMER writes `# [ok]` as the last line of a table file only when the search completed; a
file without it is not trusted, even if jackhmmer exited with code 0. A job that fails or does not finish makes every query in
it `NOT_RUN`, and those queries are left out of the `.searched` list so re-reading the raw files keeps them `NOT_RUN`. A file
you made by hand has no such list: if it ends with `# [ok]` it is taken to cover every gene of the input, otherwise none.
Output folders are protected as in the other modules (`--resume`, `--force`).

## What has been tested

- **Logic, on mock and real results** (`tests/test_homology_jackhmmer.py`, 51 checks): every rule above on 14 hand-made genes
  written in the exact `--domtblout` format of HMMER 3.4, the settings, an unfinished file, a missing file, and real HMMER
  output kept in `example_data/homology-jackhmmer/real_sample` (3 real proteins and 2 shuffled ones; at most 8 rows kept per
  query). The tests were checked by deliberately breaking the code in nine ways; each break made them fail. (Two of the nine
  were missed at first; the mock genes for overlapping domains and for a failing whole-sequence E-value were added because of
  that.)
- **The run path with a stand-in `jackhmmer`** (Linux and WSL): the command line, a failed job giving `NOT_RUN`, output without
  `# [ok]` not being trusted, a BLAST prefix refused, an analysis species inside the database, and several databases with
  `--sequential` (one broken).
- **Real HMMER 3.4, network cut off in the test suite** (2026-09-26, WSL): 200 random proteins from a real moss proteome against
  7 real bryophyte proteomes (227,344 sequences; the query species left out with `make-db`). None `NOT_RUN`; 188 genes hit and
  12 passed, and the statuses were identical with chunks of 10 and chunks of 1. Compared with `blastp` on the same genes and
  database: 187 hit in both, the same 12 passed in both, and one gene that BLAST passed has a jackhmmer hit (jackhmmer is
  the more sensitive). About 150 CPU-minutes of work in total.
- **Not yet tested:** databases and query sets of the size the full pipeline uses; a real search that fails part-way; the
  wall-clock effect of `chunk_size=1`.
