# combine

Merges the results of several runs (`domain-interpro`, `domain-cdd`, `homology-blast`, `homology-jackhmmer`, in any mix) into one
per-gene table with a final verdict, and lists the genes that passed everything. It reads only the folders the modules wrote; it
searches nothing and needs no other program.

## Run it

```bash
bash scripts/homolog_exclusion_pipeline.sh combine \
    --runs out/domain-cdd out/domain-interpro out/mosses out/jackhmmer -o out/combine
```

`--runs` lists run folders **in the order the runs were made**. That order is the chain order for iterative filtering (see below);
for runs that each searched every gene it does not matter. If you give a folder that contains run folders (for example the
output of `--sequential`, which makes `out/<database>/` folders), they are read in chain order: `--sequential` runs by their recorded
position, all others by the date they were made (a chained run is always made after the run it follows). When in doubt, list the
runs yourself. Each run needs its own
`--label` (the module name is the default), because the label is the column name in the table.

## The verdict

For every gene of every species:

| Verdict | Meaning |
|---|---|
| `REMOVED` | at least one run removed it (a counted hit) |
| `INCOMPLETE` | no run removed it, but at least one run did not finish it (`NOT_RUN`) or never searched it and nothing explains why |
| `PASS` | every run searched it and none removed it |

`INCOMPLETE` is deliberate. A gene that was not searched cannot be said to have no homolog, so it is never counted as a pass. A
gene that some run removed is `REMOVED` even if another run failed on it: one counted hit is enough.

## Iterative filtering, and `NOT_SEARCHED`

When runs are chained (each run searches only what the previous one passed), a gene removed early never reaches the later runs.
In the later columns it is written `NOT_SEARCHED`, and `dropped_by` names the run that stopped it. It is **not** written as `0`:
a `0` would say "searched, nothing found", which is false. This is the case the lab's earlier merge scripts hid, by filling every
missing value with 0.

A gene missing from a later run for a reason no earlier run explains (for example, the runs were given different input files) is
also `NOT_SEARCHED`, with `dropped_by` = `NA`, a warning, and the verdict `INCOMPLETE`.

## What you get

| File | Content |
|---|---|
| `gene_matrix.tsv` | one row per gene: one status column per run (`HIT`, `NO_HIT`, `EXCLUDED_ONLY`, `SPURIOUS`, `IN_PHYLOGENY_ONLY`, `NOT_RUN`, `NOT_SEARCHED`), then `verdict`, `removed_by` (runs joined by `;`), `dropped_by`, `n_runs_removed`, `n_runs_searched`, `n_hit_species`, `n_hit_genera` |
| `summary.tsv` | genes per species and in total by verdict |
| `run_summary.tsv` | per run: genes received, removed, passed, `NOT_RUN`, not searched, and `removed_first` (removed here and by no earlier run: how much each run adds beyond the ones before it) |
| `overlaps.tsv` | how many removed genes were removed by each exact combination of runs |
| `passed_all_ids.tsv`, `passed_all/<species>.faa` | the `PASS` genes, with their original headers |
| `gene_matrix_binary.tsv` | only with `--binary`: 1 = the run removed the gene, 0 = it searched the gene and did not, and `NA` where it did not search it or did not finish |
| `run.json`, `citations.txt` | the runs read, their order, warnings, and the tools used (for citing) |

`n_hit_species` and `n_hit_genera` count the distinct species and genera among the **counted** hits of the homology runs, a
measure of how widely a gene has homologs: one species and forty species are different biology. They are `NA` (not 0) when no
homology run searched the gene. Species are compared as the modules compare them (subspecies and strains are the same species).
Domain runs (`domain-interpro`, `domain-cdd`) do not add to these counts.

**No blank cells, ever.** A value that does not exist is `NA`. `--na X` changes that placeholder in the tables, for scripts that
want a number (for example `--na -1`); a blank is refused. In `gene_matrix_binary.tsv` do not treat `NA` as 0.

## Reading `removed_first` and `overlaps.tsv`

- With runs that each searched every gene, `overlaps.tsv` shows which searches agree and which found genes the others did not.
- With a chain, later runs never see earlier removals, so overlaps are mostly empty; read `run_summary.tsv` instead: `removed_first`
  is the marginal gain of each step, which is how you judge whether an extra database is worth its time (see "Checking that it is
  good enough" in `docs/databases.md`).

## Exit status, and what goes on to the next stage

`combine` exits 2 when any gene is `INCOMPLETE` (the tables are still written), 0 when none is, and 1 on an error. Only the `PASS`
genes (`passed_all_ids.tsv`, `passed_all/<species>.faa`) are candidates for confirmation in a later analysis. `INCOMPLETE` genes must
be searched again first; they are not passes and not removals.

## Refusals

`combine` stops, with a message, if a folder is not a run folder or has a file missing; if two runs have the same label; if a run's
own files disagree (a gene `NO_HIT` in `gene_status.tsv` but absent from `passed_ids.tsv`); if the sequence of a passed gene cannot
be found in the runs' `passed/` folders; or if the output folder already holds results (`--force` overwrites).

## What has been tested

- **Logic** (`tests/test_combine.py`, 45 checks): three mock runs in the exact files the modules write (10 genes, one rule each,
  answers worked out by hand), every rule above, the binary table and `--na`, ordering from a folder (by chain position, then by the date each run was made), and each refusal. The tests
  were checked by deliberately breaking the code in eleven ways; each break made them fail (two more tests were added for the folder-order fix).
- **Real runs** (network cut off): 200 real moss proteins through the whole chain, `domain-cdd` (real CDD 3.21), then `homology-blast` on
  its 38 survivors, then `homology-jackhmmer` on that step's 13 survivors, merged with `combine`: 12 `PASS`, 188 `REMOVED`, 0
  `INCOMPLETE`; the genes each earlier step removed are `NOT_SEARCHED` in the later columns (162 in the second, 187 in the third),
  and the whole chain took 51 seconds. The same 12 genes pass when `jackhmmer` searches all 200 genes (about 37 minutes) and when
  `blastp` and then `jackhmmer` are chained. A folder of run folders gave the same table as listing the runs in order. Trimmed run
  folders are in `example_data/combine/real_sample`.
- **Not yet tested**: real `domain-interpro` run folders in a merge, and more than about 200 genes.
