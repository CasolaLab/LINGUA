# domain-cdd

Searches each candidate protein against the NCBI Conserved Domain Database (CDD) with RPS-BLAST and removes candidates
that match a known domain. A young gene should not carry a known conserved domain, so a domain hit counts against it.

## Before you start

You install everything; the repository ships no programs or data (see `docs/setup.md`).

- **`rpsblast`** from BLAST+ (`setup/environment.yml` installs it).
- **CDD, prebuilt for rpsblast.** NCBI provides it ready to use, so no formatting step is needed. The file list is in
  `setup/data_sources.txt`. Point `--db` at the database *prefix*, for example `/data/cdd/Cdd`. The archive unpacks to
  25 numbered parts (`Cdd.00.*` to `Cdd.24.*`) and a `Cdd.pal` file that ties them together.
- The archive does not say which CDD release it is. NCBI publishes that in a separate small file, `cdd.info`, in the
  same NCBI folder. If you save it next to the database, each run records the release (for example "cdd version 3.21")
  in `run.json`. Without it the release is recorded as unknown.

## Run it

```bash
bash scripts/homolog_exclusion_pipeline.sh domain-cdd -i proteins/ -o out/ --db /data/cdd/Cdd --threads 4
```

`-i` is one FASTA file or a folder of per-species FASTA files. The species is the file name without its extension and
without `_final` (so `Arabis_alpina_final.faa` is `Arabis_alpina`); `--species-map` overrides this.

To re-filter results you already have (for example a search run on a cluster), skip the search:

```bash
bash scripts/homolog_exclusion_pipeline.sh domain-cdd -i proteins/ -o out/ --parse-only earlier_results/
```

The folder must hold `<species>.rpsblast.tsv` files made with
`rpsblast -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"`.
Any other column layout is refused, not guessed at.

## What counts as a hit

An alignment line counts toward a hit when it passes the E-value. All the alignments of one candidate to one CDD model
are then merged into a single hit, and the hit counts when **all** of these hold:

| Rule | Setting | Default | Why |
|---|---|---|---|
| The E-value is at most this | `evalue` | 1e-2 | The NCBI CD-Search default |
| The alignments cover at least this percent of the **domain model** | `min_domain_cov` | 50 | A partial match to a domain model is weaker evidence of a domain |
| The alignments cover at least this percent of the **query** | `min_qcov` | 0 (off) | Optional |
| The model's source is not in `ignore_sources` | `ignore_sources` | empty | Optional |

Coverage counts every position once: two alignments that cover 45 and 46 of a 100-position model, without overlapping,
cover 91%. Any CDD hit that counts marks the candidate `HIT`, and the candidate does not pass.

Change any setting with a flag (`--min-domain-cov 70`), `--set key=value`, or a config file. Each run records what it
used in `run.json`, and a warning names every setting that differs from the default.

## What you get

In the output folder:

| File | Content |
|---|---|
| `gene_status.tsv` | one row per candidate: `HIT`, `NO_HIT`, `EXCLUDED_ONLY` (only hits from `ignore_sources`), or `NOT_RUN` |
| `hits.tsv` | one row per candidate and CDD model: accession (for example `pfam00069`), source (`PFAM`, `SMART`, `COG`, `KOG`, `CD`, `PRK`, ...), E-value, bit score, query and model coverage, whether it counts, and a note saying why not |
| `passed/<species>.faa` | the candidates that pass, one FASTA per species, with their original headers |
| `summary.tsv` | one row per species (and a total): genes in, genes passed, and how many of each status |
| `passed_ids.tsv`, `run.json`, `citations.txt` | the passing IDs; the settings and tool versions used; the end-of-run note |
| `raw/<species>.rpsblast.tsv` | the raw RPS-BLAST results, kept so the run can be re-parsed with different settings |

**Output folders are protected.** The tool will not write into a folder that already holds files it did not create, or
one that holds another module's results, and it will not silently overwrite an earlier run. For a folder holding the
same run, use `--resume` to continue an interrupted run (finished chunks are skipped) or `--force` to overwrite the
results; otherwise choose a new folder.

**`NOT_RUN` is never "no hit".** A candidate is `NOT_RUN` if its search job failed, if RPS-BLAST's output for it never
reached its completion line (a killed or truncated run), or if it is missing from the results. Such candidates do not
pass, and a warning names them, so a crashed run cannot make genes look novel.

How completion is checked: RPS-BLAST writes one `# BLAST processed N queries` line at the end of a finished run. A
group of queries is trusted only if that line follows it and N equals the number of queries actually seen. If a run was
cut short and another finished run follows in the same file, the counts disagree and none of those queries is trusted,
because the file alone cannot say which ones were lost.

## Read this before interpreting hits

- **CDD is built from many species' proteins.** If your query species was among the sequences that models were built
  from, a hit can mean the gene is already a known, annotated gene, which is exactly what this screen is meant to catch.
  The tool still counts it. The source column lets you look at where hits come from.
- **`ignore_sources` is for a deliberate choice**, for example if you decide to discount a whole category of models.
  Discounted hits are still listed in `hits.tsv`, with the reason.
- **CDD overlaps with InterProScan.** CDD is also one of InterProScan's member databases, and Pfam and SMART models
  appear in both. Running both is fine (`combine` shows where they agree); do not expect them to be independent.
- The default 50% model coverage is a choice, not a law. A short candidate can match a large part of a small domain, and
  a long candidate can match a small part of a large one. Check borderline hits in `hits.tsv` (`tcov` and `qcov`).

## What has been tested

- **Logic, on mock results** (`tests/test_domain_cdd.py`, no databases needed): the E-value and coverage rules, merging
  of several alignments, `ignore_sources`, `NOT_RUN` handling of unfinished and malformed results, output-folder
  protection, and refusing a wrong format. The tests were also checked by deliberately breaking the code in six ways;
  each break made them fail.
- **A real run on real data** (2026-09-26, WSL, BLAST+ 2.17.0, CDD 3.21): 200 proteins from a real moss proteome plus
  two well-known domain proteins (human HRAS and CDK2) and a shuffled copy of each as negative controls. Both real
  domain proteins were `HIT`, both shuffled copies `NO_HIT`, 188 of the 200 moss proteins had a counted CDD hit, and
  none were `NOT_RUN`. About 13 seconds for 200 proteins with 4 jobs of 2 threads each.
- **The same run again with the network cut off** (`unshare -rn`): identical results (188 of 200 moss proteins hit, both real
  domain proteins hit, both shuffled copies no hit, none `NOT_RUN`). The module needs no internet connection.
- **What that real run taught:** real RPS-BLAST output differs from what was first assumed (a header before every query,
  one footer at the end of the file). The first parser read every gene as `NOT_RUN`; it failed safe, and it was fixed
  and covered by `example_data/domain-cdd/real_sample/`, four real proteins with their real RPS-BLAST output.
- **Not tested on real data:** a search job that fails part-way, `--resume` after a real interruption, and full
  proteomes at scale. Those paths are tested on mock data only.
