# domain-interpro

Runs InterProScan on each candidate protein and removes candidates that carry a known protein domain or family. A young
gene should not carry a known conserved domain, so a domain match counts against it.

## Everything runs on your machine

The module always runs InterProScan with `-dp` (no pre-calculated lookup). Without it, InterProScan sends your sequences
to EBI's web service before computing anything. With it, every result is computed locally: **nothing leaves your
machine and no internet connection is needed.** This cannot be switched off. The `applications` list (below) also needs
it: InterProScan's lookup ignores the list and returns everything it knows.

## Before you start

You install everything; the repository ships no programs or data (see `docs/setup.md`).

- **InterProScan 5.x** with the data of the analyses you will run, and **Java 11 or newer**. Tested with InterProScan
  5.78-109.0 and Java 25.
- Put `interproscan.sh` on your `PATH`, or pass `--iprscan-bin /path/to/interproscan.sh`.
- InterProScan's data is large (36.5 GB for everything in 5.78). Each analysis needs only its own data folder, so you can
  install only the analyses you will run (list them in `applications`).

## Run it

```bash
bash scripts/homolog_exclusion_pipeline.sh domain-interpro -i proteins/ -o out/ --threads 4
```

`-i` is one FASTA file or a folder of per-species FASTA files. The species is the file name without its extension and
without `_final`; `--species-map` overrides this. Sequences are split into chunks (`chunk_size`, longest first) and up to
`parallel_jobs` chunks run at once, each in its own temporary folder. `*` stop symbols are removed from the sequences
before searching, because InterProScan rejects them, and you are told how many were affected.

To re-filter results you already have, skip the search:

```bash
bash scripts/homolog_exclusion_pipeline.sh domain-interpro -i proteins/ -o out/ --parse-only earlier_results/
```

The folder holds `<species>.interproscan.tsv` (InterProScan's TSV output). See "A file that lists only matches" below for
the one extra thing this needs.

## Which matches count

The decision is made by **which analysis produced the match**, never by whether the match has an InterPro accession
(real domain matches from CDD, SFLD, PIRSF and others often have none).

| The match comes from | What it means | Effect |
|---|---|---|
| A domain or family database (Pfam, SMART, CDD, PRINTS, PIRSF, PIRSR, SFLD, Hamap, PROSITE profiles, Gene3D, SUPERFAMILY, NCBIfam, ...) | domain evidence | **counts**: the gene is `HIT` and does not pass |
| MobiDB-Lite, Coils, PROSITE patterns (`exclude_analyses`) | disorder, coiled coils, short patterns that match by chance: not domains | never counts: a gene with only these is `EXCLUDED_ONLY` and passes |
| AntiFam (`spurious_analyses`) | flags a **spurious gene model**, not a domain | never counts as a domain: a gene with only this is `SPURIOUS`, and passes unless `antifam_action=remove` |
| A domain listed in `ignore_domains` | a known lineage-specific domain for your clade | does not count: a gene with only these is `EXCLUDED_ONLY` |

If a gene has a real domain match as well as one of the ignored signals, the domain wins and the gene is `HIT`.

MobiDB-Lite is the reason for this rule. On a real test of 100 moss proteins it produced a large share of all matches;
counting it as a domain would remove genes that have no domain at all.

## Settings

Change any setting with a flag (`--max-evalue 1e-5`), `--set key=value`, or a config file. Each run records what it used
in `run.json`, and a warning names every setting that differs from the default.

| Setting | Default | Meaning |
|---|---|---|
| `applications` | the Casola lab's 12 databases plus MobiDBLite, Coils, ProSitePatterns, AntiFam | analyses to run (`-appl`). This is the lab's selection, not InterProScan's own default set. The four signals are run so they show up in the results and can be reported, not silently absent; ignoring them gives the same pass/fail answer as not running them |
| `exclude_analyses` | MobiDBLite, Coils, ProSitePatterns | matches that are not domain evidence |
| `spurious_analyses` | AntiFam | matches that flag a spurious gene model |
| `antifam_action` | `flag` | `flag`: report, gene still passes; `remove`: gene does not pass |
| `ignore_domains` | empty | Pfam or InterPro accessions or names of domains that do not count (case-insensitive) |
| `max_evalue` | empty (off) | optional E-value cut-off. It is applied only to analyses that report an E-value (Pfam, SMART, CDD, PRINTS, PIRSF, PIRSR, SFLD, AntiFam, ...). Hamap and PROSITE profiles report a score, and Coils, MobiDB-Lite and PROSITE patterns report nothing, so the cut-off never applies to them. With it off, each database's own curated thresholds decide |
| `chunk_size`, `parallel_jobs` | 1000, 4 | sequences per InterProScan job; jobs at once |

**Analysis names must be ones your InterProScan accepts.** They are not case-sensitive, but `MobiDB-Lite` with a hyphen is
rejected (use `MobiDBLite`), and some analyses (for example PANTHER) may not exist or may be deactivated in your
installation. If InterProScan rejects the list, the run stops and shows InterProScan's own message.

**Run the check before a long run.** `bash scripts/homolog_exclusion_pipeline.sh check --modules domain-interpro` confirms
that InterProScan and Java are found, that InterProScan accepts every analysis name in your list, and that each analysis
has its data folder. The last check matters: InterProScan will accept an analysis name even when its data is not
installed, and the run then fails. The name test starts InterProScan briefly (up to a minute; `--skip-analysis-test` skips
it). Use `--applications LIST` to check a different list.

## A file that lists only matches

InterProScan's TSV contains one row per match and **nothing for proteins that had no matches**. A protein absent from the
file was either searched and had no matches, or never searched, and the file cannot say which. So the module needs proof
that the run finished before it will call anything `NO_HIT`:

- **In a run,** each chunk is trusted only if InterProScan exited normally **and** wrote its `100% done` line. A chunk that
  fails or stops early gives `NOT_RUN` for its genes.
- **With `--parse-only`,** proof is InterProScan's own log next to the TSV (`<species>.interproscan.log`), the status file
  this module writes (`<species>.interproscan.status`), or you can state it with `--assume-complete`. Without any of these,
  only genes with a counted domain are decided (a match is a match); every other gene is `NOT_RUN`, with a warning.

## What you get

| File | Content |
|---|---|
| `gene_status.tsv` | one row per candidate: `HIT`, `NO_HIT`, `EXCLUDED_ONLY`, `SPURIOUS`, or `NOT_RUN` |
| `hits.tsv` | one row per match: analysis, signature, E-value or score, how much of the query it covers, whether it counts, and why not. Counted rows say `InterPro IPR...` or `no InterPro entry` |
| `summary.tsv` | one row per species and a total, with the count of each status |
| `passed/<species>.faa`, `passed_ids.tsv` | the candidates that pass, with their original headers |
| `run.json`, `citations.txt` | settings, InterProScan's version and the version of every analysis it ran; the end-of-run note |
| `raw/<species>.interproscan.tsv` and `.status` | the merged raw results and whether all chunks finished, so a run can be re-read later |

**`NOT_RUN` is never "no hit".** Output folders are protected: the tool will not overwrite an earlier run or write into a
folder it did not create (`--resume` continues an interrupted run, `--force` overwrites the same run).

## What has been tested

- **Logic, on mock and real results** (`tests/test_domain_interpro.py`): the rules in the table above, on 11 hand-made genes
  (one rule each) and on a small real InterProScan output (6 proteins, 78 rows); the requirement for proof of completion;
  refusing a wrong format (a header line, the wrong number of columns).
- **The search step with a stand-in `interproscan.sh`** (Linux and WSL): `-dp` always present, a separate temporary folder per
  job, a failed chunk giving `NOT_RUN` for its genes, exit code 0 without the `100% done` line not being trusted, an invalid
  analysis name stopping the run with InterProScan's message, and stop symbols being removed.
- **A real end-to-end run** (2026-09-26, InterProScan 5.78-109.0, Java 25, WSL, **network cut off**): 100 real moss proteins
  plus InterProScan's own 6-protein test set as two species, chunks of 50 with 2 jobs at a time. 106 genes, none `NOT_RUN`;
  the module's 538 matches on the moss proteins were identical to those of a direct InterProScan run. The whole test took
  454 seconds. This used 13 analyses: the lab's selection without Gene3D, SUPERFAMILY and NCBIfam, whose data was not
  installed for the test.
- **A finding worth knowing:** in that moss sample every protein with a match had a real domain and two had no match at all,
  so the "ignore" rules changed nothing there. They matter in the six-protein test set (one protein with only AntiFam and
  MobiDB-Lite rows; one with 22 MobiDB-Lite rows and a single real Pfam domain), which the unit tests use.
- **Not yet tested:** Gene3D, SUPERFAMILY and NCBIfam; complete proteomes and throughput at scale (start-up costs a minute or
  more per InterProScan job, so use large chunks, and expect this module to be the slowest); a search that fails part-way on
  real data. These belong on a workstation.
