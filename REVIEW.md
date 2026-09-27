# Review guide

For a person who wants to check this tool without reading every line of code. About ten minutes.

**Status:** every module is built and tested: the shared core, `domain-interpro`, `domain-cdd`, `homology-blast`,
`homology-jackhmmer`, `make-db`, `combine`, the optional `db-report`, and the optional, general-purpose `make-bed`. Each search
module was also run for real on real data with the network cut off. What is still open is in section 4. This guide has one
section per part; read the ones you care about.

## 1. Run the tests (2 minutes)

```bash
python3 tests/test_core.py
python3 tests/test_domain_cdd.py
python3 tests/test_domain_interpro.py
python3 tests/test_homology_blast.py
python3 tests/test_homology_jackhmmer.py
python3 tests/test_make_homology_db.py
python3 tests/test_combine.py
python3 tests/test_db_report.py
python3 tests/test_make_cds_bed.py
bash tests/run_offline.sh        # the same tests with the network cut off (Linux and WSL)
```

They need no external tools and no data. Each line printed is one check. The last line of each must say `ALL PASSED`.
The tests use small made-up FASTA, config and search-result files and write only to a temporary folder.

## 2. Read the settings (5 minutes): the part only a scientist can judge

`defaults/*.conf` hold every cut-off and rule, each with a short comment saying why. This is where a wrong biological
choice would show. Look at:

| File | What to check |
|---|---|
| `defaults/domain-interpro.conf` | which InterProScan databases run; which signals do not count as a domain |
| `defaults/domain-cdd.conf` | E-value and how much of the domain model must be covered |
| `defaults/homology-blast.conf` | E-value; how much of the query must be covered; low-complexity masking; blastp vs tblastn |
| `defaults/homology-jackhmmer.conf` | E-values; number of iterations |
| `presets/organelle_tblastn.conf` | an example of an organelle screen |

Every value can be changed by a flag, a config file, or `--set`, and each run writes what it used to `run.json`.

## 3. What the core promises, and which check proves it

You do not need to read the code if these promises are what you want and the check passes.

| Promise | Check in `tests/test_core.py` |
|---|---|
| A job that fails is never counted as "no hit"; its genes become `NOT_RUN` | `run_jobs failure reported`, `tally` |
| A misspelled setting stops the run instead of being ignored | `unknown override rejected`, `unknown section key rejected` |
| A setting changed from its default is recorded and warned about | `non_default recorded`, `start_module non_default` |
| Duplicate or empty gene IDs stop the run | `duplicate ids rejected` |
| Two input files for the same species stop the run | `species clash rejected` |
| A gene ID matches whether or not it carries the Stage 1 `__Species_Name` suffix | `split_id ...`, `load ids canonical` |
| Survivors are written with their original headers, and only genes that pass | `passed fasta keeps original header, only passers` |
| When several statuses apply to one gene, the strictest wins (`HIT` beats `SPURIOUS`, and so on) | `classify`, `passes` |
| A stopped run can be resumed without redoing finished chunks | `resume skips done, retries failed` |
| The organelle example turns on the analysis species' own organelle genome | `organelle example` |
| The end-of-run note names only the tools the run used, never invents a version, and can be silenced (the file is still written) | `note names the tool and version`, `missing tool version is stated, not invented`, `quiet prints nothing but still writes citations.txt` |

## 3b. `domain-cdd` (built)

Read `docs/domain-cdd.md` (what counts as a hit, and what to be careful about), then the checks in
`tests/test_domain_cdd.py`, whose header lists each mock gene and the answer worked out by hand.

| Promise | Check |
|---|---|
| A hit needs the E-value and enough of the domain covered; two alignments to one domain add up, overlaps counted once | `two alignments to one model ...`, `min_domain_cov=90`, `evalue=1e-10 ...` |
| Each hit reports its CDD accession and source; a rejected hit says why | `a hit reports CDD accession and source`, `a rejected hit says why ...` |
| `ignore_sources` discounts a source but still lists the hit | `ignore_sources=COG ...`, `ignored hit is still reported ...` |
| A result that never finished is `NOT_RUN`, never "no hit" (at the end of the file, or in the middle) | `result without its completion footer ...`, `an unfinished block in the middle ...` |
| A wrong file format is refused | `wrong column count stops the run ...` |
| Real RPS-BLAST output (header before every query, one footer at the end) is read correctly | the `real output: ...` checks, using `example_data/domain-cdd/real_sample/` |
| A result cut short is never trusted, even when a finished run follows it in the same file | `the first run was cut short ...`, `a footer whose query count does not match ...` |
| An earlier run's results are not overwritten by accident, and unrelated folders are never written into | `same folder, same run, again: refused`, `a folder that already holds other files is refused ...`, `running again into the same output folder is refused ...` |
| `summary.tsv` counts match the hand-worked answers | `summary.tsv ...` |
| The tests notice if the rules are broken | done by hand: six deliberate bugs, each made tests fail |

## 3c. `domain-interpro` (built)

Read `docs/domain-interpro.md`, then the checks in `tests/test_domain_interpro.py`, whose header lists each mock gene and
the answer worked out by hand.

| Promise | Check |
|---|---|
| Domain and family databases count; MobiDB-Lite, Coils and PROSITE patterns never do; AntiFam is flagged, not counted | `default statuses`, `MobiDB-Lite and Coils rows are listed but do not count ...` |
| Decided by analysis, not by InterPro accession: real domain matches with no accession still count | `a domain match with no InterPro accession still counts ...`, `real output: domain matches without an InterPro accession ...` |
| A real domain outweighs an AntiFam flag; 22 MobiDB-Lite rows plus one real domain is a hit | `the AntiFam row on a gene with a real domain ...`, `real output: 22 MobiDB-Lite rows ...` |
| The optional E-value cut-off never touches score-type analyses (Hamap, PROSITE profiles) | `max_evalue applies to E-value analyses ...`, `Hamap's column 9 is a score ...` |
| InterProScan's TSV lists only proteins with matches, so "no matches" needs proof the run finished; without proof it is `NOT_RUN`, never `NO_HIT` | `no evidence the run finished ...`, `--assume-complete ...`, `a status file saying 'incomplete' is not proof` |
| `-dp` (no lookup: local and offline) is always on and cannot be turned off | `the command line always has -dp`, `no setting can switch -dp off`, `run path: every call has -dp` |
| A chunk that fails, or exits 0 without InterProScan's `100% done` line, gives `NOT_RUN` for its genes only | `a chunk that fails ...`, `exit code 0 without InterProScan's '100% done' line is NOT trusted ...` |
| Parallel chunks never share a temporary folder | `run path: every call has its own temp folder` |
| A wrong file (header line, wrong column count) is refused | `wrong column count ...`, `a header line is refused ...` |

## 3d. `homology-blast` (built)

Read `docs/homology-blast.md`, then the checks in `tests/test_homology_blast.py`, whose header lists each mock gene and the
answer worked out by hand.

| Promise | Check |
|---|---|
| Any hit counts, subject to the E-value and (optional) coverage and identity settings; several alignments to one target add up | `default statuses`, `two alignments to one target are merged ...`, `min_qcov=0`, `min_pident=40`, `evalue=1e-12 ...` |
| Hits to the species you are analysing (or to `ignore_taxa`) are listed but do not count; a gene whose only hits are to itself is `IN_PHYLOGENY_ONLY` and passes | `a hit to the analysis species is listed but not counted ...`, `ignore_in_phylogeny=no ...` |
| A hit of unknown species is counted, never exempted | `a hit to a target with no species tag counts ...` |
| Subspecies and strains are the same species by default; `exact` and `genus` levels exist | `ignore_taxa: a subspecies is treated as the species`, `match_level=exact ...`, the `taxon_key` and `matcher` checks in `test_core.py` |
| Organelle search: the analysis species' own organelle genome counts under the preset; no target coverage for `tblastn` | `the organelle preset ...`, `tblastn does not report target coverage ...`, `min_tcov with tblastn is refused ...` |
| BLAST+ refuses `-ungapped` without `-comp_based_stats 0`; the module says so before running | `ungapped without comp_based_stats=0 is refused ...` |
| A cut-short result is never trusted, and a gene is never "no hit" by accident | `result without BLAST's completion footer ...`, `the first run was cut short ...`, `a footer whose query count does not match ...`, `output without BLAST's completion footer is NOT trusted ...` |
| A database of the wrong type, built without `-parse_seqids`, or with untagged IDs is caught | `a nucleotide database with blastp is refused ...`, `a database built without -parse_seqids ...`, `a database whose IDs carry no species tag ...` |
| Several databases: a broken one does not stop the others | `several databases: a broken one is reported ...` |
| `--sequential`: each database searches only what passed the previous one; a failed database stops the chain so nothing passes unchecked | `--sequential: the second database searches only the 3 survivors ...`, `--sequential with a broken database: the chain STOPS ...`, `WITHOUT --sequential every database searches every gene ...` |
| No table the tool writes has an empty cell; missing values are `NA` | `no output table has an empty cell ...` (in every module's tests) |
| Real BLAST+ output is read correctly | the `real blastp output ...` and `real tblastn output ...` checks, using `example_data/homology-blast/real_sample*` |
| The tests notice if the rules are broken | done by hand: eight deliberate bugs, each made tests fail |

## 3e. `homology-jackhmmer` (built)

Read `docs/homology-jackhmmer.md`, then `tests/test_homology_jackhmmer.py`, whose header lists each mock gene and the answer
worked out by hand.

| Promise | Check |
|---|---|
| The whole-sequence E-value decides first; then coverage from the domains whose own E-value also passes (overlaps merged, not added) | `default settings give the hand-worked status for all 14 genes`, `query coverage: ...` |
| jackhmmer writes nothing for a query without hits, so a finished file (`# [ok]`) is the only proof a query was searched | `a result file without '# [ok]' is NOT trusted ...`, `output without HMMER's '# [ok]' is NOT trusted ...` |
| A failed job makes only its own queries `NOT_RUN`, and they stay `NOT_RUN` when the raw files are re-read | `a chunk that fails ...`, `a failed chunk's genes are NOT in the searched-list ...` |
| A hit of unknown species counts; analysis species and `ignore_taxa` are exempt; subspecies and strains match their species | `the target species comes from the ID tag ...`, the `ignore_taxa` checks |
| A BLAST database prefix is refused, and an analysis species inside the database is named in a warning | `a BLAST database prefix is refused ...`, `an analysis species INSIDE the database is named ...` |
| `--sequential` works and stops on a failed database | the `--sequential: ...` checks |
| Real HMMER 3.4 output is read correctly | `real HMMER output: ...` (uses `example_data/homology-jackhmmer/real_sample`) |
| The tests notice if the rules are broken | done by hand: nine deliberate bugs, each made tests fail (two were missed at first and their tests were added) |

## 3f. `make-db` (built)

Builds the species-tagged FASTA that both homology modules use as a database. `tests/test_make_homology_db.py` checks: IDs are
tagged and existing tags kept; `*` removed; duplicate IDs refused with no output left; analysis species (and their subspecies)
left out; no overwrite without `--force`; a species table with no blank cell.

## 3g. `combine` (built)

Read `docs/combine.md`, then `tests/test_combine.py`, whose header lists three mock runs and each gene's expected verdict.

| Promise | Check |
|---|---|
| Verdicts: `REMOVED` (any run removed it), `INCOMPLETE` (not removed, but some run did not finish or never searched it), `PASS` | `verdicts: g1 REMOVED ... g10 INCOMPLETE` |
| A gene an earlier run removed is `NOT_SEARCHED` later, never 0 or blank, and `dropped_by` says which run stopped it | `a gene an earlier run removed is NOT_SEARCHED ...`, `dropped_by names the run ...` |
| A gene never given to a run for no known reason is `INCOMPLETE` with a warning | `g10: never given to cdd ...` |
| Breadth counts (species, genera) come from counted homology hits only, and are `NA` (not 0) when no homology run searched the gene | the `breadth: ...` checks |
| No table has an empty cell; `--na` changes the placeholder but a blank is refused | `no table has an empty cell`, `--na changes the placeholder ...`, `a blank --na is refused` |
| Run files that disagree with each other are refused, not guessed at | `run files that disagree ... are refused` |
| A passed gene whose sequence cannot be found is an error, never silently dropped | `a passed gene whose sequence cannot be found is an error ...` |
| Works on real run folders | the `real runs: ...` checks (use `example_data/combine/real_sample`) |
| The tests notice if the rules are broken | done by hand: eleven deliberate bugs, each made tests fail |

## 3h. `db-report` (built, optional)

Read `docs/db-report.md`, then `tests/test_db_report.py`, whose header draws the small taxonomy the answers come from.

| Promise | Check |
|---|---|
| Species inside your analysis clade are found and flagged | `Arabidopsis thaliana and Brassica rapa (strain FPsc) are INSIDE ...` |
| Tiers are the clade's ancestors, nearest first; empty and thin tiers are flagged | `tier summary: ...`, `tier order: ...`, `a tier below --thin is THIN ...` |
| A name that cannot be placed, or is ambiguous, is `UNMATCHED`, never guessed; common names are not species names | `an abbreviated name ... is UNMATCHED ...`, `a common name ... is NOT accepted ...`, `an ambiguous clade name is refused ...` |
| An accidental `__` inside a random ID is not read as a species | `untagged IDs are counted and the accidental '__' ... is NOT read as a species` |
| Works without the taxonomy, from your own tiers table | `--tiers with --species-list: ...` |
| No blank cells | `no table has an empty cell` |
| Works on the real NCBI taxonomy | done by hand on 2026-09-26 (bryophyte test database and a real plants database); the tests use a small made-up taxonomy in the real file layout |

## 3i. `make-bed` (built, optional, general purpose)

Read `docs/make-bed.md`, then `tests/test_make_cds_bed.py`, whose header explains the two hand-built GFF3 files (mirroring
Phytozome and MAKER, the two real annotation sources this was checked against) the answers come from.

| Promise | Check |
|---|---|
| Primary matching (mRNA ID/Name/gene/protein_id) and the CDS-only fallback for GFFs without a usable mRNA Name both find the right ID | `AT1G01010.1: 2 CDS segments ...`, `fallback match ... still found`, `MAKER-style file ...` |
| A duplicated CDS line (a real annotation-export quirk) is written once and counted once, in both the BED and the tally | `the identical duplicate CDS line is NOT double-counted ...`, `the BED itself also has exactly one row ...` |
| An ID that matches nothing is reported in `not_found.tsv`, never silently dropped, and counts the run as incomplete (exit 2) | `AT9G99999.1 is reported in not_found.tsv ...` |
| A species whose ID list or GFF cannot even be read is warned loudly (not just logged) and also counts as incomplete | `an unknown --id-column is loud (warned) and counted as incomplete ...` |
| Batch pairing by `--species-map` is exact; guessed pairing refuses rather than picks when two files could match | `two GFFs matching the same guessed prefix is refused ...`, `two id files that both single-match the same GFF ... is refused ...` |
| Gene length and gene ID are `NA`, never a fake 0, when nothing resolves them (even through the fallback match, which still finds the real parent-of-parent gene when one exists) | `gene_length is still resolved for a fallback-matched CDS ...`, `orphanRNA has no Parent= at all ... NA, never a fake 0` |
| No blank cells; output filenames are filesystem-safe even when the species label has a space | `no table has an empty cell`, batch-mode filename check |
| Real GFF3 data, network cut off | 2026-09-26: all 23 species of a real Brassicaceae project (about 2.5 GB), both annotation styles, `--species-map` and guessed pairing: 79,631 of 79,631 IDs found, 0 warnings, exact parity with the script this was rewritten from |
| The tests notice if the rules are broken | done by hand: eleven deliberate bugs, each made tests fail |

## 4. Not covered yet

- **Real data at scale.** Every search module was run for real on 200 proteins of one moss and on databases of about 250,000
  sequences (InterProScan on a subset of its analyses). Not yet tried: Gene3D, SUPERFAMILY and NCBIfam data, full-size
  databases and candidate sets, a search that fails part-way on real data, and `combine` on real `domain-interpro` folders.
- **No benchmark of the defaults.** The cut-offs are standard choices (`docs/thresholds.md` says what was checked against the
  tools' own help text), not tuned. Whether they suit a given clade is for the project to judge.
- **The species exemption is one global set, by design.** The modules exempt the species of your input files (plus `ignore_taxa`) for every
  gene; the Casola Lab's earlier filter exempted, per gene, only the species inside the clade of the node where the gene originated.
  The two agree when the database holds no species of your analysis tree. If it holds one (an outgroup, say) that has no input file,
  exclude it from the database or list it in `ignore_taxa` on purpose; `docs/interpreting_hits.md` and `docs/databases.md` explain.
- **Platforms.** Only Linux and WSL (Ubuntu 24.04, Python 3.12) and Python 3.11 on Windows have been used. macOS and the workstation
  have not.
- **`make-bed` checked on four real annotation sources now** (Phytozome, MAKER, NCBI RefSeq, Ensembl): 150/150 and 150/150
  found on real NCBI and Ensembl protein accessions (through the CDS-only fallback in both cases); 134/150 on real NCBI
  locus tags, with all 16 misses confirmed by hand to be non-coding genes with no CDS to find.
- `check` was tested with real BLAST+ and HMMER installs in WSL and with stand-in programs for the rest; it has no automated test.

## 5. Where the code is

| File | Lines (about) | Role |
|---|---|---|
| `bin/hep_common.py` | 950 | shared helpers: settings, FASTA, species, output files, running jobs, the BLAST table reader |
| `bin/domain_interpro.py`, `bin/domain_cdd.py` | 380, 280 | the two domain modules |
| `bin/homology_blast.py`, `bin/homology_jackhmmer.py` | 540, 400 | the two homology modules |
| `bin/make_homology_db.py`, `bin/combine.py`, `bin/db_report.py` | 140, 350, 530 | build a database, merge runs, check a database (optional) |
| `bin/make_cds_bed.py` | 400 | look up CDS coordinates for an ID list in a GFF3 (optional, general purpose) |
| `scripts/homolog_exclusion_pipeline.sh` | 60 | picks which module to run |
| `scripts/check_tools.sh` | 230 | reports which tools and databases are found; changes nothing |
| `setup/get_data.sh`, `setup/get_taxonomy.sh` | 130, 90 | optional downloads (the only steps that use the internet) |

For an independent second opinion on the code, run `/code-review` in Claude Code on this folder.
