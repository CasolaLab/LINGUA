# Thresholds: every default, where it comes from, and how far it has been checked

Every number below is a **default**. Change any of them with a flag (`--evalue 1e-5`), `--set key=value`, or a config file; each
run records what it used in `run.json`, and prints a warning naming every setting that differs from its default. The defaults are
standard, common choices for filtering, not the settings of any one project. Your data may need others: a divergent clade may
want a looser E-value; a project that wants fewer false removals may want a stricter one.

"Checked" below means: read against the help text of the installed program (BLAST+ 2.17.0, HMMER 3.4) on 2026-09-26, or, where
stated, taken from the documentation named. Nothing here has been tuned on a benchmark.

## All modules

| Setting | Default | Meaning |
|---|---|---|
| `species_strip_suffixes` | `_final` | removed from the end of an input file's name to get the species (Stage 1 writes `<species>_final.faa`); `--species-map` overrides |

## domain-interpro

| Setting | Default | Why | Checked |
|---|---|---|---|
| `applications` | Pfam, PRINTS, SMART, NCBIfam, SFLD, SUPERFAMILY, Gene3D, Hamap, ProSiteProfiles, CDD, PIRSR, PIRSF, plus MobiDBLite, Coils, ProSitePatterns, AntiFam | the 12 member databases the Casola Lab uses, plus four run for visibility (below). This is a lab selection, not InterProScan's own default set; other labs may choose differently | names confirmed against a real 5.78 installation by `check` |
| `exclude_analyses` | MobiDBLite, Coils, ProSitePatterns | disorder, coiled coils and short patterns that match by chance are not domains; on a real test of 100 moss proteins MobiDB-Lite alone produced a large share of all matches | tested on real output |
| `spurious_analyses` / `antifam_action` | AntiFam / `flag` | AntiFam flags a spurious gene model, not a domain; the gene still passes unless you set `remove` | decided with the user |
| `ignore_domains` | empty | known lineage-specific domains of your clade that should not count | (yours) |
| `max_evalue` | empty (off) | each member database has its own curated cut-offs, which InterProScan applies; an extra E-value applies only to analyses that report one (Hamap and PROSITE profiles report a score, Coils/MobiDB-Lite/PROSITE patterns report `-`) | column meaning read from real output |
| `chunk_size`, `parallel_jobs` | 1000, 4 | speed only | |

The four extra analyses do not change the answer: ignoring them and not running them give the same pass or fail. They are run so
they can be reported.

## domain-cdd

| Setting | Default | Why | Checked |
|---|---|---|---|
| `evalue` | 1e-2 | the NCBI CD-Search default, and the value in the Casola Lab's methods | not verified locally (NCBI's web documentation; rpsblast's own default is 10) |
| `min_domain_cov` | 50 | percent of the domain **model** covered; a partial match to a model is weaker evidence of a domain | a proposal |
| `min_qcov` | 0 | off; the domain, not the query, is what matters here | |
| `ignore_sources` | empty | CDD source databases (for example COG, KOG) whose hits should not count; every hit's accession and source is written to `hits.tsv` | |
| `comp_based_stats` | 1 | rpsblast's own default (`-comp_based_stats D` is equivalent to 1) | checked against rpsblast 2.17.0 |

## homology-blast

| Setting | Default | Why | Checked |
|---|---|---|---|
| `program` | blastp | `blastp` (protein database) or `tblastn` (nucleotide database, for the organelle screen) | |
| `evalue` | 1e-3 | a common filtering cutoff; BLAST's own default (10) is a search default, not a filter | BLAST default confirmed as 10; 1e-3 is a choice |
| `min_qcov` | 50 | percent of the **query** covered by the union of a hit's alignments; a short local match should not count as homology | a proposal |
| `min_tcov` | 0 (off) | target coverage would discard a short candidate that matches one domain of a long protein (`blastp` only) | |
| `min_pident` | empty (off) | identity depends on how diverged the species are | |
| `seg` | yes | mask low-complexity regions in the query; they produce spurious hits. BLAST+'s command-line `blastp` does not do this by default, the BLAST web page does | wording from the BLAST+ manual |
| `max_target_seqs` | 5000 | BLAST+ keeps the **first** N hits it finds, not the best N, so a small value can hide hits | help text says "5 or more is recommended"; the first-N behaviour is documented in the BLAST literature |
| `comp_based_stats` | empty | BLAST's own default (blastp: `D` is equivalent to 2) | checked |
| `max_hsps`, `ungapped`, `dbsize` | empty, no, empty | passed through to BLAST+ when set | |
| `ignore_in_phylogeny` | yes | hits from analysis species and `ignore_taxa` are reported but not counted | (`no` for organelle) |
| `match_level` | species | the first two words: a subspecies or strain counts as the same species | |

Organelle screen (`presets/organelle_tblastn.conf`, an example): `tblastn`, ungapped, `comp_based_stats=0`, `max_hsps=5`,
`min_qcov=0`, `ignore_in_phylogeny=no`. Any hit to an organelle genome counts, including one to the analysis species' own.

## homology-jackhmmer

| Setting | Default | Why | Checked |
|---|---|---|---|
| `evalue` | 1e-3 | a common filtering cutoff; HMMER's own default is 10 | HMMER default confirmed as 10 |
| `incE` | 1e-3 | the inclusion threshold for building the profile between iterations; this is HMMER's own default | |
| `iterations` | 3 | 3 is a common depth for remote-homology searches; HMMER's own default is 5; more iterations raise sensitivity and the risk of profile drift | HMMER default confirmed as 5 |
| `min_qcov`, `min_tcov` | 50, 0 | as for BLAST, but from the domains whose own E-value also passes | a proposal |
| `dbsize` | empty | `-Z`: fix the database size so E-values compare across databases | |
| `chunk_size` | 1 | one query per job, because what a query costs depends on how many hits it finds, not its length | measured: about 7% more CPU, wall-clock gain not measured |

## A note on stricter or looser settings

The Casola Lab's earlier scripts used stricter jackhmmer settings (an E-value of 1e-4 or 1e-5) than these defaults, and RPS-BLAST
cut-offs of 1e-4 or 1e-5 in different script versions. A stricter E-value removes fewer genes; a looser one removes more. The
choice belongs to the project, not to the tool: state the values you used in your methods (`run.json` records them).

## Not yet checked

Whether these defaults are the best ones for any clade has not been tested. The tool was run for real on one moss proteome
against bryophyte proteomes, on one chloroplast genome, and on CDD; no benchmark of false and missed removals exists yet.
