# Interpreting the results

What the tool concludes, and what it cannot conclude. Stage 2 asks one question of each candidate gene: **is there evidence that
this gene is not really specific to your lineage?** A gene that has evidence against it is removed; a gene without such evidence
passes. Passing is the absence of evidence in the searches you ran, never proof.

## The statuses

Every module gives each gene one status:

| Status | Meaning | Passes? |
|---|---|---|
| `HIT` | counted evidence against the gene: a domain (`domain-*`) or a homolog (`homology-*`) | no |
| `NO_HIT` | the search ran and found nothing that counts | yes |
| `EXCLUDED_ONLY` | only signals that are not domains (disorder, coiled coils, short patterns), or domains you chose to ignore | yes |
| `SPURIOUS` | flagged by AntiFam as a probable spurious gene model; not a domain | yes, unless `antifam_action=remove` |
| `IN_PHYLOGENY_ONLY` | its only hits are to your own analysis species (or `ignore_taxa`), which are reported but not counted | yes |
| `NOT_RUN` | the search failed or did not finish for this gene | no |

`NOT_RUN` is never `NO_HIT`. A gene without a finished search is not claimed to have no homolog. `combine` adds `NOT_SEARCHED` (a
later run never received the gene because an earlier run removed it) and the verdicts `PASS`, `REMOVED` and `INCOMPLETE`; see
`docs/combine.md`.

## Exit status

Every search module, and `combine`, ends with one of three exit statuses, so a script can tell them apart:

| Exit status | Meaning |
|---|---|
| 0 | finished, and every gene was searched |
| 1 | could not be done (an error, such as a missing database or a wrong file); nothing usable was written |
| 2 | finished and the results are written, but some genes are `NOT_RUN` (for `combine`: some genes are `INCOMPLETE`); those genes do not pass |

With `set -e` or a workflow manager, status 2 stops the script instead of passing partial results to the next step. The results of the
finished genes are still in the output folder; re-run the failed part (`--resume` skips the batches that finished) before you use the
pass list.

## Reading a `HIT`

Open `hits.tsv`: one row per candidate and target, with the target, its species, the E-value, bit score, coverage, whether it counted
and why not. Before trusting a removal that matters, look at the hit itself:

- **Domain hit (`domain-*`):** the accession and source are given. A domain match is evidence of an old protein family. Check that it
  is a domain and not a repeat or a low-complexity region; MobiDB-Lite, Coils and PROSITE patterns are already excluded.
- **Homolog hit (`homology-*`):** the target's species tells you how far the homolog lies. A hit in a close relative outside your
  clade is different from a hit in a bacterium. `combine` counts the distinct species and genera a gene hits (`n_hit_species`,
  `n_hit_genera`): one weak hit and forty strong hits are different biology.
- **Coverage:** a hit that covers a small part of the query can be a shared domain, not shared ancestry of the gene. The default
  requires half of the query (`min_qcov=50`); change it if you want to catch partial homologs.

## Why a gene may pass when it should not (false passes)

- **The database is thin where the homolog lives.** Public sequence databases are biased toward well-studied clades; "no hit" means
  "no homolog in the sequences available". See "Designing a database for your clade" in `docs/databases.md`.
- **The homolog is too diverged for the search.** BLAST misses remote homologs that jackhmmer finds; jackhmmer misses some that
  a structure-based method would. A very young gene that has diverged beyond recognition looks the same as a gene without a
  homolog.
- **Gene models.** A candidate that is a fragment or a mis-annotated gene may not match even if the real gene would.
- **A relative is missing from the database.** The tool counts a hit to a species only if that species is in the database.

## Why a gene may be removed when it should not (false removals)

- **Your database contains your own clade or a close relative** and the species check did not exempt it. The species check needs
  species tags (`ID__Species_Name`) and full analysis-species names (`--species-map`); each run warns if an analysis species is
  found inside its database. A species of your analysis tree that has no input file (an outgroup) is not an analysis species to the
  modules: leave it out of the database or list it in `ignore_taxa`.
- **Contamination or a horizontally transferred gene.** A hit in bacteria can mean the candidate is a contaminant, not that
  it is old; look at the hit.
- **A weak or partial hit** that passes the E-value. Tighten `evalue`, or raise `min_qcov`.
- **A domain shared by chance**, or a lineage-specific domain of your clade. Use `ignore_domains` for the latter.

## What the tool does not do

It does not decide that a gene is a de novo gene. Passing Stage 2 means that no evidence of a homolog or a known domain was found
by these searches in these databases; other analyses are needed for the rest. Report the databases, the tiers they cover, the tools' versions and the settings
you used (`run.json` records them) so that a reader can judge how much a pass is worth.

## Checks you can run on your own data

- Positive controls: genes you know are old should mostly be `HIT`. If many are not, the database or the search is too weak.
- Negative controls: shuffled sequences should not hit (the tests use them).
- Compare two searches (`combine` on `homology-blast` and `homology-jackhmmer`): genes that only the more sensitive search removed
  show what the fast search would have missed.
