# hog-rule

Runs after `combine`. A candidate can pass every search and still sit in the same OrthoFinder hierarchical orthogroup (HOG) as a
gene that a search removed. OrthoFinder grouped them because the two proteins are similar, so the candidate may share the very
region that gave the removed gene its homolog or domain. `hog-rule` tests that directly.

## What it does

1. Reads `gene_matrix.tsv` and `run.json` from a `combine` folder. Evidence members are genes that an **evidence run** removed
   (`evidence_runs`, by label; leave out screens such as TE or organelle that are not homology to another lineage).
2. Each evidence member's **evidence region** is the union of the query intervals of its counted hits in those runs (column
   `q_intervals` of `hits.tsv`, written by every module since 0.2.0).
3. Every `PASS` gene that shares a HOG (from OrthoFinder's `Phylogenetic_Hierarchical_Orthogroups/N0.tsv`) with an evidence
   member is aligned to those members with BLASTp (`-seg yes`, `evalue`, a fixed `dbsize` so E-values do not depend on how
   many members are searched).
4. For each member, the candidate's alignments are merged on the member's coordinates and intersected with the member's region.
   The candidate is `REMOVED_HOG` when the best overlap is at least `min_overlap_aa` amino acids.

| Verdict | Meaning |
|---|---|
| `PASS` | not in a HOG with an evidence member, no alignment to one, or overlap below `min_overlap_aa` |
| `REMOVED_HOG` | overlap of at least `min_overlap_aa` with an evidence member's region |
| `REGION_UNKNOWN` | aligned to a member whose region has no coordinates, and overlapped no known region enough; does not pass unless `unknown_region=pass`; the run exits with status 2 |

Only genes that `combine` passed are judged; the others keep their `combine` verdict.

## Run it

```bash
bash scripts/homolog_exclusion_pipeline.sh hog-rule -i candidates/ -o out/hog --combine out/combine \
    --hogs OrthoFinder/Results/Phylogenetic_Hierarchical_Orthogroups/N0.tsv --set evidence_runs=blast,jackhmmer,cdd,interpro
```

`-i` is the same candidate FASTA folder the search modules read (the members' sequences come from there). `--parse-only FILE`
reads an earlier within-HOG BLASTp output (outfmt `7 qseqid sseqid evalue bitscore qstart qend sstart send`) instead of running
BLAST+.

## Output

| File | Content |
|---|---|
| `hog_rule.tsv` | one row per `PASS` gene of combine: HOG, verdict, best member, overlap in aa, that member's region |
| `passed_ids.tsv`, `passed/<species>.faa` | the genes that still pass |
| `summary.tsv` | per species: combine passes, `PASS`, `REMOVED_HOG`, `REGION_UNKNOWN` |
| `run.json` | settings, inputs, evidence runs, counts, BLAST+ version |

## Settings

| Setting | Default | Why |
|---|---|---|
| `evidence_runs` | blank = every run | runs whose hits count as homology evidence |
| `min_overlap_aa` | 20 | the Casola Lab's Brassicaceae choice; 1 = any overlap. There is no natural value: on that dataset 1, 10, 20 and 30 aa kept 33,207, 33,614, 34,020 and 34,782 of 44,306 genes |
| `evalue`, `dbsize`, `seg`, `max_target_seqs` | 1e-3, 1e8, yes, 5000 | within-HOG BLASTp |
| `unknown_region` | fail | a pass needs evidence of absence |

## Checked

On the Brassicaceae data, the decision function reproduced the project's own calls for all 16,279 genes in HOGs with an evidence
member (10,286 removed, best overlap identical for every gene), and a run with BLAST+ 2.17.0 on four of those genes reproduced
their verdicts and overlaps. The tests (`tests/test_stage2_rules.py`) cover the decision, the HOG table, the BLAST parser and a
full run with `--parse-only`.
