# paralog-rule

Runs after `hog-rule` (or after `combine`). Every screen before it looks for homologs outside the clade, and `hog-rule` looks
inside the candidate's OrthoFinder HOG. Neither sees a similar copy among the established genes of the candidate's own genome
when OrthoFinder left the pair unclustered, which happens with short or low-complexity proteins. Such a candidate is a diverged
duplicate of an established gene rather than a gene without homologs. `paralog-rule` searches for these copies directly.

## What it does

1. Reads `passed_ids.tsv` from the previous step (`--passed`, a `hog-rule` or `combine` folder). Only these genes are judged.
2. Aligns each of them with BLASTp (`-seg no`, so cysteine-rich or other low-complexity proteins are not masked) to the full
   proteome of its species (`--proteomes`, one FASTA per species, named like the `-i` files).
3. Ignores hits to the gene itself and to every other candidate (any gene in the `-i` FASTA files). A copy that is itself a
   candidate makes a lineage-specific family, not a duplicate of an established gene.
4. Removes the gene (`REMOVED_PARALOG`) when one alignment to an established gene has at least `min_pident` % identity and spans
   at least `min_qcov` % of the gene's length, at E ≤ `evalue`.

| Verdict | Meaning |
|---|---|
| `PASS` | no alignment to an established gene meets both cuts |
| `REMOVED_PARALOG` | at least one does |
| `NOT_RUN` | the BLASTp output has no finished record of the gene (possible with `--parse-only`); does not pass; exit status 2 |

## Run it

```bash
bash scripts/homolog_exclusion_pipeline.sh paralog-rule -i candidates/ -o out/paralog --passed out/hog \
    --proteomes proteomes/ --threads 16
```

`--parse-only DIR` reads earlier BLASTp outputs (`DIR/<species>.tsv`, outfmt `7 qseqid sseqid pident evalue bitscore qstart qend
qlen`) instead of running BLAST+.

## Output

| File | Content |
|---|---|
| `paralog_rule.tsv` | one row per judged gene: verdict, the reported established paralog, its identity, the share of the gene's length the alignment spans, E-value |
| `passed_ids.tsv`, `passed/<species>.faa` | the genes that still pass |
| `summary.tsv` | per species: judged, `PASS`, `REMOVED_PARALOG`, `NOT_RUN` |
| `run.json` | settings, inputs, proteomes, counts, BLAST+ version |

The reported paralog is the qualifying alignment with the highest identity, or, when none qualifies, the alignment with the
highest identity.

## Settings

| Setting | Default | Why |
|---|---|---|
| `min_pident`, `min_qcov` | 50, 50 | the Casola Lab's Brassicaceae choice (2026-11); both must hold for one alignment |
| `evalue` | 1e-5 | the E-value at which the study counts BLASTp hits |
| `seg` | no | the study searched without low-complexity masking, so that short cysteine-rich proteins such as *MT1C* are searched over their full length |
| `max_target_seqs` | 50 | targets per gene |

## Checked

On the Brassicaceae data, a full run with BLAST+ 2.17.0 (the 34,020 genes that passed `hog-rule`, the 192,193 candidates as
`-i`, the 23 longest-isoform proteomes) removed 954 genes and passed 33,066, the same 954 genes that the project's own analysis
of a separate BLASTp run had found. The *MT1C* (AT1G07610) control finds *MT1A* (AT1G07600) at 95.6% identity over the full
length; *MT1A* is itself a candidate, so *MT1C* passes. The tests (`tests/test_stage2_rules.py`) cover the hit filter, the
decision and a full run with `--parse-only`.
