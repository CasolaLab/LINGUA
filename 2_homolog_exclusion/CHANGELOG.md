# Changelog

## 0.3.1 (2026-11-10)

- `paralog-rule`: `max_target_seqs` default 50 → 1000. With 50, genes in large lineage-specific families filled every slot with
  other candidates and an old paralog could fall outside the list; on the Brassicaceae data 1,000 targets remove 52 more genes
  (1,006 instead of 954; 33,014 lineage-specific genes) and no gene reaches the new cap. Preset, protocol and docs updated.

## 0.3.0 (2026-11-10)

- New module `paralog-rule` (`docs/paralog-rule.md`): after `hog-rule`, removes passing genes with a similar established gene
  (not a candidate) in their own proteome, namely one BLASTp alignment (SEG off, E ≤ 1e-5) with at least 50% identity over at
  least 50% of the gene's length. It catches pairs OrthoFinder left unclustered, such as *MT1C*/*MT1A* in *Arabidopsis thaliana*.
- `presets/brassicaceae_2026.conf`: `[paralog-rule]` section; `STAGE2_PROTOCOL.md`: the rule, its command and the study's counts
  (34,020 → 33,066).
- Tests for the module in `tests/test_stage2_rules.py`.

## 0.2.1 (2026-11-09)

Fixes from a review of 0.2.0.

- `hog-rule`: a candidate missing from the within-HOG BLASTp output (possible with `--parse-only`) is now `NOT_RUN` and does not
  pass; before, it passed. `evidence_runs` accepts shell wildcards, so the preset's `blast*` matches the per-database runs that
  several `--db` create (before, the protocol's command stopped with "evidence_runs not found"). `min_overlap_aa` must be at least 1.
- `domain-interpro`: a PANTHER hit whose family is missing from `panther_family_species` is now counted, with a warning, instead of
  being silently excluded.
- `STAGE2_PROTOCOL.md`: states that the study's ten targets per database and genera counted across databases cannot both be
  reproduced, and what each option gives.

## 0.2.0 (2026-11-09)

Added for the Casola Lab's stage-2 re-screen of Brassicaceae candidates (`STAGE2_PROTOCOL.md`). Every new setting is off by
default, so runs made with the defaults give the same answers as 0.1.0.

- `homology-blast`, `homology-jackhmmer`: `filter_evalue` (count hits at a stricter E-value than the search), `min_genera` (a gene
  counts as having homologs only when its counted hits come from at least N genera).
- `homology-jackhmmer`: `max_bias_frac` (a hit whose full-sequence composition bias is at least this fraction of its bit score
  does not count).
- `homology-blast`: `qcov_mode=single_hsp` (coverage from the longest single alignment, the TE rule), `word_size`, `matrix`.
- `domain-cdd`: `ignore_domains` (CDD accessions that do not count).
- `domain-interpro`: `require_evalue`, `max_evalue_exempt`, and the PANTHER family rule (`panther_family_species`,
  `panther_exclude_species`); new helper `panther-families` builds the family table from PANTHER classification files.
- Every module writes the merged query intervals of each hit (`q_intervals`, last column of `hits.tsv`).
- New module `hog-rule` (`docs/hog-rule.md`): after `combine`, removes passing genes whose alignment to a removed member of
  their HOG overlaps that member's hit region by at least `min_overlap_aa`.
- Presets `presets/brassicaceae_2026.conf` and `presets/brassicaceae_2026_te.conf`.
- Tests: `tests/test_stage2_rules.py`.
