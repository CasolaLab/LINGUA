# Stage 2 toy run

The hand-built example of `tests/test_run_stage2.py` (same seed, same sequences, same expected fate of each candidate),
written out as files for the tutorial page [`docs/tutorial/03_stage2_homolog_exclusion.md`](../../../docs/tutorial/03_stage2_homolog_exclusion.md).
`fake_interproscan.sh` stands in for InterProScan and never finds a domain. The BLAST databases are built by the
tutorial with `makeblastdb`, and git ignores the files it writes in `db/`.
