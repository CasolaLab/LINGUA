# Contributing

Bug reports, questions and improvements are welcome. Please open an issue first for anything larger than a small fix.

## Before you change code

1. Read `REVIEW.md` (what each part promises and which check proves it) and `docs/thresholds.md`.
2. Run the tests; every file must say `ALL PASSED`:

```bash
for t in tests/test_*.py; do python3 "$t" | tail -1; done
bash tests/run_offline.sh          # Linux/WSL: the same, with the network cut off
```

## What a change should keep true

- **Standard library only, Python 3.11+.** No packages to install. The tool ships no programs, databases or data.
- **Fail closed.** A search that did not finish, or a gene that was not searched, is never reported as "no hit".
- **No blank cells.** A missing value is `NA`; a `0` must never stand for "not searched".
- **Real output first.** Before writing or changing a parser, run the real program and read its actual output; keep a small real
  sample as a test. Mock data must copy the real format exactly.
- **Test the tests.** After adding a check, break the code on purpose and confirm the check fails.
- **Dataset-agnostic.** No project-specific species, databases or thresholds in the tool. Those belong in your project's own notes.
- **No network in the modules.** Only `setup/get_*.sh` may download, and only when you run them.

## Reporting a bug

Please include the command, the tool versions (`run.json` records them), and the smallest input that shows the problem. Do not
attach unpublished sequences; a few made-up or public sequences that reproduce it are enough.
