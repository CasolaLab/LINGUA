# Setup

You install everything yourself. The repository ships no programs, no databases and no data; it is Python 3.11+ (standard library
only) and shell scripts. Read this once, then use `bash scripts/homolog_exclusion_pipeline.sh check` to confirm the result.

## What you need

| To run | You need |
|---|---|
| any module | Python 3.11 or newer and Linux (or Windows with WSL2). macOS has not been tested. |
| `domain-interpro` | InterProScan 5.x with its data, and Java 11 or newer |
| `domain-cdd` | BLAST+ (`rpsblast`) and the CDD database prebuilt for rpsblast |
| `homology-blast` | BLAST+ (`blastp`, `tblastn`, `makeblastdb`, `blastdbcmd`) and a database you build (`docs/databases.md`) |
| `homology-jackhmmer` | HMMER 3 (`jackhmmer`) and a database you build |
| `combine`, `make-db` | nothing beyond Python |
| `db-report` (optional) | the NCBI taxonomy dump, or your own tiers table |

You need only what the modules you will run need. The modules are independent.

## 1. The programs (conda)

```bash
conda env create -f setup/environment.yml
conda activate hep
```

This installs Python 3.12, BLAST+, HMMER and OpenJDK from **conda-forge and bioconda**. Use Miniforge (it uses conda-forge only).
Do not use Anaconda's `defaults` channel: its terms restrict use at larger organisations. Tested with BLAST+ 2.17.0, HMMER 3.4,
Python 3.12 (the modules also ran on Python 3.11) and Java 25. Python 3.7 is not supported.

## 2. The data

`setup/data_sources.txt` lists every source, URL, size and checksum note. `setup/get_data.sh` can download them; doing it by hand
is equally fine.

```bash
bash setup/get_data.sh cdd --dry-run            # show the size first, download nothing
bash setup/get_data.sh cdd data/cdd             # 1.75 GB download, 3.8 GB unpacked; use --db data/cdd/Cdd
bash setup/get_data.sh interproscan data        # 7 GB download, about 36.5 GB unpacked
```

Downloading is the only step in the tool that needs the internet. The searches never do: `domain-interpro` always runs
InterProScan with `-dp` (no lookup at EBI), and no module sends sequences anywhere.

**Disk space.** Plan for the unpacked size plus the download while unpacking. InterProScan is the big one (36.5 GB for everything
in 5.78; the largest data folders are FunFam 14.9 GB, Gene3D 5.3 GB and PANTHER 4.9 GB). Each analysis needs only its own data
folder, so you can install a subset and list only those analyses in `applications`. On Windows, WSL2 keeps its disk as one
virtual file on `C:`; the free space on `C:` is the real limit, and the file does not shrink when you delete things.

**Slow InterProScan download.** EBI can throttle a single connection heavily (0.3 MB/s was measured; 8 parallel connections gave
about 3 MB/s). If yours is slow, use a multi-connection downloader on the same URL and check the `.md5` next to it.

## 3. Your homology databases

`homology-blast` and `homology-jackhmmer` search databases you build: one species-tagged protein FASTA each. See
`docs/databases.md` for how (`make-db`) and for how to design a database that suits your clade. Organelle databases are nucleotide
FASTA files searched with `tblastn` (`presets/organelle_tblastn.conf`).

## 4. Check it

```bash
bash scripts/homolog_exclusion_pipeline.sh check --cdd-db data/cdd/Cdd --hmmer-db db/mosses.faa \
    --iprscan-bin data/interproscan-5.78-109.0/interproscan.sh
```

`check` is read-only. It reports each tool and its version, InterProScan's analysis names (and that each analysis's data folder
exists), the CDD and BLAST databases, and the jackhmmer database. It exits with 1 if anything you asked for is missing, and says
what to do about it. It downloads and installs nothing.

## 5. Run the tests

```bash
for t in tests/test_*.py; do python3 "$t" | tail -1; done      # every file should say ALL PASSED
bash tests/run_offline.sh                                       # the same, with the network cut off (Linux/WSL)
```

The tests need no InterProScan, BLAST+ or HMMER: they use example output in each tool's exact format and stand-in programs. On
Windows (without WSL) the checks that need stand-in programs are skipped and say so.

## 6. First run

Each module has a folder in `example_data/` with a small input and matching results in the tool's real format. For example, to see
the whole flow without any search program:

```bash
python3 bin/homology_blast.py -i example_data/homology-blast/real_sample/input -o /tmp/try \
    --parse-only example_data/homology-blast/real_sample/blast
```

Then read `docs/workflows.md`.

## What has and has not been tested

Tested on Ubuntu 24.04 under WSL2 with real InterProScan 5.78-109.0 (a subset of its analyses), CDD 3.21, BLAST+ 2.17.0 and HMMER
3.4, with the network cut off. **Not yet tested:** Gene3D, SUPERFAMILY and NCBIfam data (the InterProScan subset did not include
them), full-size databases, macOS, Python 3.11 with the run path on Linux, a search that fails part-way on real data.
