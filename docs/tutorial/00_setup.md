# 00 Setup

The worked examples need a terminal, conda and about 1 GB of disk for the environment. They run on Linux (including a
shared workstation or an HPC login node), on macOS and on Windows through WSL.

## 1. Get the code

```bash
git clone https://github.com/CasolaLab/LINGUA.git
cd LINGUA
```

Every command in the tutorial is run from this top folder.

## 2. Install conda, if you do not have it

[Miniforge](https://github.com/conda-forge/miniforge) is the recommended installer, because it uses the conda-forge
channel by default. On a shared machine, check first whether conda is already installed (`conda --version`).

## 3. Create the tutorial environment

```bash
conda env create -f docs/tutorial/environment.yml
conda activate lingua-tutorial
```

This one environment covers every worked example. It holds Python 3.12 with pandas and Biopython, BLAST+, HMMER,
GNU parallel, a recent bash and GNU coreutils. Check it with

```bash
python3 -c "import pandas, Bio; print('python ok')"
blastp -version | head -1
jackhmmer -h | sed -n 2p
parallel --version | head -1
bash --version | head -1
```

Each line should print a version, and bash must be 4.3 or later. On a Mac, `/bin/bash` is version 3.2, so always run the
scripts as `bash script.sh` with the environment active, which picks up the newer bash from the environment.

## What the examples do not need

Real runs need four more programs, none of which the worked examples call, because the examples start from their
precomputed output.

| Program | Used by | Install notes |
|---|---|---|
| OrthoFinder | stage 1, to build the orthogroups | `1_gene_stratigraphy/README.md` |
| Progressive Cactus (with HAL tools) | before stage 3, to build the whole-genome alignment | the Cactus documentation; run it on a cluster or a large workstation |
| mafTools (`mafExtractor`) | stage 3 step 1 | the mafTools documentation (ComparativeGenomicsToolkit on GitHub) |
| samtools | stage 3, `build_cds_from_genome.py` | conda (`samtools`, bioconda) |
| InterProScan and the CDD database | stage 2 domain screens | `2_homolog_exclusion/docs/setup.md` |

## Folder for the runs

All examples write into `tutorial_runs/` in the repository folder. Git ignores it, and deleting it resets the tutorial.

Next, [01 Proteins](01_stage0_proteins.md).
