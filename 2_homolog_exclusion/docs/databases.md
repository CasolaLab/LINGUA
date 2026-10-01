# Homology databases: building, tagging and ordering

`homology-blast` and `homology-jackhmmer` search databases that YOU build. The tool ships none and assumes nothing about
their names, number or order. This page is about the technical form; which species go in which database is the project's decision.

## The one thing to make: a species-tagged FASTA

One protein FASTA per database, in which every ID is `<ID>__<Species_Name>` (Stage 0 `--add-species` writes this form).
The species tag is how the tool knows where a hit came from (so it can leave out hits from your analysis species).

```bash
bash scripts/homolog_exclusion_pipeline.sh make-db -i per_species_fasta/ -o db/mosses.faa --exclude-inputs candidates/
```

`make-db` reads a folder of per-species FASTA files (species = file name without the extension; `--species-map` overrides),
adds the tag to IDs that lack one, keeps a tag that is there only if it looks like a species name (`Genus_species...`; a `__`
inside a random-looking ID, such as `Acid_ENSB:xy__FwSs-9Ag`, is not a species, so `make-db` adds its own tag after it), strips `*`,
refuses duplicate IDs, and writes `db/mosses.faa.species.tsv` (species, source file, sequence count, already-tagged count,
retagged count, INCLUDED/EXCLUDED; no empty cells).
It never overwrites without `--force`, and leaves no half-written file if it fails.

**Leave out the analysis species and their close relatives.** `--exclude-inputs DIR` drops the species of the FASTA files in DIR;
`--exclude-species NAME ...` drops named ones. Comparison is at species level, so a subspecies or strain counts as the same
species; `--match-level exact|genus` changes that. Each search run also warns, by name, when an analysis species is found inside
its database.

**Species of your analysis tree that have no input file** (an outgroup, for example) are not known to the modules as analysis species,
so hits to them count. If such a species is in your database, decide on purpose: leave it out (`--exclude-species`), or keep it and
exempt its hits with `ignore_taxa` (a species name, or a single word for a whole genus). `db-report --clade-species` lists the
database species that sit inside your tree's clade.

## Using it

- `homology-jackhmmer --db db/mosses.faa` reads the FASTA directly.
- `homology-blast --db db/mosses.faa` builds a BLAST index from it the first time (`makeblastdb -parse_seqids`, local) and
  reuses it while the FASTA is unchanged. The index goes in the run's output folder, or a shared folder with `--index-dir`.
  It is a disposable cache: delete it and it is rebuilt. A database already made with `makeblastdb -parse_seqids` can be
  given by its prefix instead.

Several databases: repeat `--db`. Default: each is searched against all genes (separate labelled runs). With `--sequential`
each searches only the survivors of the one before, in the order given, and stops if one fails.

**Using several databases of different size: set `dbsize`.** An E-value depends on the size of the database searched, so with one
fixed E-value the same alignment passes more easily in a small database than in a large one (a 70 MB database against a 900 MB
one is a factor of about thirteen in search space). `dbsize` (BLAST `-dbsize`, jackhmmer `-Z`) fixes the size used for the E-value, so
the databases are held to the same standard. It is off by default, because the right value is a choice: pick one number, for
example the size of your largest database (in residues for BLAST, in sequences for jackhmmer), and give it to every database's run
(`--dbsize N`). With a single database, or databases of similar size, you can leave it off. `run.json` records what was used.

## Size and speed (measured, 8 real moss and liverwort proteomes)

| item | size / time |
|---|---|
| tagged FASTA (248,879 sequences) | 106 MB |
| BLAST index | 154 MB, built in 6 s |
| `blastp -subject FASTA` (no index) | same hits, 1.7 to 2.7 times slower, one thread: not used |
| `make-db` on the same proteomes (one excluded) | 1.4 s |
| `blastp` of 200 proteins against 227,344 sequences | 5 min on one thread |

## Designing a database for your clade

This is one way to do it, based on how an experienced evolutionary biologist reasons about the problem. It is guidance, not a
standard: the tool accepts any species-tagged FASTA, and your own knowledge of your clade should override anything here.

**What the database is for.** The candidate genes come from an orthology analysis that only saw the species you gave it. A gene
can look restricted to your clade because the lineages that hold its homolog were never in that analysis, or because the match was
too diverged for the grouping to catch. The database stands in for everything outside your clade that the analysis never saw. A
gene that has a homolog there is not lineage-specific, so it is removed.

### Three properties of a good database

1. **Depth: a ladder of tiers, from just outside your clade to the most distant.** The closer a relative, the easier a homolog
   is to detect, so the nearest outside lineages remove the most genes. But genes can also match very distant lineages (ancient
   genes, horizontal transfer, contamination), so the ladder should reach the far end. A sketch for a plant clade: close outside
   relatives (the same order or family group, outside your clade), the rest of the same major group, other land plants, algae,
   fungi and animals, other eukaryotes (protists), prokaryotes. For another kingdom, build the equivalent ladder: what is just
   outside your clade, then progressively further out.
2. **Breadth: many independent lineages at each tier, not many species of one lineage.** Sampling one clade densely adds size and
   redundancy but little new sequence diversity: a database where a few species hold a large share of all sequences mostly
   re-tests the same protein families. As many species as exist at the nearest tier, and one to a few per family or order at the
   far tiers, is a reasonable starting point. (Any numbers here are proposals that have not been tested; the checks below are how
   you find out whether yours is enough.)
3. **Quality: one protein per gene, reasonably complete proteomes, and no analysis species or close relatives inside.** Use one
   (usually the longest) isoform per gene, as the Stage 0 preprocessing does; leave out fragmentary or badly annotated proteomes
   and ones likely to be contaminated. Leave out every species of your analysis and every relative that could belong in your
   tree: this is what `make-db --exclude-inputs` is for, and each search run warns by name when an analysis species is found
   inside a database.

### The search matters as much as the database

A database only detects what the search can reach. `homology-blast` is fast and suits the nearest tiers; `homology-jackhmmer`
builds a profile and reaches remote homologs that BLAST misses, at a much higher cost. A common arrangement is to put the
databases most likely to hit first and search them one after another with `--sequential`, so every gene removed early is a gene
the slower searches never see.

### A caveat: public databases are uneven

NCBI, Ensembl, UniProt and similar sources are biased toward the most-studied clades (model organisms, crops, medically
important species). Some lineages have hundreds of proteomes and others have none. So "no hit" means "no homolog **in the
sequences available**": it is weaker evidence for a poorly sampled lineage than for a well-sampled one. Work with what exists,
sample as widely as it allows, and say in your methods which tiers were thin.

### Organelle databases

The same reasoning applies on a smaller scale, and the search is `tblastn` (see `docs/homology-blast.md`). The nearest tier is
the organelle genomes of the closest available relatives; then a wider spread of the same kingdom's organelle genomes; then a few
distant ones (for a plastid: algae and cyanobacteria; for a mitochondrion: other eukaryotes and alpha-proteobacteria). Organelle
genomes are small, so a database of a few hundred is cheap. Keep the analysis species' own organelle genome **in** it: a hit to
it is the signal that a "nuclear" gene is really organellar (use `ignore_in_phylogeny=no`).

### Checking that it is good enough

No database is complete. These checks show how good yours is:

- **Look at what you have, by tier:** how many species and how many sequences sit in each tier, which tiers are empty or thin,
  and whether one species holds a large share. (`db-report` will do this; it is optional and planned. A spreadsheet works too.)
- **Positive controls:** genes you already know have homologs outside your clade (for example the analysis species' own genes
  from old, widely shared gene families). Most should hit. If many do not, the database or the search is too weak.
- **Negative controls:** shuffled sequences (the tests use them) should not hit.
- **Saturation:** add species stepwise, or search the database in parts, and look at how many additional candidates each part
  removes. If a tier removes almost nothing that the others did not, it is cheap insurance; if each added species keeps removing
  new candidates, the database is not yet saturated and needs more lineages of that kind. (`combine`'s overlap table will show
  which database removed genes that no other did.)
