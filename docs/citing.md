# Citing

If you use this pipeline, please cite it, and cite the programs and databases your runs used. Each run ends with a note naming the
tools it used and their versions, and writes the same note to `citations.txt` (`combine` lists the tools of all the runs it read).
The full references are here; the note does not repeat them.

> The ten references below were checked against Crossref (the DOI registry) on 2026-09-26: authors, title, journal, year, volume and
> pages or article number all match, and each carries its DOI. The citation of this pipeline itself is completed when a release is
> archived, and the LINGUA framework paper is added once it is available.

## This pipeline

Owoyemi A, Casola Lab (Protein Evolution Lab, Texas A&M University). *homolog-exclusion-pipeline*: Stage 3 of the LINGUA
comparative genomics framework. Version 0.1.0. See `CITATION.cff`. *(A DOI will be added when a release is archived.)*

Stage 1 (`protein-preprocessing-isoform-pipeline`) and Stage 2 (`gene-stratigraphy-pipeline`) of LINGUA produce this pipeline's input
and have their own citations.

## What to cite for what you ran

| You ran | Cite |
|---|---|
| `domain-interpro` | InterProScan (Jones et al. 2014) and InterPro (Blum et al. 2025), plus the member databases whose matches you counted (for example Pfam, SMART, PANTHER) |
| `domain-cdd` | CDD (Wang et al. 2023), CD-Search (Marchler-Bauer & Bryant 2004), BLAST+ (Camacho et al. 2009) |
| `homology-blast` | BLAST+ (Camacho et al. 2009) |
| `homology-jackhmmer` | HMMER (Eddy 2011) and the jackhmmer iteration procedure (Johnson et al. 2010) |
| `db-report` | NCBI Taxonomy (Schoch et al. 2020) |
| your candidate genes | OrthoFinder (Emms & Kelly 2019), if that is how Stage 2 grouped them |

## References

Author lists are given in full for papers with up to four authors and shortened to the first author otherwise; give the full list
from the journal when you cite.

- Blum M, et al. InterPro: the protein sequence classification resource in 2025. *Nucleic Acids Research* 2025; 53(D1): D444–D456.
  doi:10.1093/nar/gkae1082
- Camacho C, et al. BLAST+: architecture and applications. *BMC Bioinformatics* 2009; 10: 421. doi:10.1186/1471-2105-10-421
- Eddy SR. Accelerated profile HMM searches. *PLoS Computational Biology* 2011; 7(10): e1002195. doi:10.1371/journal.pcbi.1002195
- Emms DM, Kelly S. OrthoFinder: phylogenetic orthology inference for comparative genomics. *Genome Biology* 2019; 20: 238.
  doi:10.1186/s13059-019-1832-y
- Johnson LS, Eddy SR, Portugaly E. Hidden Markov model speed heuristic and iterative HMM search procedure. *BMC Bioinformatics*
  2010; 11: 431. doi:10.1186/1471-2105-11-431
- Jones P, et al. InterProScan 5: genome-scale protein function classification. *Bioinformatics* 2014; 30(9): 1236–1240.
  doi:10.1093/bioinformatics/btu031
- Marchler-Bauer A, Bryant SH. CD-Search: protein domain annotations on the fly. *Nucleic Acids Research* 2004; 32: W327–W331.
  doi:10.1093/nar/gkh454
- Schoch CL, et al. NCBI Taxonomy: a comprehensive update on curation, resources and tools. *Database* 2020; baaa062.
  doi:10.1093/database/baaa062
- Shah N, Nute MG, Warnow T, Pop M. Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows.
  *Bioinformatics* 2019; 35(9): 1613–1614. doi:10.1093/bioinformatics/bty833 (Why `max_target_seqs` is set large.)
- Wang J, et al. The conserved domain database in 2023. *Nucleic Acids Research* 2023; 51(D1): D384–D388. doi:10.1093/nar/gkac1096

The pipeline's own version and the versions of the tools are recorded in each run's `run.json`.
