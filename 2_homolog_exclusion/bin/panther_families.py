#!/usr/bin/env python3
"""
panther_families.py - build the PANTHER family/species table used by domain-interpro's PANTHER family rule.

Reads the PANTHER sequence classification files (one per genome, e.g. PTHR19.0_arabidopsis from
https://data.pantherdb.org/ftp/sequence_classifications/<version>/PANTHER_Sequence_Classification_files/), whose first column
starts with a five-letter species code (ARATH|...) and whose fourth column is the family:subfamily (PTHR10000:SF1).
Writes family <TAB> species code <TAB> number of genes, one row per pair, sorted.

  python3 panther_families.py PTHR19.0_* > panther_family_species.tsv
Then: domain-interpro --set panther_family_species=panther_family_species.tsv --set panther_exclude_species=ARATH,BRANA
"""
import collections
import os
import re
import sys

FAM = re.compile(r"^(PTHR\d+)")


def table(paths):
    n = collections.Counter()
    for p in paths:
        if not os.path.isfile(p):
            sys.exit("not a file: %s" % p)
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                c = line.rstrip("\r\n").split("\t")
                if len(c) < 4:
                    continue
                m = FAM.match(c[3])
                sp = c[0].split("|", 1)[0].strip().upper()
                if m and sp:
                    n[(m.group(1), sp)] += 1
    return n


def main(argv=None):
    paths = (argv if argv is not None else sys.argv[1:])
    if not paths or paths[0] in ("-h", "--help"):
        sys.stdout.write(__doc__)
        return 0 if paths else 2
    n = table(paths)
    if not n:
        sys.exit("no family rows found: are these PANTHER sequence classification files?")
    w = sys.stdout
    w.write("family\tspecies\tn_genes\n")
    for (f, sp), k in sorted(n.items()):
        w.write("%s\t%s\t%d\n" % (f, sp, k))
    sys.stderr.write("%d families, %d species\n" % (len({f for f, _ in n}), len({s for _, s in n})))
    return 0


if __name__ == "__main__":
    sys.exit(main())
