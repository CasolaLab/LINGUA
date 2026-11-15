#!/usr/bin/env python3
"""
filter_stop_events.py  (stage 4, step 6b; method options drop_terminal_stop, stops_focal_frame_only)

Remove stop events that are not premature stops of the focal ORF from the step 6 event list, before steps 6c and 7
(Brassicaceae, 2026-11-14). The focal CDS length is read from the focal row of the alignment: the reference_length
column of the stage 3 stop tables is not the CDS length in every dataset (Brassicaceae: values 1-8), which is why the
tree-config key stops_focal_frame_only could not be used there.
  drop_terminal_stop       a stop at the focal ORF's own stop codon (CDS position L-2) is the stop the focal gene also
                           has, not a disabler (Brassicaceae locus KAF8101894.1: TGA at 169 of 171 used as an event).
                           Only when the focal row ends with a stop codon: where the alignments carry the CDS without
                           it (Brassicaceae Arabis alpina, 0 of 3,501), L-2 is the last sense codon (2026-11-15)
  stops_focal_frame_only   only stops in the focal reading frame ((position - 1) % 3 == 0) are events; a stop in a
                           shifted frame is usually the same triplet read after a frameshift, which is counted as an
                           indel event already (T2T rule of 2026-10-09)
Start-codon events (position 1) are in the focal frame and are kept. The step 6 list is kept as
all_stops_nodes_<F>.unfiltered.tsv (read from there on reruns) and every removed event is written to
stops_dir/removed_stop_events_<F>.tsv with its reason.

Usage: python3 filter_stop_events.py STOPS_DIR ALN_DIR --focal-species Arabidopsis_thaliana
"""
import argparse
import collections
import csv
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lingua_config as cfg  # noqa: E402
from orf_length_node import read_aln  # noqa: E402


def die(m):
    sys.exit(f"ERROR: {m}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stops_dir")
    ap.add_argument("aln_dir")
    ap.add_argument("--focal-species", required=True)
    a = ap.parse_args()
    F = a.focal_species
    term_on, frame_on = cfg.opt("drop_terminal_stop"), cfg.opt("stops_focal_frame_only")
    events_path = os.path.join(a.stops_dir, f"all_stops_nodes_{F}.tsv")
    backup = events_path[:-4] + ".unfiltered.tsv"
    if not os.path.isdir(a.aln_dir):
        die(f"not found: {a.aln_dir}")
    if not os.path.exists(backup):
        if not os.path.exists(events_path):
            die(f"not found: {events_path}")
        shutil.copy2(events_path, backup)
    events = list(csv.DictReader(open(backup), delimiter="\t"))
    fields = list(events[0].keys()) if events else ["geneID", "StopCodonStart_RefPos", "Nodes_stops"]

    cds_len = {}

    def L(g):
        if g not in cds_len:
            gk = g[:-2] if g.endswith(".p") else g  # Phytozome IDs carry ".p" in the protein tables, not in alignments
            p = glob.glob(os.path.join(a.aln_dir, f"{F}_{gk}_*.aln"))
            if len(p) != 1:
                die(f"{len(p)} alignments for {g} in {a.aln_dir}")
            cds = read_aln(p[0])[F].replace("-", "").upper()
            cds_len[g] = (len(cds), cds[-3:] in ("TAA", "TAG", "TGA"))
        return cds_len[g]

    kept, removed, reasons = [], [], collections.Counter()
    for r in events:
        pos = int(r["StopCodonStart_RefPos"])
        why = []
        if term_on and L(r["geneID"])[1] and pos == L(r["geneID"])[0] - 2:
            why.append("focal_terminal_stop")
        if frame_on and (pos - 1) % 3:
            why.append("outside_focal_frame")
        if why:
            removed.append({**r, "reason": ",".join(why)})
            reasons[",".join(why)] += 1
        else:
            kept.append(r)
    with open(events_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(kept)
    with open(os.path.join(a.stops_dir, f"removed_stop_events_{F}.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields + ["reason"], delimiter="\t")
        w.writeheader()
        w.writerows(removed)
    assert len(kept) + len(removed) == len(events)
    print(f"{F}: {len(events)} stop events, {len(removed)} removed {dict(reasons)}, {len(kept)} kept", file=sys.stderr)


if __name__ == "__main__":
    main()
