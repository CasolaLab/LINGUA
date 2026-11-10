"""
step6b_orf_bearing_nodes_filter.py  (LINGUA stage 4, new step 6b; his rule of 2026-10-29)

A reconstructed ancestral node whose sequence carries the ORF (tBLASTn AncSeq table, per-node 0/1: an ORF of at least
80% of the CDS) cannot be a disabler node: a frameshift or stop reconstructed at that node does not disable that ORF.
The step 6 node tables (all_indels_nodes_<species>.tsv, all_stops_nodes_<species>.tsv) are rewritten with every
ORF-bearing node removed from each event's node list; an event left without nodes is dropped. Step 7 and the step 7d
evidence then both read the filtered tables, so the chosen enablers, the ORF-clade check and the any-valid-disabler
support all see the same events. The per-species detail tables (living species and ancestor rows with lengths and
positions) are not changed: the reconstructed mutation remains a fact of the ancestral sequence, it just is not a
disabler there.

Two readings (his choice 2026-10-29: the strict one, --drop-whole-event):
  per node (default)      the ORF-bearing nodes are removed from the event's node list; the event stays at its other nodes
  per event (--drop-whole-event)   a mutation that coexists with the ORF at ANY node does not disable that ORF, so the
                          whole event is dropped (this keeps step 7's old exclusion of events at the ORF node and only
                          removes the contamination)

Usage:
  python3 step6b_orf_bearing_nodes_filter.py <ancseq_tsv> <all_indels_nodes_tsv> <all_stops_nodes_tsv> <out_dir> [--log LOG_TSV] [--drop-whole-event]
      writes <out_dir>/<same basenames>, and a log with one row per modified event
  python3 step6b_orf_bearing_nodes_filter.py --self-test
Inputs may be gzipped (.gz); outputs are plain TSV. geneIDs are matched exactly, then with a trailing ".p" stripped.
"""
import argparse
import csv
import gzip
import os
import sys
import tempfile

csv.field_size_limit(10 ** 9)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def gk(g):
    return g[:-2] if g.endswith(".p") else g


def read_orf_nodes(path):
    """geneID -> set of node labels whose ancestral sequence carries the ORF (column value '1')."""
    out = {}
    with open_text(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        node_cols = [c for c in rd.fieldnames if c.startswith("N") and c[1:].isdigit()]
        if not node_cols or "Node_AncSeq_ORF" not in rd.fieldnames:
            sys.exit(f"ERROR: {path} has no per-node ORF columns")
        for r in rd:
            g = r["geneID"].strip()
            nodes = {c for c in node_cols if r[c].strip() == "1"}
            out[g] = nodes
            out.setdefault(gk(g), nodes)      # ".p" may be present in one table and absent in the other
    return out


def filter_table(path, nodes_col, orf_nodes, out_path, log_rows, kind, whole_event=False):
    n_in = n_out = n_touched = n_dropped = n_nodes_removed = 0
    with open_text(path) as fh, open(out_path, "w", newline="") as oh:
        rd = csv.DictReader(fh, delimiter="\t")
        if nodes_col not in rd.fieldnames:
            sys.exit(f"ERROR: {path} lacks column {nodes_col}")
        wr = csv.DictWriter(oh, delimiter="\t", fieldnames=rd.fieldnames, lineterminator="\n")
        wr.writeheader()
        pos_col = [c for c in rd.fieldnames if c.endswith("_RefPos")][0]
        for r in rd:
            n_in += 1
            g = r["geneID"].strip()
            orf = orf_nodes.get(g)
            if orf is None:
                orf = orf_nodes.get(gk(g), set())
            nodes = [n for n in r[nodes_col].split("|") if n not in ("", "NA")]
            removed = [n for n in nodes if n in orf]
            kept = [] if (whole_event and removed) else [n for n in nodes if n not in orf]
            if removed:
                n_touched += 1
                n_nodes_removed += len(removed)
                log_rows.append([kind, g, r[pos_col], "|".join(nodes), "|".join(kept) or "-", "|".join(removed)])
            if not kept:
                n_dropped += 1
                continue
            r[nodes_col] = "|".join(kept)
            wr.writerow(r)
            n_out += 1
    assert n_out == n_in - n_dropped
    return n_in, n_out, n_touched, n_dropped, n_nodes_removed


def run(ancseq, indels, stops, out_dir, log_path, whole_event=False):
    for p in (ancseq, indels, stops):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: missing {p}")
    os.makedirs(out_dir, exist_ok=True)
    orf_nodes = read_orf_nodes(ancseq)
    log_rows = []
    stats = {}
    for kind, path, col in (("indel", indels, "Nodes_indels"), ("stop", stops, "Nodes_stops")):
        base = os.path.basename(path)[:-3] if path.endswith(".gz") else os.path.basename(path)
        out_path = os.path.join(out_dir, base)
        if os.path.abspath(out_path) == os.path.abspath(path):
            sys.exit(f"ERROR: output would overwrite input {path}")
        stats[kind] = filter_table(path, col, orf_nodes, out_path, log_rows, kind, whole_event)
    if log_path:
        with open(log_path, "w") as lh:
            lh.write("kind\tgeneID\tposition\tnodes_before\tnodes_after\tnodes_removed_orf_bearing\n")
            for row in log_rows:
                lh.write("\t".join(row) + "\n")
    for kind, (n_in, n_out, n_touched, n_dropped, n_rm) in stats.items():
        print(f"  {kind}: events in {n_in:,}, out {n_out:,}; {n_touched:,} events had {n_rm:,} ORF-bearing node(s); "
              f"{n_dropped:,} events dropped ({'whole event' if whole_event else 'no node left'})", file=sys.stderr)
    return stats


def self_test():
    d = tempfile.mkdtemp(prefix="step6b_")
    anc = os.path.join(d, "x_cdngs_tblastn_3_AncSeq.tsv")
    with open(anc, "w") as fh:
        fh.write("geneID\tNode_AncSeq_ORF\tN21\tN20\tN19\tN18\tN17\tN16\tN15\tN14\tN13\tN12\tN11\tN10\tN9\tN8\tN7\tN6\tN5\tN4\tN3\tN2\tN1\tN0\n")
        # G1: ORF at N3 and every younger lineage node N8..N20 (the 441-PDNG pattern)
        fh.write("G1\tN3\t0\t1\t1\t1\t0\t0\t1\t0\t0\t0\t1\t0\t0\t1\t0\t0\t0\t0\t1\t0\t0\t0\n")
        # G2.p: ORF at N9 only
        fh.write("G2.p\tN9\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\n")
        # G3: no ORF anywhere
        fh.write("G3\tX-specific\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n")
    ind = os.path.join(d, "all_indels_nodes_X.tsv")
    with open(ind, "w") as fh:
        fh.write("geneID\tIndelStart_RefPos\tNodes_indels\n")
        fh.write("G1\t417\tN3|N8|N11|N15|N18|N19|N20\n")       # every node ORF-bearing -> dropped
        fh.write("G1\t366\tN2|N5\n")                            # untouched
        fh.write("G2\t104\tN2|N3|N4|N6|N9\n")                   # N9 removed (matched through G2.p)
        fh.write("G3\t5\tN9\n")                                 # untouched
    stp = os.path.join(d, "all_stops_nodes_X.tsv")
    with open(stp, "w") as fh:
        fh.write("geneID\tStopCodonStart_RefPos\tNodes_stops\n")
        fh.write("G1\t419\tN3|N8|N11|N15|N18|N19|N20|N21\n")   # N21 not ORF-bearing -> kept as N21
        fh.write("G2\t111\tN3|N9\n")                            # -> N3
    out = os.path.join(d, "out")
    log = os.path.join(d, "log.tsv")
    stats = run(anc, ind, stp, out, log)
    got_i = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "all_indels_nodes_X.tsv"))][1:]
    got_s = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, "all_stops_nodes_X.tsv"))][1:]
    assert got_i == [["G1", "366", "N2|N5"], ["G2", "104", "N2|N3|N4|N6"], ["G3", "5", "N9"]], got_i
    assert got_s == [["G1", "419", "N21"], ["G2", "111", "N3"]], got_s
    assert stats["indel"] == (4, 3, 2, 1, 8) and stats["stop"] == (2, 2, 2, 0, 8), stats
    logrows = [l.rstrip("\n").split("\t") for l in open(log)][1:]
    assert len(logrows) == 4 and logrows[0] == ["indel", "G1", "417", "N3|N8|N11|N15|N18|N19|N20", "-", "N3|N8|N11|N15|N18|N19|N20"], logrows
    out2 = os.path.join(d, "out_strict")
    stats = run(anc, ind, stp, out2, None, whole_event=True)
    got_i = [l.rstrip("\n").split("\t") for l in open(os.path.join(out2, "all_indels_nodes_X.tsv"))][1:]
    got_s = [l.rstrip("\n").split("\t") for l in open(os.path.join(out2, "all_stops_nodes_X.tsv"))][1:]
    assert got_i == [["G1", "366", "N2|N5"], ["G3", "5", "N9"]], got_i          # G2's event touched N9 -> whole event dropped
    assert got_s == [], got_s
    assert stats["indel"] == (4, 2, 2, 2, 8) and stats["stop"] == (2, 0, 2, 2, 8), stats
    print("  self-test passed (per node and per event)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("ancseq", nargs="?")
    ap.add_argument("indels", nargs="?")
    ap.add_argument("stops", nargs="?")
    ap.add_argument("out_dir", nargs="?")
    ap.add_argument("--log", default=None)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--drop-whole-event", action="store_true", help="strict reading: drop the whole event if any node is ORF-bearing")
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return
    if not all((a.ancseq, a.indels, a.stops, a.out_dir)):
        ap.error("ancseq, indels, stops and out_dir are required (or --self-test)")
    run(a.ancseq, a.indels, a.stops, a.out_dir, a.log, a.drop_whole_event)


if __name__ == "__main__":
    main()
