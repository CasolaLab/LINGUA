#!/usr/bin/env python3
"""
==============================================================================
Script: step9_exon_coordinate_trim.py
Description: Step 9 trimming by exon coordinates, for any species. Works out two things
from the data instead of assuming them:

  1. BED convention, once per run: 0-based half-open (exon length = end - start, as in
     Athaliana.bed) or 1-based inclusive (end - start + 1, as in reformatted_beds/),
     whichever makes BED exon lengths equal the FNA lengths for more transcripts.
  2. What step 1 extracted, per file, from the reference row's length:
       CDS + 1 base per exon -> the extra base is removed from every exon block (last
                                base of each block on + strand genes, first base on -
                                strand genes), located from the BED coordinates, so
                                every species is cut at the same, correct columns
       CDS length exactly    -> nothing to remove

Either way, a file is written only if the trimmed reference row equals its FNA sequence
exactly. A reference row that has the CDS length but every exon shifted one base along
the genome (one exon base missing, one outside base added) cannot be repaired here,
because the missing base was never extracted; those files are refused as
EXTRACTION_SHIFTED and step 1 has to be rerun with corrected coordinates.

Columns where the reference row has '-' are removed, so every species is projected
onto the reference CDS (as in the earlier step 9 versions). The report records how many
such columns, and how many non-reference bases in them, were removed per file.

A file is written only if, after trimming, the reference row equals its FNA sequence
exactly. Every other file is not written and is listed with its reason in
<output_dir>/step9_trimming_report.tsv, together with every written file.

The reference row is the record whose species name (header text before the first '.')
equals the species prefix of the step 8 file name
(<species>_<transcriptID>_<chr>_<start>_<end>.aln). The transcript is identified by an
ID in the file name that is in the BED and whose exon start/end match the file name.

Usage: python3 scripts/step9_exon_coordinate_trim.py <step8_stitched_dir> <output_dir> <reference_nucleotide_fna_file> <bed_file> <threads>
==============================================================================
"""

import concurrent.futures
import datetime
import hashlib
import os
import sys

REPORT_NAME = "step9_trimming_report.tsv"
STATUSES = [
    "OK",
    "DUPLICATE_FILENAME",
    "READ_ERROR",
    "NO_RECORDS",
    "TEXT_BEFORE_FIRST_HEADER",
    "ROW_LENGTHS_DIFFER",
    "NO_TRANSCRIPT_MATCH",
    "NOT_IN_FNA",
    "DUPLICATE_FNA_ID",
    "BED_FNA_LENGTH_MISMATCH",
    "NO_REFERENCE_ROW",
    "MULTIPLE_REFERENCE_ROWS",
    "UNEXPECTED_REFERENCE_LENGTH",
    "EXTRACTION_SHIFTED",
    "REFERENCE_MISMATCH_AFTER_TRIM",
    "WRITE_ERROR",
]
REPORT_COLUMNS = [
    "file", "transcript", "strand", "n_exons", "n_species", "status", "detail",
    "bed_convention", "reference_extraction", "alignment_columns",
    "extra_base_columns_removed_1based", "ref_gap_columns_removed",
    "nonref_bases_in_removed_ref_gap_columns", "nonref_bases_in_extra_base_columns",
]

_bed = {}
_fna = {}
_fna_ambiguous = set()
_convention = "0-based"


def log(msg):
    print(msg, file=sys.stderr)


def load_bed(path):
    tx = {}
    with open(path) as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 7:
                sys.exit(f"Error: {path} line {n}: expected at least 7 tab-separated columns")
            try:
                start, end = int(c[1].strip()), int(c[2].strip())
            except ValueError:
                if n == 1:
                    continue
                sys.exit(f"Error: {path} line {n}: chromStart/chromEnd are not integers")
            chrom, strand, tid = c[0].strip(), c[5].strip(), c[6].strip()
            if strand not in ("+", "-"):
                sys.exit(f"Error: {path} line {n}: strand must be + or -, got {strand!r}")
            if end <= start:
                sys.exit(f"Error: {path} line {n}: chromEnd <= chromStart")
            t = tx.setdefault(tid, {"chrom": chrom, "strand": strand, "exons": []})
            if (t["chrom"], t["strand"]) != (chrom, strand):
                sys.exit(f"Error: {path}: transcript {tid} has exons on different chromosomes or strands")
            t["exons"].append((start, end))
    for t in tx.values():
        # Same exon order step 8 used to stitch: ascending start on +, descending on -.
        t["exons"].sort(reverse=(t["strand"] == "-"))
    return tx


def load_fna(path):
    seqs, ambiguous = {}, set()
    seq_id, parts = None, []

    def store():
        if seq_id is None:
            return
        s = "".join(parts).upper()
        if seq_id in seqs and seqs[seq_id] != s:
            ambiguous.add(seq_id)
        seqs.setdefault(seq_id, s)

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                store()
                fields = line[1:].split()
                seq_id, parts = (fields[0] if fields else ""), []
            elif seq_id is not None:
                parts.append("".join(line.split()))
    store()
    return seqs, ambiguous


def init_worker(bed, fna, fna_ambiguous, convention):
    global _bed, _fna, _fna_ambiguous, _convention
    _bed, _fna, _fna_ambiguous, _convention = bed, fna, fna_ambiguous, convention


def exon_lengths(exons, convention):
    return [e - s + (1 if convention == "1-based" else 0) for s, e in exons]


def detect_convention(bed, fna):
    """Count transcripts whose BED exon lengths equal the FNA length under each reading."""
    fits = {"0-based": 0, "1-based": 0}
    for tid, t in bed.items():
        if tid in fna:
            for conv in fits:
                if sum(exon_lengths(t["exons"], conv)) == len(fna[tid]):
                    fits[conv] += 1
    convention = "1-based" if fits["1-based"] > fits["0-based"] else "0-based"
    return convention, fits


def is_shifted(ref_bases, cds, lengths, strand):
    """True if every exon block is the exon moved one base toward higher genome
    coordinates: on + it lacks the exon's first base and ends with an outside base; on -
    (reverse complemented) it starts with an outside base and lacks the last base."""
    o = 0
    for L in lengths:
        b, x = ref_bases[o:o + L], cds[o:o + L]
        o += L
        if not (b[1:] == x[:-1] if strand == "-" else b[:-1] == x[1:]):
            return False
    return True


def match_transcript(stem):
    """Returns (species, transcript_id) or (None, reason)."""
    toks = stem.split("_")
    if len(toks) < 5:
        return None, "file name is not <species>_<transcriptID>_<chr>_<start>_<end>.aln"
    try:
        f_start, f_end = int(toks[-2]), int(toks[-1])
    except ValueError:
        return None, "file name does not end in _<start>_<end>"
    found = []
    for i in range(1, len(toks) - 3):
        for j in range(i + 1, len(toks) - 2):
            tid = "_".join(toks[i:j])
            t = _bed.get(tid)
            if t is None:
                continue
            starts = [s for s, _ in t["exons"]]
            ends = [e for _, e in t["exons"]]
            if min(starts) == f_start and max(ends) == f_end:
                found.append(("_".join(toks[:i]), tid))
    if len(found) == 1:
        return found[0]
    if not found:
        return None, (f"no BED transcript in the file name has exons spanning {f_start}-{f_end} "
                      "(not in this BED, or step 8 stitched only some of its exons)")
    return None, "file name matches more than one BED transcript: " + ", ".join(t for _, t in found)


def read_alignment(path):
    with open(path, "rb") as f:
        text = f.read().decode("latin-1")
    headers, rows, preamble = [], [], False
    for line in text.split("\n"):
        line = line.rstrip("\r")
        if line.startswith(">"):
            headers.append(line)
            rows.append([])
        elif rows:
            rows[-1].append("".join(line.split()))
        elif line.strip():
            preamble = True
    return headers, ["".join(r) for r in rows], preamble


def wrap60(seq):
    if not seq:
        return [""]
    return [seq[i:i + 60] for i in range(0, len(seq), 60)]


def trim_file(path, output_dir):
    name = os.path.basename(path)
    row = dict.fromkeys(REPORT_COLUMNS, "NA")
    row["file"] = name

    def refuse(status, detail):
        row["status"], row["detail"] = status, detail
        return row

    try:
        headers, rows, preamble = read_alignment(path)
    except OSError as e:
        return refuse("READ_ERROR", e.strerror)
    if not headers:
        return refuse("NO_RECORDS", "no '>' records")
    row["n_species"] = str(len(headers))
    if preamble:
        return refuse("TEXT_BEFORE_FIRST_HEADER", "non-blank text before the first '>' record")
    width = len(rows[0])
    row["alignment_columns"] = str(width)
    if any(len(r) != width for r in rows):
        return refuse("ROW_LENGTHS_DIFFER", "rows have lengths " + ",".join(sorted({str(len(r)) for r in rows})))

    species, tid = match_transcript(name[:-len(".aln")])
    if species is None:
        return refuse("NO_TRANSCRIPT_MATCH", tid)
    t = _bed[tid]
    exons = t["exons"]
    row.update(transcript=tid, strand=t["strand"], n_exons=str(len(exons)))
    if tid not in _fna:
        return refuse("NOT_IN_FNA", f"{tid} not in FNA")
    if tid in _fna_ambiguous:
        return refuse("DUPLICATE_FNA_ID", f"{tid} appears in the FNA more than once with different sequences")
    cds = _fna[tid]
    row["bed_convention"] = _convention
    lengths = exon_lengths(exons, _convention)
    if sum(lengths) != len(cds):
        return refuse("BED_FNA_LENGTH_MISMATCH",
                      f"BED exon lengths ({_convention}) sum to {sum(lengths)}, FNA sequence is {len(cds)}")

    ref_rows = [k for k, h in enumerate(headers) if (h[1:].split() or [""])[0].split(".")[0] == species]
    if not ref_rows:
        return refuse("NO_REFERENCE_ROW", f"no record named {species}.<chr>")
    if len(ref_rows) > 1:
        return refuse("MULTIPLE_REFERENCE_ROWS", f"{len(ref_rows)} records named {species}.<chr>")
    ref = rows[ref_rows[0]]

    base_cols = [i for i, ch in enumerate(ref) if ch != "-"]
    if len(base_cols) == len(cds) + len(exons):
        row["reference_extraction"] = "one_extra_per_exon"
        extra_idx, offset = [], 0
        for length in lengths:
            block = length + 1
            extra_idx.append(offset if t["strand"] == "-" else offset + block - 1)
            offset += block
    elif len(base_cols) == len(cds):
        row["reference_extraction"] = "cds_length"
        extra_idx = []
    else:
        return refuse("UNEXPECTED_REFERENCE_LENGTH",
                      f"reference has {len(base_cols)} bases; expected {len(cds)} (CDS) or "
                      f"{len(cds) + len(exons)} (CDS + one extra per exon)")
    extra_cols = [base_cols[u] for u in extra_idx]
    extra_set = set(extra_cols)
    keep_cols = [c for c in base_cols if c not in extra_set]

    trimmed_ref = "".join(ref[c] for c in keep_cols).upper()
    if trimmed_ref != cds:
        diffs = [i for i, (a, b) in enumerate(zip(trimmed_ref, cds)) if a != b]
        if not extra_idx and is_shifted(trimmed_ref, cds, lengths, t["strand"]):
            return refuse("EXTRACTION_SHIFTED",
                          f"every exon is shifted one base along the genome ({len(diffs)} of {len(cds)} "
                          "positions differ); one exon base per exon was never extracted, so this "
                          "cannot be trimmed - rerun step 1 with corrected coordinates")
        return refuse("REFERENCE_MISMATCH_AFTER_TRIM",
                      f"{len(diffs)} of {len(cds)} positions differ from the FNA, first at CDS position "
                      f"{diffs[0] + 1}")

    ref_gap_cols = [i for i, ch in enumerate(ref) if ch == "-"]
    others = [r for k, r in enumerate(rows) if k != ref_rows[0]]
    row["extra_base_columns_removed_1based"] = ",".join(str(c + 1) for c in extra_cols)
    row["ref_gap_columns_removed"] = str(len(ref_gap_cols))
    row["nonref_bases_in_removed_ref_gap_columns"] = str(sum(r[c] != "-" for r in others for c in ref_gap_cols))
    row["nonref_bases_in_extra_base_columns"] = str(sum(r[c] != "-" for r in others for c in extra_cols))

    out = []
    for h, r in zip(headers, rows):
        out.append(h)
        out.extend(wrap60("".join(r[c] for c in keep_cols)))
    out_path = os.path.join(output_dir, name)
    tmp_path = out_path + ".partial"
    try:
        with open(tmp_path, "wb") as f:
            f.write(("\n".join(out) + "\n").encode("latin-1"))
        os.replace(tmp_path, out_path)
    except OSError as e:
        return refuse("WRITE_ERROR", e.strerror)
    row["status"], row["detail"] = "OK", "trimmed reference == FNA"
    return row


def find_aln_files(root):
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        for fn in sorted(filenames):
            p = os.path.join(dirpath, fn)
            if fn.endswith(".aln") and os.path.isfile(p) and not os.path.islink(p):
                found.append(p)
    return found


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    if len(sys.argv) != 6:
        print("Usage: python3 step9_exon_coordinate_trim.py <step8_stitched_dir> <output_dir> <reference_nucleotide_fna_file> <bed_file> <threads>")
        sys.exit(1)
    input_dir, output_dir, fna_file, bed_file = sys.argv[1:5]
    try:
        threads = int(sys.argv[5])
        assert threads >= 1
    except (ValueError, AssertionError):
        sys.exit("Error: threads must be a positive integer")

    if not os.path.isdir(input_dir):
        sys.exit(f"Error: step 8 directory not found: {input_dir}")
    for p in (bed_file, fna_file):
        if not os.path.isfile(p):
            sys.exit(f"Error: input file not found: {p}")
    if os.path.abspath(output_dir) == os.path.abspath(input_dir):
        sys.exit("Error: output_dir must differ from step8_stitched_dir")
    if os.path.isdir(output_dir):
        existing = [f for f in os.listdir(output_dir) if f.endswith(".aln") or f == REPORT_NAME]
        if existing:
            sys.exit(f"Error: {output_dir} already contains {len(existing)} .aln/report files; "
                     "use an empty or new output directory so earlier results are not mixed in or overwritten")
    os.makedirs(output_dir, exist_ok=True)

    log("Loading BED...")
    bed = load_bed(bed_file)
    log("Loading FNA...")
    fna, fna_ambiguous = load_fna(fna_file)
    log(f"BED: {len(bed)} transcripts; FNA: {len(fna)} sequences")
    convention, fits = detect_convention(bed, fna)
    log(f"BED convention: {convention} (transcripts whose exon lengths match the FNA: "
        f"0-based {fits['0-based']}, 1-based {fits['1-based']})")

    files = find_aln_files(input_dir)
    by_name = {}
    for p in files:
        by_name.setdefault(os.path.basename(p), []).append(p)
    rows, to_run = [], []
    for name, paths in by_name.items():
        if len(paths) > 1:
            for p in paths:
                r = dict.fromkeys(REPORT_COLUMNS, "NA")
                r.update(file=os.path.relpath(p, input_dir), status="DUPLICATE_FILENAME",
                         detail=f"{len(paths)} files named {name} under {input_dir}; none written")
                rows.append(r)
        else:
            to_run.append(paths[0])
    log(f"Trimming {len(to_run)} alignment files using {threads} threads...")

    chunk = max(1, len(to_run) // (threads * 8))
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads, initializer=init_worker,
                                                initargs=(bed, fna, fna_ambiguous, convention)) as ex:
        for k, r in enumerate(ex.map(trim_file, to_run, [output_dir] * len(to_run), chunksize=chunk), 1):
            rows.append(r)
            if k % 1000 == 0:
                log(f"Progress: {k} / {len(to_run)}")

    counts = {s: 0 for s in STATUSES}
    for r in rows:
        counts[r["status"]] += 1
    assert sum(counts.values()) == len(files) == len(rows)
    written = sum(f.endswith(".aln") for f in os.listdir(output_dir))
    assert written == counts["OK"], f"{written} .aln files in output_dir but {counts['OK']} OK"

    rows.sort(key=lambda r: (r["status"] != "OK", r["status"], r["file"]))
    with open(os.path.join(output_dir, REPORT_NAME), "w") as f:
        f.write(f"# script: {os.path.abspath(__file__)}\n")
        f.write(f"# date: {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"# step8_stitched_dir: {os.path.abspath(input_dir)}\n")
        f.write(f"# bed_file: {os.path.abspath(bed_file)} sha256={sha256(bed_file)}\n")
        f.write(f"# reference_nucleotide_fna_file: {os.path.abspath(fna_file)} sha256={sha256(fna_file)}\n")
        f.write(f"# bed_convention: {convention} (exon lengths matching FNA: 0-based {fits['0-based']}, "
                f"1-based {fits['1-based']})\n")
        f.write(f"# files_found: {len(files)}; " + "; ".join(f"{s}: {n}" for s, n in counts.items() if n) + "\n")
        f.write("\t".join(REPORT_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r[c] for c in REPORT_COLUMNS) + "\n")

    print(f"Trimming complete: {len(files)} .aln files found, {counts['OK']} written to {output_dir}")
    print(f"  BED convention used: {convention}")
    for mode in ("one_extra_per_exon", "cds_length"):
        n = sum(r["status"] == "OK" and r["reference_extraction"] == mode for r in rows)
        if n:
            print(f"  written, reference extraction {mode}: {n}")
    for s in STATUSES[1:]:
        if counts[s]:
            print(f"  not written, {s}: {counts[s]}")
    gap_files = sum(r["ref_gap_columns_removed"] not in ("NA", "0") for r in rows)
    if gap_files:
        print(f"  note: {gap_files} written files had columns where the reference has '-' (removed; see report)")
    if counts["EXTRACTION_SHIFTED"]:
        print(f"  WARNING: {counts['EXTRACTION_SHIFTED']} files have every exon shifted one base by step 1; "
              "they cannot be trimmed correctly - rerun step 1 with corrected coordinates")
    print(f"Per-file details: {os.path.join(output_dir, REPORT_NAME)}")


if __name__ == "__main__":
    main()
