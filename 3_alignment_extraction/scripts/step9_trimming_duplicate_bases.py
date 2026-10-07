#!/usr/bin/env python3
"""
==============================================================================
Script: step9_trimming_duplicate_bases.py
Description: Open Reading Frame Protector. Trims extraneous padding and intron bleed-over from stitched sequences by tightly anchoring the alignment against the true Reference FNA sequence. Guarantees the 3-bp codon frame is preserved.
Usage: python3 step9_trimming_duplicate_bases.py <step8_stitched_dir> <output_dir> <reference_nucleotide_fna_file> <bed_file> <threads>
==============================================================================
"""

import concurrent.futures
import fnmatch
import os
import re
import stat
import sys

AWK_FS_CHARS = " \t\n"
AWK_SPACE_CLASS = " \t\n\v\f\r"
ASCII_LOWER = str.maketrans("ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz")

_truth_by_id = {}


def read_records(path):
    # latin-1 maps each byte to one character, matching mawk's byte-oriented strings
    with open(path, "rb") as f:
        lines = f.read().decode("latin-1").split("\n")
    if lines[-1] == "":
        lines.pop()
    return lines


def awk_fields(line):
    stripped = line.strip(AWK_FS_CHARS)
    return re.split(r"[ \t\n]+", stripped) if stripped else []


def nth(fields, n):
    return fields[n - 1] if len(fields) >= n else ""


def ref_lookup_lines(ref_fasta):
    """The (id, sequence) lines the bash script wrote to ref_lookup.tmp."""
    try:
        lines = read_records(ref_fasta)
    except OSError as e:
        print(f'awk: cannot open "{ref_fasta}" ({e.strerror})', file=sys.stderr)
        return []
    out, seq_id, seq = [], "", []
    for line in lines:
        if line.startswith(">"):
            if seq:
                out.append(seq_id + "\t" + "".join(seq))
            seq_id = nth(awk_fields(line), 1)[1:]
            seq = []
            continue
        seq.append(line)
    out.append(seq_id + "\t" + "".join(seq))
    return out


def build_truth_index(lookup_lines):
    """ID -> lowercased sequence; for a repeated ID the first record wins, as in bash."""
    by_id = {}
    for line in lookup_lines:
        f = awk_fields(line)
        key = nth(f, 1)
        if key:
            by_id.setdefault(key, nth(f, 2).translate(ASCII_LOWER))
    return by_id


def find_transcript_id(file_name):
    """Returns (FNA ID found in the file name, None) or (None, reason).

    Equivalent to testing every FNA ID for `id in name` with whole-ID boundaries, but
    looks up only the name's substrings that start and end at a boundary."""
    stem = file_name[:-len(".aln")] if file_name.endswith(".aln") else file_name
    starts = [i for i in range(len(stem)) if i == 0 or not stem[i - 1].isalnum()]
    ends = [j for j in range(1, len(stem) + 1) if j == len(stem) or not stem[j].isalnum()]
    hits = {stem[i:j] for i in starts for j in ends if j > i and stem[i:j] in _truth_by_id}
    if not hits:
        return None, "no FNA ID found in the file name"
    longest = max(hits, key=len)
    if all(h in longest for h in hits):
        return longest, None
    return None, "file name contains more than one FNA ID: " + ", ".join(sorted(hits))


def init_worker(by_id):
    global _truth_by_id
    _truth_by_id = by_id


def find_aln_files(start):
    """Same files, same order as `find start -type f -name "*.aln"` (no symlink following)."""
    try:
        st = os.lstat(start)
    except OSError as e:
        print(f"find: '{start}': {e.strerror}", file=sys.stderr)
        return
    if stat.S_ISREG(st.st_mode):
        if fnmatch.fnmatchcase(os.path.basename(start.rstrip("/")), "*.aln"):
            yield start
    elif stat.S_ISDIR(st.st_mode):
        yield from walk_preorder(start.rstrip("/") or "/")


def walk_preorder(directory):
    try:
        entries = list(os.scandir(directory))
    except OSError as e:
        print(f"find: '{directory}': {e.strerror}", file=sys.stderr)
        return
    for entry in entries:
        path = directory + "/" + entry.name
        if entry.is_dir(follow_symlinks=False):
            yield from walk_preorder(path)
        elif entry.is_file(follow_symlinks=False) and fnmatch.fnmatchcase(entry.name, "*.aln"):
            yield path


def trim_file(aln_file, output_dir):
    """Returns (status, stderr_message or None)."""
    file_name = os.path.basename(aln_file)
    target_id, lookup_problem = find_transcript_id(file_name)
    out_path = output_dir + "/" + file_name

    # The shell truncated the output file before awk read the input.
    try:
        out = open(out_path, "wb")
    except OSError as e:
        return "OUTPUT_ERROR", f"bash: {out_path}: {e.strerror}"
    with out:
        truth_seq = _truth_by_id.get(target_id, "") if target_id else ""
        try:
            lines = read_records(aln_file)
        except OSError as e:
            return "INPUT_ERROR", f'awk: cannot open "{aln_file}" ({e.strerror})'

        seq_count, headers, seqs = 0, {}, {}
        for line in lines:
            if line.startswith(">"):
                seq_count += 1
                headers[seq_count] = line
                continue
            seqs.setdefault(seq_count, []).append(line)
        seqs = {k: "".join(v) for k, v in seqs.items()}

        if lookup_problem:
            status = "AMBIGUOUS_ID" if "more than one" in lookup_problem else "NO_TRUTH_MATCH"
            return status, f"Error: No truth sequence found for {file_name} ({lookup_problem})"
        if not truth_seq:
            return "NO_TRUTH_MATCH", f"Error: No truth sequence found for {target_id} (FNA sequence is empty)"

        ref_aln = seqs.get(1, "").rstrip(AWK_SPACE_CLASS)
        keep_cols = []
        t_idx = 0
        for i, ch in enumerate(ref_aln):
            if ch == "-":
                continue
            if ch.translate(ASCII_LOWER) == truth_seq[t_idx]:
                keep_cols.append(i)
                t_idx += 1
            if t_idx >= len(truth_seq):
                break

        parts = []
        for s in range(1, seq_count + 1):
            parts.append(headers[s] + "\n")
            current = seqs.get(s, "")
            final_str = "".join(current[i] for i in keep_cols if i < len(current))
            while len(final_str) > 60:
                parts.append(final_str[:60] + "\n")
                final_str = final_str[60:]
            parts.append(final_str + "\n")
        out.write("".join(parts).encode("latin-1"))

    # Greedy matching keeps only bases equal to the FNA in order, so the trimmed
    # reference always equals the start of the FNA; the only failure is running short.
    if t_idx < len(truth_seq):
        return f"TRUNCATED: Output ({t_idx}) != FNA ({len(truth_seq)})", None
    return "SUCCESS", None


def trim_group(paths, output_dir):
    """Files sharing a basename write the same output path, so they run in order."""
    return [(p, trim_file(p, output_dir)) for p in paths]


def main():
    if len(sys.argv) != 6:
        print("Usage: python3 step9_trimming_duplicate_bases_bash_method.py <step8_stitched_dir> <output_dir> <reference_nucleotide_fna_file> <bed_file> <threads>")
        sys.exit(1)

    input_dir, output_dir, reference_nucleotide_fna_file, bed_file = sys.argv[1:5]
    try:
        threads = int(sys.argv[5])
    except ValueError:
        print("Error: threads must be an integer")
        sys.exit(1)

    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        print(f"mkdir: cannot create directory '{output_dir}': {e.strerror}", file=sys.stderr)

    try:
        open(bed_file, "rb").close()
    except OSError as e:
        print(f'awk: cannot open "{bed_file}" ({e.strerror})', file=sys.stderr)

    print("Loading Truth Sequences from FNA...")
    by_id = build_truth_index(ref_lookup_lines(reference_nucleotide_fna_file))

    aln_files = list(find_aln_files(input_dir))
    groups = {}
    for p in aln_files:
        groups.setdefault(os.path.basename(p), []).append(p)
    total_files = len(aln_files)
    print(f"Loaded {len(by_id)} truth sequences.")
    print(f"Processing {total_files} alignment files using {threads} threads...")
    if len(groups) < total_files:
        print(f"Note: {total_files - len(groups)} files share a name with another file; "
              f"only the last of each name (in find order) remains in {output_dir}")

    counts = {"SUCCESS": 0, "TRUNCATED": 0, "NO_TRUTH_MATCH": 0, "AMBIGUOUS_ID": 0, "OTHER_ERROR": 0}
    results, next_i, completed = {}, 0, 0
    chunk = max(1, len(groups) // (max(threads, 1) * 8))
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads, initializer=init_worker,
                                                initargs=(by_id,)) as executor:
        for group_result in executor.map(trim_group, groups.values(),
                                         [output_dir] * len(groups), chunksize=chunk):
            results.update(group_result)
            # Report in find order, as the bash loop did.
            while next_i < total_files and aln_files[next_i] in results:
                path = aln_files[next_i]
                status, err = results.pop(path)
                next_i += 1
                completed += 1
                if err:
                    print(err, file=sys.stderr)
                if status == "SUCCESS":
                    counts["SUCCESS"] += 1
                elif status.startswith("TRUNCATED"):
                    counts["TRUNCATED"] += 1
                    print(f"[{os.path.basename(path)}]: {status}")
                elif status in ("NO_TRUTH_MATCH", "AMBIGUOUS_ID"):
                    counts[status] += 1
                else:
                    counts["OTHER_ERROR"] += 1
                if completed % 1000 == 0:
                    print(f"Progress: {completed} / {total_files}")
    assert completed == total_files == sum(counts.values())

    print("-" * 50)
    print(f"Trimming complete! {total_files} files processed.")
    print(f"  trimmed reference == FNA: {counts['SUCCESS']}")
    print(f"  trimmed reference shorter than FNA (written anyway, as bash did): {counts['TRUNCATED']}")
    print(f"  no FNA match (empty output file): {counts['NO_TRUTH_MATCH']}")
    print(f"  more than one FNA ID in the file name (empty output file): {counts['AMBIGUOUS_ID']}")
    print(f"  could not read input or write output: {counts['OTHER_ERROR']}")


if __name__ == "__main__":
    main()
