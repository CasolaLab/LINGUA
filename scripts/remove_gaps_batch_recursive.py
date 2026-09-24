#!/usr/bin/env python3
"""

python3 remove_gaps_batch_recursive.py -i brassicaceae_23species_LINGUA_results -o brassicaceae_23species_LINGUA_nogaps

remove_gaps_batch_recursive.py (recursive, per-folder outputs)

Given a master input folder that contains multiple subfolders (one per species),
this script finds folders containing *.aln files and, for each such folder,
writes ungapped FASTA (*.fna) files into a corresponding output folder under a
master output folder.

Example folder naming
---------------------
Input subfolder name:
  Bstricta_DNG_45k_cds_stitched_exons
Output subfolder name (created under -o output_masterfolder):
  bstricta_stitched_nogaps

Rule: output folder name = <prefix>_stitched_nogaps
where <prefix> is the part of the input folder name before the first underscore,
lowercased (e.g., "Bstricta" -> "bstricta", "Ahalleri" -> "ahalleri").

File naming
-----------
- Input files are expected to end with: .aln
- Output filename = original name with "_nogaps" appended before the extension
- The ".aln" extension is replaced with ".fna"

Example:
  Brassica_oleracea_Bo00285s020.1_Scaffold00285_7031_8541.aln
->Brassica_oleracea_Bo00285s020.1_Scaffold00285_7031_8541_nogaps.fna

Usage
-----
python3 remove_gaps_batch_recursive.py -i masterfolder -o output_masterfolder
"""

from __future__ import annotations

from pathlib import Path
import argparse
import sys
from typing import Iterator, Tuple, List, Set


def iter_fasta_records(text: str) -> Iterator[Tuple[str, str]]:
    header = None
    seq_chunks: List[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_chunks)
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


def wrap_seq(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def process_file(in_path: Path, out_dir: Path) -> None:
    out_name = in_path.name[:-4] + "_nogaps.fna"  # strip ".aln"
    out_path = out_dir / out_name

    text = in_path.read_text(errors="replace")

    out_lines: List[str] = []
    for hdr, seq in iter_fasta_records(text):
        seq = "".join(seq.split())
        seq = seq.replace("-", "")
        out_lines.append(f">{hdr}")
        out_lines.append(wrap_seq(seq))

    if not out_lines:
        raise ValueError(f"No FASTA records found in {in_path}")

    out_path.write_text("\n".join(out_lines) + "\n")


def output_subfolder_name(input_folder_name: str) -> str:
    prefix = input_folder_name.split("_", 1)[0].strip()
    if not prefix:
        prefix = "unknown"
    return f"{prefix.lower()}_stitched_nogaps"


def find_aln_parent_dirs(master_in: Path) -> List[Path]:
    """
    Return directories that contain at least one *.aln file.
    Keep only the shallowest such directories (avoid double-processing nested dirs).
    """
    parents: Set[Path] = set()
    for aln in master_in.rglob("*.aln"):
        parents.add(aln.parent)

    if not parents:
        return []

    parents_list = sorted(parents, key=lambda p: (len(p.parts), str(p)))
    shallow: List[Path] = []
    for p in parents_list:
        # skip p if it's inside an already-chosen shallow dir
        if any(p.is_relative_to(s) and p != s for s in shallow):
            continue
        shallow.append(p)
    return shallow


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Remove '-' gaps from FASTA alignment files (*.aln) in all subfolders under a master folder."
    )
    ap.add_argument("-i", "--input", required=True, help="Master input folder containing subfolders with .aln files.")
    ap.add_argument("-o", "--output", required=True, help="Master output folder to create per-subfolder outputs.")
    ap.add_argument(
        "--recursive-within-species",
        action="store_true",
        help="Also process *.aln files in nested subfolders within each detected species folder (default: only that folder).",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files if present (default: skip if already exists).",
    )
    args = ap.parse_args()

    master_in = Path(args.input).expanduser().resolve()
    if not master_in.is_dir():
        print(f"ERROR: not a directory: {master_in}", file=sys.stderr)
        sys.exit(2)

    master_out = Path(args.output).expanduser().resolve()
    master_out.mkdir(parents=True, exist_ok=True)

    species_dirs = find_aln_parent_dirs(master_in)
    if not species_dirs:
        print(f"No .aln files found under {master_in}", file=sys.stderr)
        sys.exit(1)

    total_aln = 0
    ok = 0
    fail = 0
    skipped = 0

    for sp_dir in species_dirs:
        out_sub = master_out / output_subfolder_name(sp_dir.name)
        out_sub.mkdir(parents=True, exist_ok=True)

        aln_files = list(sp_dir.rglob("*.aln")) if args.recursive_within_species else list(sp_dir.glob("*.aln"))
        if not aln_files:
            continue

        for f in aln_files:
            total_aln += 1
            out_name = f.name[:-4] + "_nogaps.fna"
            out_path = out_sub / out_name
            if out_path.exists() and not args.overwrite:
                skipped += 1
                continue
            try:
                process_file(f, out_sub)
                ok += 1
            except Exception as e:
                fail += 1
                print(f"[FAIL] {f}: {e}", file=sys.stderr)

    print("Gap removal completed.")
    print(f"Master input folder : {master_in}")
    print(f"Master output folder: {master_out}")
    print(f"Species folders found: {len(species_dirs)}")
    print(f".aln files seen      : {total_aln}")
    print(f"Processed           : {ok}")
    print(f"Skipped (exists)    : {skipped}")
    print(f"Failed              : {fail}")

    if fail > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()