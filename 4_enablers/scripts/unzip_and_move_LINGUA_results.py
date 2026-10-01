#!/usr/bin/env python3
"""
python unzip_and_move_LINGUA_results.py _all_species_data

unzip_and_move_LINGUA_results.py

Unzip all archive files in an input directory, then move selected LINGUA result
folders and summary TSV files into target directories created NEXT TO the input dir.

Usage:
    python unzip_and_move_LINGUA_results.py input_zip_dir

Optional:
    python unzip_and_move_LINGUA_results.py input_zip_dir --workdir unzipped_tmp --overwrite
"""

import argparse
import shutil
import zipfile
from pathlib import Path


RESULT_SUFFIX = "_DNG_45k_cds_LINGUA_results_version3"
INDELS_SUFFIX = "_DNG_45k_cds_stitched_exons_nomafft_final_indels_summary.tsv"
STOPS_SUFFIX = "_DNG_45k_cds_stitched_exons_nomafft_final_stop_codon_summary.tsv"


def unzip_file(zip_path: Path, outdir: Path) -> Path:
    extract_dir = outdir / zip_path.stem
    extract_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nUnzipping: {zip_path}")
    print(f"  To: {extract_dir}")

    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)

    return extract_dir


def move_path(src: Path, dest_dir: Path, overwrite: bool = False):
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / src.name

    if dest.exists():
        if overwrite:
            if dest.is_dir():
                shutil.rmtree(dest)
            else:
                dest.unlink()
            print(f"  Overwriting existing: {dest}")
        else:
            print(f"  SKIP, destination exists: {dest}")
            return

    shutil.move(str(src), str(dest))
    print(f"  Moved: {src} -> {dest}")
def process_extracted_dir(extract_dir: Path, base_output_dir: Path, overwrite: bool = False):
    lingua_outdir = base_output_dir / "brassicaceae_23species_LINGUA_results"
    indels_outdir = base_output_dir / "indels_dir"
    stops_outdir = base_output_dir / "stops_dir"

    outer_result_dirs = [
        p for p in extract_dir.rglob(f"*{RESULT_SUFFIX}*")
        if p.is_dir()
    ]

    if not outer_result_dirs:
        print(f"  WARNING: no result folders found in {extract_dir}")
        return

    for result_dir in outer_result_dirs:
        print(f"\nProcessing outer folder: {result_dir}")

        # Move the final stitched-exons folder
        final_dirs = [
            p for p in result_dir.rglob("*_DNG_45k_cds_stitched_exons_version3_nomafft_final")
            if p.is_dir()
        ]

        if final_dirs:
            for final_dir in final_dirs:
                move_path(final_dir, lingua_outdir, overwrite=overwrite)
        else:
            print(f"  WARNING: final stitched-exons folder not found in {result_dir}")

        # Move indels summary file
        indels_files = [
            p for p in result_dir.rglob(f"*{INDELS_SUFFIX}")
            if p.is_file()
        ]

        if indels_files:
            for indels_file in indels_files:
                move_path(indels_file, indels_outdir, overwrite=overwrite)
        else:
            print(f"  WARNING: indels summary not found in {result_dir}")

        # Move stop codon summary file
        stops_files = [
            p for p in result_dir.rglob(f"*{STOPS_SUFFIX}")
            if p.is_file()
        ]

        if stops_files:
            for stops_file in stops_files:
                move_path(stops_file, stops_outdir, overwrite=overwrite)
        else:
            print(f"  WARNING: stop codon summary not found in {result_dir}")

def main():
    parser = argparse.ArgumentParser(
        description="Unzip LINGUA result archives and move selected species result files/folders."
    )
    parser.add_argument("input_dir", help="Directory containing .zip files")
    parser.add_argument(
        "--workdir",
        default="unzipped_LINGUA_results_tmp",
        help="Extraction directory name (created inside input_dir parent)"
    )
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    base_output_dir = input_dir.parent.resolve()
    workdir = base_output_dir / args.workdir

    if not input_dir.exists():
        raise SystemExit(f"ERROR: input directory does not exist: {input_dir}")

    zip_files = sorted(input_dir.glob("*.zip"))

    if not zip_files:
        raise SystemExit(f"ERROR: no .zip files found in {input_dir}")

    print(f"Input dir:       {input_dir}")
    print(f"Output base dir:{base_output_dir}")
    print(f"Work dir:        {workdir}")

    workdir.mkdir(parents=True, exist_ok=True)

    for zip_path in zip_files:
        extract_dir = unzip_file(zip_path, workdir)
        process_extracted_dir(extract_dir, base_output_dir, overwrite=args.overwrite)

    print("\nDone.")


if __name__ == "__main__":
    main()