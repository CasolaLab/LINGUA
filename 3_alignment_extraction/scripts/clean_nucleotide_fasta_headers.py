import argparse
import sys

def load_bed_ids(bed_path):
    """Extract transcript IDs from the last column of a tab-separated BED file."""
    bed_ids = set()
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith(('#', 'track', 'browser')):
                cols = line.split('\t')
                bed_ids.add(cols[-1].strip())
    return bed_ids

def process_fasta(fasta_path, bed_path=None, output_path=None):
    bed_ids = load_bed_ids(bed_path) if bed_path else None
    out_handle = open(output_path, 'w') if output_path else sys.stdout

    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header_body = line[1:].strip()
                    first_token = header_body.split()[0]

                    # Match against BED IDs if provided
                    if bed_ids:
                        if first_token in bed_ids:
                            target_id = first_token
                        else:
                            # Fallback: search for a BED ID token within the full header
                            target_id = next((b_id for b_id in bed_ids if b_id in header_body), first_token)
                    else:
                        target_id = first_token

                    out_handle.write(f">{target_id}\n")
                else:
                    out_handle.write(line)
    finally:
        if output_path and out_handle != sys.stdout:
            out_handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim FASTA headers down to transcript IDs.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to input FASTA file")
    parser.add_argument("-b", "--bed", help="Optional path to BED file (uses last column for ID validation)")
    parser.add_argument("-o", "--output", help="Path to output FASTA file (defaults to terminal stdout)")
    
    args = parser.parse_args()
    process_fasta(args.fasta, args.bed, args.output)