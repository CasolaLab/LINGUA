#!/usr/bin/env bash
# homolog_exclusion_pipeline.sh - entry point for LINGUA Stage 2 (homolog exclusion).
# A thin dispatcher: it picks a module and passes the remaining arguments to it.
# All real work is done by the Python scripts in bin/ (Python 3.11+, standard library only).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
BIN_DIR="$ROOT_DIR/bin"

usage() {
    cat <<'EOF'
Usage: homolog_exclusion_pipeline.sh <command> [options]

Commands:
  check                 Check that the required tools and databases are available (read-only)
  make-db               Build the species-tagged FASTA used as a homology database
  domain-interpro       Domain search with InterProScan
  domain-cdd            Domain search with RPS-BLAST against CDD
  homology-blast       Protein homology search with BLASTp (fast)
  homology-jackhmmer    Protein homology search with jackhmmer (sensitive)
  combine               Combine the results of several runs into one per-gene table and a verdict
  db-report             OPTIONAL: check how well a database covers the tree of life outside your clade
  make-bed              Look up CDS coordinates for a list of gene/protein IDs in a GFF3 (general purpose)

Run '<command> --help' for the options of one command.
Each module also runs on its own:  python3 bin/<module>.py --help
EOF
}

if [[ $# -lt 1 ]]; then
    usage
    exit 1
fi

cmd="$1"
shift

case "$cmd" in
    -h|--help|help)
        usage
        exit 0
        ;;
    -V|--version)
        exec python3 "$BIN_DIR/combine.py" --version
        ;;
    check)
        exec bash "$SCRIPT_DIR/check_tools.sh" "$@"
        ;;
    domain-interpro|domain-cdd|homology-blast|homology-jackhmmer|combine|make-db|db-report|make-bed)
        module_file="$BIN_DIR/${cmd//-/_}.py"
        [[ "$cmd" == "make-db" ]] && module_file="$BIN_DIR/make_homology_db.py"
        [[ "$cmd" == "make-bed" ]] && module_file="$BIN_DIR/make_cds_bed.py"
        if [[ ! -f "$module_file" ]]; then
            echo "Error: '$cmd' is not implemented yet ($module_file not found)." >&2
            exit 3
        fi
        command -v python3 >/dev/null 2>&1 || { echo "Error: python3 not found in PATH." >&2; exit 1; }
        exec python3 "$module_file" "$@"
        ;;
    *)
        echo "Error: unknown command '$cmd'." >&2
        usage >&2
        exit 1
        ;;
esac
