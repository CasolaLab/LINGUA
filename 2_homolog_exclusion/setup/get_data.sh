#!/usr/bin/env bash
# get_data.sh - download the public reference data some modules need. OPTIONAL: you can get the same files by hand
# (see setup/data_sources.txt). Nothing here is part of this repository, and nothing else in the tool calls this script.
# Downloading is the only thing in the whole tool that needs the internet; the searches themselves run offline.
#
# Usage:   bash setup/get_data.sh cdd          [DEST] [--dry-run] [--force] [--keep-archive]
#          bash setup/get_data.sh interproscan [DEST] [--dry-run] [--force] [--keep-archive] [--version X.Y-Z.0]
#          bash setup/get_data.sh taxonomy     [DEST] [--force] [--keep-archive]      (see setup/get_taxonomy.sh)
#
#   cdd           NCBI Conserved Domain Database, prebuilt for rpsblast (domain-cdd). 1.75 GB download (2026-09-26),
#                 3.8 GB unpacked. Default DEST: ./data/cdd. Point --db at DEST/Cdd.
#   interproscan  InterProScan with all its data (domain-interpro). 7.04 GB download for 5.78-109.0, about 36.5 GB
#                 unpacked. Default DEST: ./data. Needs Java 11 or newer to run (not installed by this script).
#                 EBI can throttle a single connection heavily; if the download is very slow, use a
#                 multi-connection downloader on the same URL (see setup/data_sources.txt).
#   taxonomy      NCBI taxonomy dump for the optional db-report (80 MB download).
#   --dry-run     show what would be downloaded and how big it is, and download nothing
# A partly downloaded file is resumed if you run the same command again.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WHAT="${1:-}"
[[ -z "$WHAT" || "$WHAT" == "-h" || "$WHAT" == "--help" ]] && { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0; }
shift
DEST=""
DRY=0; FORCE=0; KEEP=0; IPS_VERSION="5.78-109.0"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY=1 ;;
        --force) FORCE=1 ;;
        --keep-archive) KEEP=1 ;;
        --version) shift; IPS_VERSION="${1:?--version needs a value, for example 5.78-109.0}" ;;
        -*) echo "Error: unknown option '$1'." >&2; exit 1 ;;
        *) DEST="$1" ;;
    esac
    shift
done

command -v curl >/dev/null 2>&1 || { echo "Error: curl not found." >&2; exit 1; }
command -v tar >/dev/null 2>&1 || { echo "Error: tar not found." >&2; exit 1; }

remote_size() { { curl -fsIL --max-time 30 "$1" 2>/dev/null || true; } | tr -d '\r' | awk -F': ' 'tolower($1)=="content-length"{v=$2} END{print v}'; }
human() { awk -v b="$1" 'BEGIN{ if (b == "" || b == 0) print "size unknown"; else if (b < 1e9) printf "%.0f MB", b/1e6; else printf "%.2f GB", b/1e9 }'; }
md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1
    elif command -v md5 >/dev/null 2>&1; then md5 -q "$1"
    else echo "Error: no md5sum or md5 found; cannot verify the download." >&2; exit 1; fi
}
download() {   # download URL OUTPUT  (resumes)
    curl -fL -C - --retry 3 --retry-delay 5 -sS --progress-bar -o "$2" "$1"
}
free_gb() { df -Pk "$1" | awk 'NR==2{printf "%.0f", $4/1048576}'; }

case "$WHAT" in
taxonomy)
    args=()
    [[ -n "$DEST" ]] && args+=("$DEST")
    [[ "$FORCE" -eq 1 ]] && args+=(--force)
    [[ "$KEEP" -eq 1 ]] && args+=(--keep-archive)
    [[ "$DRY" -eq 1 ]] && { echo "Would download https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz ($(human "$(remote_size https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)"))"; exit 0; }
    exec bash "$SCRIPT_DIR/get_taxonomy.sh" "${args[@]}"
    ;;
cdd)
    DEST="${DEST:-data/cdd}"
    URL="https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz"
    INFO="https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd.info"
    size="$(remote_size "$URL")"
    [[ -n "$size" ]] || { echo "Error: could not reach $URL (no internet, or the file moved). See setup/data_sources.txt." >&2; exit 1; }
    echo "CDD (prebuilt for rpsblast): $URL"
    echo "  download: $(human "$size"); unpacked: about 3.8 GB; needs about 6 GB free while unpacking."
    if [[ "$DRY" -eq 1 ]]; then echo "  (dry run: nothing downloaded)"; exit 0; fi
    if [[ -f "$DEST/Cdd.pal" && "$FORCE" -eq 0 ]]; then
        echo "Error: $DEST already holds Cdd.pal. Use --force to download again." >&2; exit 1
    fi
    mkdir -p "$DEST"
    echo "  free space in $DEST: $(free_gb "$DEST") GB"
    archive="$DEST/Cdd_LE.tar.gz"
    download "$URL" "$archive"
    got="$(wc -c < "$archive" | tr -d ' ')"
    if [[ "$got" != "$size" ]]; then
        echo "Error: the download is $got bytes, NCBI lists $size. Run the same command again to resume." >&2; exit 1
    fi
    # NCBI publishes no checksum next to this file, so check the size (above) and that the archive reads to its end.
    tar -tzf "$archive" >/dev/null || { echo "Error: the archive is damaged; delete $archive and run again." >&2; exit 1; }
    curl -fsSL --retry 3 -o "$DEST/cdd.info" "$INFO" || echo "Warning: could not fetch cdd.info (the CDD release will be recorded as unknown)." >&2
    tar -xzf "$archive" -C "$DEST"
    [[ "$KEEP" -eq 1 ]] || rm -f "$archive"
    {
        echo "source=$URL"
        echo "downloaded_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        echo "archive_bytes=$size"
        [[ -f "$DEST/cdd.info" ]] && echo "release=$(head -n1 "$DEST/cdd.info")"
    } > "$DEST/SOURCE.txt"
    echo "Done: use  --db $DEST/Cdd"
    [[ -f "$DEST/Cdd.pal" ]] || echo "Warning: $DEST/Cdd.pal was not found after unpacking; check the folder." >&2
    ;;
interproscan)
    DEST="${DEST:-data}"
    V="$IPS_VERSION"
    BASE="https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/5/$V"
    NAME="interproscan-$V-64-bit.tar.gz"
    size="$(remote_size "$BASE/$NAME")"
    [[ -n "$size" ]] || { echo "Error: $BASE/$NAME was not found. Check the version (see setup/data_sources.txt)." >&2; exit 1; }
    echo "InterProScan $V: $BASE/$NAME"
    echo "  download: $(human "$size"); unpacked: about 36.5 GB (version 5.78); needs about 45 GB free while unpacking."
    echo "  Java 11 or newer is needed to run it (not installed by this script)."
    if [[ "$DRY" -eq 1 ]]; then echo "  (dry run: nothing downloaded)"; exit 0; fi
    if [[ -d "$DEST/interproscan-$V" && "$FORCE" -eq 0 ]]; then
        echo "Error: $DEST/interproscan-$V already exists. Use --force to download again." >&2; exit 1
    fi
    mkdir -p "$DEST"
    echo "  free space in $DEST: $(free_gb "$DEST") GB"
    archive="$DEST/$NAME"
    download "$BASE/$NAME" "$archive"
    curl -fsSL --retry 3 -o "$archive.md5" "$BASE/$NAME.md5"
    expected="$(cut -d' ' -f1 "$archive.md5" | head -n1)"
    actual="$(md5_of "$archive")"
    if [[ "$expected" != "$actual" ]]; then
        echo "Error: the checksum does not match (expected $expected, got $actual). Delete $archive and run again." >&2; exit 1
    fi
    echo "Checksum OK ($actual). Unpacking (this takes a while) ..."
    tar -xzf "$archive" -C "$DEST"
    [[ "$KEEP" -eq 1 ]] || rm -f "$archive" "$archive.md5"
    echo "source=$BASE/$NAME" > "$DEST/interproscan-$V/SOURCE.txt"
    echo "archive_md5=$actual" >> "$DEST/interproscan-$V/SOURCE.txt"
    echo "downloaded_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$DEST/interproscan-$V/SOURCE.txt"
    echo "Done: $DEST/interproscan-$V/interproscan.sh  (use --iprscan-bin, or put its folder on PATH)"
    ;;
*)
    echo "Error: unknown data set '$WHAT' (use cdd, interproscan or taxonomy)." >&2; exit 1 ;;
esac
