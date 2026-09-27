#!/usr/bin/env bash
# get_taxonomy.sh - download NCBI's taxonomy dump (names.dmp and nodes.dmp) for the OPTIONAL `db-report` command.
#
# This is the only step of the whole tool that needs the internet, and nothing else calls it. The searches never need it.
# The dump is public data from NCBI; it is not part of this repository and is not shipped with it. It changes daily, so
# download it yourself; `db-report` records which copy it used.
#
# Usage:   bash setup/get_taxonomy.sh [DESTINATION_FOLDER] [--force] [--keep-archive]
#          (default folder: ./taxonomy)
# Needs:   curl or wget, tar, and md5sum (or md5 on macOS). About 80 MB to download (79.6 MB on 2026-09-26; about 50 s on
#          a home connection); names.dmp (302 MB) and nodes.dmp (222 MB) are unpacked from it, so allow about 550 MB of
#          disk. The other files in the archive are not needed.

set -euo pipefail

BASE_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy"
DEST="taxonomy"
FORCE=0
KEEP=0
for arg in "$@"; do
    case "$arg" in
        -h|--help) sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        --force) FORCE=1 ;;
        --keep-archive) KEEP=1 ;;
        -*) echo "Error: unknown option '$arg'." >&2; exit 1 ;;
        *) DEST="$arg" ;;
    esac
done

fetch() {   # fetch URL OUTPUT
    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 --retry-delay 5 -sS -o "$2" "$1"
    elif command -v wget >/dev/null 2>&1; then
        wget -q -O "$2" "$1"
    else
        echo "Error: neither curl nor wget was found." >&2; exit 1
    fi
}
header_date() { # the archive's Last-Modified date, or "unknown"
    if command -v curl >/dev/null 2>&1; then
        curl -fsIL --max-time 30 "$1" | tr -d '\r' | awk -F': ' 'tolower($1)=="last-modified"{v=$2} END{print (v==""?"unknown":v)}'
    else
        echo "unknown"
    fi
}
md5_of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1
    elif command -v md5 >/dev/null 2>&1; then md5 -q "$1"
    else echo "Error: no md5sum or md5 found; cannot verify the download." >&2; exit 1; fi
}
command -v tar >/dev/null 2>&1 || { echo "Error: tar not found." >&2; exit 1; }

if [[ -f "$DEST/names.dmp" && -f "$DEST/nodes.dmp" && "$FORCE" -eq 0 ]]; then
    echo "Error: $DEST already holds names.dmp and nodes.dmp. Use --force to download again." >&2
    exit 1
fi
mkdir -p "$DEST"
WORK="$(mktemp -d "$DEST/.download.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

echo "Downloading $BASE_URL/taxdump.tar.gz (about 80 MB; about 550 MB once unpacked) ..."
LASTMOD="$(header_date "$BASE_URL/taxdump.tar.gz")"
fetch "$BASE_URL/taxdump.tar.gz" "$WORK/taxdump.tar.gz"
fetch "$BASE_URL/taxdump.tar.gz.md5" "$WORK/taxdump.tar.gz.md5"

expected="$(cut -d' ' -f1 "$WORK/taxdump.tar.gz.md5" | head -n1)"
actual="$(md5_of "$WORK/taxdump.tar.gz")"
if [[ "$expected" != "$actual" ]]; then
    echo "Error: the checksum does not match (expected $expected, got $actual). NCBI may have replaced the file while it" >&2
    echo "       was downloading; run this again. Nothing was unpacked." >&2
    exit 1
fi
echo "Checksum OK ($actual)."

tar -xzf "$WORK/taxdump.tar.gz" -C "$WORK" names.dmp nodes.dmp
mv -f "$WORK/names.dmp" "$DEST/names.dmp"
mv -f "$WORK/nodes.dmp" "$DEST/nodes.dmp"
if [[ "$KEEP" -eq 1 ]]; then mv -f "$WORK/taxdump.tar.gz" "$DEST/taxdump.tar.gz"; fi
cat > "$DEST/TAXONOMY_SOURCE.txt" <<EOF
source=$BASE_URL/taxdump.tar.gz
archive_last_modified=$LASTMOD
archive_md5=$actual
downloaded_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
echo "Done: $DEST/names.dmp, $DEST/nodes.dmp, $DEST/TAXONOMY_SOURCE.txt"
du -h "$DEST/names.dmp" "$DEST/nodes.dmp" | sed 's/^/  /'
