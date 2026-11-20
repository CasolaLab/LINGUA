#!/usr/bin/env bash
# Stand-in for InterProScan in the tutorial toy: accepts the real command line, writes an empty result and the
# completion lines the stage 2 module checks for. It never finds a domain.
out=""
while [[ $# -gt 0 ]]; do
  case "$1" in -o) out="$2"; shift 2;; --version) echo "InterProScan version 5.99-test"; exit 0;; *) shift;; esac
done
: > "$out"
echo "Running the following analyses:"; echo "[Pfam-38.2]"; echo "100% done:  InterProScan analyses completed"
