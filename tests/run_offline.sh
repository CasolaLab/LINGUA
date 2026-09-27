#!/usr/bin/env bash
# run_offline.sh - run the tests with the network cut off, to show the tool needs no internet.
# Linux and WSL only. Uses "unshare -rn": a private network namespace that has no connection to anything.
# It needs no sudo. If your system does not allow it, the script says so and stops.
#
# Usage: bash tests/run_offline.sh
# To run one of the modules against your own data with no network, put "unshare -rn" in front of the command:
#   unshare -rn python3 bin/domain_cdd.py -i proteins/ -o out/ --db /data/cdd/Cdd

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v unshare >/dev/null 2>&1; then
    echo "unshare not found; cannot cut off the network here." >&2
    exit 2
fi
if ! unshare -rn true 2>/dev/null; then
    echo "This system does not allow 'unshare -rn' (no permission to create a network namespace)." >&2
    exit 2
fi

echo "Checking that the network really is cut off inside the test environment:"
if unshare -rn bash -c 'exec 3<>/dev/tcp/1.1.1.1/53' 2>/dev/null; then
    echo "  a network connection succeeded, so the network is NOT cut off; stopping." >&2
    exit 3
fi
echo "  no connection possible: good."
echo

status=0
for t in "$HERE"/test_*.py; do
    printf '%s: ' "$(basename "$t")"
    out="$(unshare -rn python3 "$t" 2>&1)"; rc=$?
    echo "$(printf '%s\n' "$out" | grep -c '^PASS') checks passed, $(printf '%s\n' "$out" | grep -E 'ALL PASSED|FAILED')"
    if [[ $rc -ne 0 ]]; then
        status=1
        printf '%s\n' "$out" | grep -E '^FAIL|Traceback' | head -5
    fi
done
exit $status
