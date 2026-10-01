#!/usr/bin/env bash
# check_tools.sh - read-only check of the tools and databases the modules need.
# It never installs, downloads, or changes anything.
#
# Usage: check_tools.sh [--cdd-db PREFIX] [--blast-db PREFIX] [--hmmer-db FASTA] [--modules LIST]
#                       [--iprscan-bin PATH] [--applications LIST] [--skip-analysis-test]
#   --cdd-db PREFIX       path prefix of the formatted CDD database for rpsblast (e.g. /db/Cdd)
#   --blast-db PREFIX     path prefix of a formatted BLAST protein database
#   --hmmer-db FASTA      path of the FASTA file used as the jackhmmer database
#   --modules LIST        comma-separated subset to check (default: all)
#                         names: domain-interpro,domain-cdd,homology-blast,homology-jackhmmer
#   --iprscan-bin PATH    interproscan.sh to check (default: found on PATH)
#   --applications LIST   InterProScan analyses to check (default: the domain-interpro default list)
#   --skip-analysis-test  do not start InterProScan to test the analysis names (that test takes up to a minute)
# Exit status: 0 if everything requested was found, 1 otherwise.
# Nothing is downloaded and no network is used.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
CDD_DB=""
BLAST_DB=""
HMMER_DB=""
IPRSCAN_BIN=""
APPLICATIONS=""
SKIP_ANALYSIS_TEST=0
MODULES="domain-interpro,domain-cdd,homology-blast,homology-jackhmmer"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cdd-db) CDD_DB="${2:-}"; shift 2 ;;
        --blast-db) BLAST_DB="${2:-}"; shift 2 ;;
        --hmmer-db) HMMER_DB="${2:-}"; shift 2 ;;
        --modules) MODULES="${2:-}"; shift 2 ;;
        --iprscan-bin) IPRSCAN_BIN="${2:-}"; shift 2 ;;
        --applications) APPLICATIONS="${2:-}"; shift 2 ;;
        --skip-analysis-test) SKIP_ANALYSIS_TEST=1; shift ;;
        -h|--help) sed -n '2,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

problems=0
ok()   { printf '  [ ok ] %-22s %s\n' "$1" "$2"; }
bad()  { printf '  [MISS] %-22s %s\n' "$1" "$2"; problems=$((problems + 1)); }
info() { printf '  [info] %-22s %s\n' "$1" "$2"; }
wants() { [[ ",$MODULES," == *",$1,"* ]]; }

# InterProScan keeps each analysis's data in its own folder under data/. Coils and MobiDB-Lite have none.
ips_data_dir() {
    case "$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')" in
        pfam) echo pfam ;; prints) echo prints ;; smart) echo smart ;; ncbifam) echo ncbifam ;; sfld) echo sfld ;;
        superfamily) echo superfamily ;; gene3d) echo gene3d ;; hamap) echo hamap ;;
        prositeprofiles|prositepatterns) echo prosite ;; cdd) echo cdd ;; pirsr) echo pirsr ;; pirsf) echo pirsf ;;
        antifam) echo antifam ;; funfam) echo funfam ;; panther) echo panther ;;
        *) echo "" ;;
    esac
}

# Start InterProScan on a tiny sequence with the given analyses (-dp, so no network) and stop it as soon as it has either
# rejected a name or started the analyses. Prints the first rejected name, or nothing if the list was accepted.
ips_first_rejected() {
    local ipr="$1" list="$2" tmp log pid i
    tmp="$(mktemp -d)"; mkdir "$tmp/t"; log="$tmp/log"
    printf '>q\nMKVLAAGIVGLLLAQPAMAMKVLAAGIVGLLL\n' > "$tmp/q.faa"
    setsid "$ipr" -i "$tmp/q.faa" -f tsv -o "$tmp/o.tsv" -dp -T "$tmp/t" -appl "$list" > "$log" 2>&1 &
    pid=$!
    for i in $(seq 1 90); do
        sleep 1
        grep -q "does not exist or is deactivated" "$log" 2>/dev/null && break
        grep -q "Running the following analyses" "$log" 2>/dev/null && break
        kill -0 "$pid" 2>/dev/null || break
    done
    kill -TERM -- "-$pid" 2>/dev/null; sleep 1; kill -KILL -- "-$pid" 2>/dev/null
    grep -m1 -o "Analysis [^ ]* does not exist or is deactivated" "$log" | awk '{print $2}'
    rm -rf "$tmp"
}

check_ips_analyses() {
    local ipr="$1" apps="$2" datadir a d list="$2" name attempt
    local -a missing=() rejected=()
    datadir="$(dirname "$(readlink -f "$ipr")")/data"
    for a in ${apps//,/ }; do
        d="$(ips_data_dir "$a")"
        if [[ -n "$d" && ! -d "$datadir/$d" ]]; then missing+=("$a"); fi
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        bad "analysis data" "no data folder in $datadir for: ${missing[*]} (a run that includes these analyses would fail; an analysis name can be accepted even when its data is missing)"
    else
        ok "analysis data" "a data folder exists for every analysis that needs one"
    fi
    if [[ $SKIP_ANALYSIS_TEST -eq 1 ]]; then
        info "analysis names" "not tested (--skip-analysis-test)"
        return
    fi
    for attempt in $(seq 1 25); do
        name="$(ips_first_rejected "$ipr" "$list")"
        [[ -z "$name" ]] && break
        rejected+=("$name")
        list="$(printf '%s' "$list" | tr ',' '\n' | awk -v b="$name" 'tolower($0)!=tolower(b)' | paste -sd, -)"
        [[ -z "$list" ]] && break
    done
    if [[ ${#rejected[@]} -gt 0 ]]; then
        bad "analysis names" "not accepted by this InterProScan: ${rejected[*]} (they do not exist or are deactivated; names are case-insensitive, but MobiDB-Lite must be written MobiDBLite)"
    else
        ok "analysis names" "all accepted by this InterProScan"
    fi
}

echo "homolog-exclusion-pipeline: tool check"
echo

echo "Python"
if command -v python3 >/dev/null 2>&1; then
    pyver="$(python3 -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')"
    if python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else 1)'; then
        ok "python3" "$pyver ($(command -v python3))"
    else
        bad "python3" "$pyver found, but 3.11 or newer is required"
    fi
else
    bad "python3" "not found in PATH"
fi
echo

if wants domain-interpro; then
    echo "domain-interpro"
    ipr="${IPRSCAN_BIN:-$(command -v interproscan.sh 2>/dev/null || command -v iprscan 2>/dev/null || true)}"
    if [[ -n "$ipr" && -x "$ipr" ]]; then
        ver="$("$ipr" --version 2>&1 | grep -m1 -i -E 'interproscan|version' || true)"
        ok "interproscan" "${ver:-version unknown} ($ipr)"
        if command -v java >/dev/null 2>&1; then
            ok "java" "$(java -version 2>&1 | head -n1) (InterProScan needs 11 or newer)"
        else
            bad "java" "not found in PATH (InterProScan needs Java 11 or newer)"
        fi
        apps="${APPLICATIONS:-$(grep -m1 '^applications=' "$ROOT_DIR/defaults/domain-interpro.conf" 2>/dev/null | cut -d= -f2-)}"
        if [[ -n "$apps" ]]; then
            check_ips_analyses "$ipr" "$apps"
        fi
    elif [[ -n "$ipr" ]]; then
        bad "interproscan" "not executable: $ipr"
    else
        bad "interproscan" "neither interproscan.sh nor iprscan found in PATH (or give --iprscan-bin)"
    fi
    echo
fi

if wants domain-cdd; then
    echo "domain-cdd"
    if command -v rpsblast >/dev/null 2>&1; then
        ok "rpsblast" "$(rpsblast -version 2>&1 | head -n1) ($(command -v rpsblast))"
    else
        bad "rpsblast" "not found in PATH (BLAST+)"
    fi
    if [[ -n "$CDD_DB" ]]; then
        if compgen -G "${CDD_DB}.*" >/dev/null; then
            ok "CDD database" "$CDD_DB"
        else
            bad "CDD database" "no files matching ${CDD_DB}.*"
        fi
    else
        info "CDD database" "not checked (pass --cdd-db PREFIX)"
    fi
    echo
fi

if wants homology-blast; then
    echo "homology-blast"
    for t in blastp tblastn makeblastdb; do
        if command -v "$t" >/dev/null 2>&1; then
            ok "$t" "$("$t" -version 2>&1 | head -n1) ($(command -v "$t"))"
        else
            bad "$t" "not found in PATH (BLAST+)"
        fi
    done
    if [[ -n "$BLAST_DB" ]]; then
        dbdesc=""
        for dbt in prot nucl; do
            dbout="$(blastdbcmd -db "$BLAST_DB" -dbtype "$dbt" -info 2>/dev/null)" || continue
            n="$(printf '%s' "$dbout" | grep -m1 -o '[0-9,]* sequences')"
            [[ -n "$n" ]] && { dbdesc="$([[ $dbt == prot ]] && echo protein || echo nucleotide) database, $n"; break; }
        done
        if [[ -n "$dbdesc" ]]; then
            ok "BLAST database" "$BLAST_DB ($dbdesc; blastp needs protein, tblastn needs nucleotide)"
        elif compgen -G "${BLAST_DB}.[pn]*" >/dev/null; then
            ok "BLAST database" "$BLAST_DB (files found; blastdbcmd could not read them)"
        else
            bad "BLAST database" "no formatted BLAST database matching ${BLAST_DB} (build it with makeblastdb -parse_seqids)"
        fi
    else
        info "BLAST database" "not checked (pass --blast-db PREFIX)"
    fi
    echo
fi

if wants homology-jackhmmer; then
    echo "homology-jackhmmer"
    if command -v jackhmmer >/dev/null 2>&1; then
        # capture first: grep -m1 closes the pipe early, so its exit status must not decide the fallback
        hver="$(jackhmmer -h 2>&1 | grep -m1 '^# HMMER')"
        ok "jackhmmer" "${hver:-version unknown} ($(command -v jackhmmer))"
    else
        bad "jackhmmer" "not found in PATH (HMMER 3)"
    fi
    if [[ -n "$HMMER_DB" ]]; then
        if [[ -s "$HMMER_DB" ]]; then
            ok "HMMER database" "$HMMER_DB"
        else
            bad "HMMER database" "file missing or empty: $HMMER_DB"
        fi
    else
        info "HMMER database" "not checked (pass --hmmer-db FASTA)"
    fi
    echo
fi

echo "Optional"
if command -v parallel >/dev/null 2>&1; then
    info "GNU parallel" "found (not required; the modules run jobs themselves)"
else
    info "GNU parallel" "not found (not required)"
fi
echo

if [[ $problems -eq 0 ]]; then
    echo "All requested tools were found."
    exit 0
fi
echo "$problems item(s) missing. Install them (see docs/setup.md) and run check again."
exit 1
