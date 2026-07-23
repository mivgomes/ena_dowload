#!/usr/bin/env bash
# Download the FASTQ files for an ENA study from the ENA filereport, with
# resume (-c) and md5 verification against ENA's published checksums.
#
# For PRJEB18629 each run has three files on ENA:
#   ERR*.fastq.gz     <- the generated merged file
#   ERR*_1.fastq.gz   <- paired forward
#   ERR*_2.fastq.gz   <- paired reverse
#
# Default is to downlaod all, but you can use WHICH=merged or WHICH=paired if only interested on those
#
# Usage:
#   bash download_ena.sh
#   bash download_ena.sh PRJEB18629 input merged
#
# Requires: curl, aria2, md5sum 

set -euo pipefail

ACCESSION="${1:-${ACCESSION:-PRJEB18629}}"
OUTDIR="${2:-${OUTDIR:-input}}"
WHICH="${3:-${WHICH:-all}}"   # merged | paired | all

mkdir -p "$OUTDIR"
REPORT="$OUTDIR/.ena_report.${ACCESSION}.tsv"

echo ">> Fetching ENA filereport for $ACCESSION"
curl -sSf "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=run_accession,fastq_md5,fastq_ftp&format=tsv" \
    -o "$REPORT"

# want_file <run> <basename> -> 0 if this file should be downloaded per WHICH
want_file() {
    local run="$1" base="$2"
    case "$WHICH" in
        all) return 0 ;;
        merged) [[ "$base" == "${run}.fastq.gz" ]] ;;
        paired) [[ "$base" == "${run}_1.fastq.gz" || "$base" == "${run}_2.fastq.gz" ]] ;;
        *) echo "Unknown WHICH='$WHICH' (use merged|paired|all)" >&2; exit 2 ;;
    esac
}

fail=0
# Skip header; each line: run_accession \t md5;md5;... \t url;url;...
tail -n +2 "$REPORT" | while IFS=$'\t' read -r run md5field ftpfield; do
    [[ -z "$ftpfield" ]] && { echo "!! $run: no fastq_ftp (skipped)"; continue; }
    IFS=';' read -ra md5s <<< "$md5field"
    IFS=';' read -ra ftps <<< "$ftpfield"
    for i in "${!ftps[@]}"; do
        url="${ftps[$i]}"; md5_exp="${md5s[$i]}"
        base="$(basename "$url")"
        want_file "$run" "$base" || continue
        dest="$OUTDIR/$base"

        # already present and correct? skip.
        if [[ -f "$dest" ]] && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5_exp" ]]; then
            echo "== $base (already verified)"
            continue
        fi

        echo ">> downloading $base"
        # wget -c -q --show-progress -O "$dest" "https://$url"
        aria2c -x16 -s16 -c -o "$dest" "https://$url"

        md5_got="$(md5sum "$dest" | cut -d' ' -f1)"
        if [[ "$md5_got" == "$md5_exp" ]]; then
            echo "OK  $base"
        else
            echo "MD5_FAIL $base (expected $md5_exp got $md5_got)"
            fail=1
        fi
    done
done

echo ">> Done. Files in $OUTDIR/"
[[ "$fail" -eq 0 ]] || { echo "Some files failed md5 - re-run to resume/re-download."; exit 1; }
