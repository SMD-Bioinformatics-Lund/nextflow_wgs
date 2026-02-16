#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 1 ]]; then
    echo "Usage: $(basename "$0") [output_directory]" >&2
    exit 1
fi

for cmd in wget zgrep awk; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "ERROR: '${cmd}' command not found in PATH." >&2
        exit 1
    fi
done

OUTDIR="${1:-.}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

wget -O GRCh38.nr_deletions.pathogenic.tsv.gz "https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/deletions/GRCh38.nr_deletions.pathogenic.tsv.gz"
wget -O GRCh38.nr_duplications.pathogenic.tsv.gz "https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/duplications/GRCh38.nr_duplications.pathogenic.tsv.gz"
wget -O GRCh38.nr_insertions.pathogenic.tsv.gz "https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/insertions/GRCh38.nr_insertions.pathogenic.tsv.gz"

zgrep -v '^#' GRCh38.nr_deletions.pathogenic.tsv.gz | awk -F'\t' -v OFS='\t' '{split($12, arr, ";"); $13=arr[2]; $12=arr[1]; print $1, $2, $3, $10, $11, $12, $13}' > dbvar_del.bed
zgrep -v '^#' GRCh38.nr_duplications.pathogenic.tsv.gz | awk -F'\t' -v OFS='\t' '{split($12, arr, ";"); $13=arr[2]; $12=arr[1]; print $1, $2, $3, $10, $11, $12, $13}' > dbvar_dup.bed
zgrep -v '^#' GRCh38.nr_insertions.pathogenic.tsv.gz | awk -F'\t' -v OFS='\t' '{split($12, arr, ";"); $13=arr[2]; $12=arr[1]; print $1, $2, $3, $10, $11, $12, $13}' > dbvar_ins.bed
