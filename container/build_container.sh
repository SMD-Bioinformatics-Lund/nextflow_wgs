#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
SINGULARITY_DEF="${SCRIPT_DIR}/Singularity"
OUTPUT_IMAGE="${1:-wgs_$(date +%Y-%m-%d).sif}"

if [[ $# -gt 1 ]]; then
    echo "Usage: $(basename "$0") [output_image.sif]" >&2
    exit 1
fi

if ! command -v singularity >/dev/null 2>&1; then
    echo "ERROR: 'singularity' command not found in PATH." >&2
    exit 1
fi

if [[ ! -f "${SINGULARITY_DEF}" ]]; then
    echo "ERROR: Missing Singularity definition file: ${SINGULARITY_DEF}" >&2
    exit 1
fi

sudo -E singularity build "${OUTPUT_IMAGE}" "${SINGULARITY_DEF}"
