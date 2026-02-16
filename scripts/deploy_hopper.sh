#!/usr/bin/env bash
set -euo pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && cd .. && pwd )"
cd "${DIR}"
echo "> Current working directory: ${DIR}"

DEST_HOST="rs-fs1.lunarc.lu.se"
PIPELINE_DEST="/fs1/pipelines/wgs_germline38"
DEST="${DEST_HOST}:${PIPELINE_DEST}"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "ERROR: ${DIR} is not inside a Git repository." >&2
    exit 1
fi

prompt_yes_no() {
    local prompt="$1"
    local response
    echo "${prompt} (y/n)"
    read -r response
    [[ "${response}" =~ ^[Yy]$ ]]
}

current_branch=$(git rev-parse --abbrev-ref HEAD)
if [[ "${current_branch}" != "master" ]]; then
    echo "> You are not on the master branch. Current branch: ${current_branch}"
    if ! prompt_yes_no "> Do you want to deploy anyway?"; then
        echo "Aborting"
        exit 0
    fi
fi

if ! prompt_yes_no "> Confirm, there are no running jobs for this pipeline?"; then
    echo "Aborting"
    exit 0
fi

echo "> Retrieving current git.hash from ${DEST}/git.hash, please wait ..."
if ! remote_hash=$(ssh "${DEST_HOST}" "if [ -f '${PIPELINE_DEST}/git.hash' ]; then cat '${PIPELINE_DEST}/git.hash'; fi"); then
    echo "ERROR: Failed to read remote git hash from ${DEST}/git.hash." >&2
    exit 1
fi
if [[ -z "${remote_hash}" ]]; then
    echo "ERROR: Remote git hash file is missing or empty: ${DEST}/git.hash" >&2
    exit 1
fi

local_hash=$(git rev-parse HEAD)

echo "> Local hash is ${local_hash}, remote hash is ${remote_hash}"
echo "> Writing local hash to ${DIR}/git.hash"
echo "${local_hash}" > "${DIR}/git.hash"

echo "> Showing diff to ${remote_hash}"
git diff "${remote_hash}" --stat

if prompt_yes_no "> Do you want to view the full diff?"; then
    git diff "${remote_hash}"
fi

if prompt_yes_no "> Do you want to proceed with deploying?"; then

    echo "Deploying to ${DEST} ..."

    required_paths=(
        "${DIR}/main.nf"
        "${DIR}/workflows"
        "${DIR}/modules"
        "${DIR}/bin"
        "${DIR}/configs/nextflow.hopper.config"
        "${DIR}/git.hash"
    )
    for required_path in "${required_paths[@]}"; do
        if [[ ! -e "${required_path}" ]]; then
            echo "ERROR: Required path missing: ${required_path}" >&2
            exit 1
        fi
    done

    # Copy pipeline script
    scp "${DIR}/main.nf" "${DEST}"

    # Copy workflows
    scp -r "${DIR}/workflows" "${DEST}"

    # Copy modules
    scp -r "${DIR}/modules" "${DEST}"

    # Copy configuration file
    scp "${DIR}/configs/nextflow.hopper.config" "${DEST}/nextflow.config"

    # Copy other files
    scp -r "${DIR}/bin" "${DEST}"

    #git rev-parse HEAD > git.hash
    scp "${DIR}/git.hash" "${DEST}"

else
    echo "Deploy aborted"
    exit 0
fi
