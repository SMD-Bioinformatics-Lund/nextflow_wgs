DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && cd .. && pwd )"
cd "${DIR}"
echo "> Current working directory: ${DIR}"

DEST_HOST="rs-fs1.lunarc.lu.se"
PIPELINE_DEST="/fs1/pipelines/wgs_germline38"
DEST="${DEST_HOST}:${PIPELINE_DEST}"

current_branch=$(git branch | grep "^*" | sed "s/^* //")
if [[ "${current_branch}" != "master" ]]; then
    echo "> You are not on the master branch. Current branch: ${current_branch}"
    echo "> Do you want to deploy anyway? (y/n)"
    read -r deploy_despite_not_master
    if ! [[ "${deploy_despite_not_master}" =~ ^[yY]$ ]]; then
        echo "Aborting"
        exit 0
    fi
fi

echo "> Confirm, there are no running jobs for this pipeline? (y/n)"
read -r no_running

if ! [[ "${no_running}" =~ ^[yY]$ ]]; then
    echo "Aborting"
    exit 0
fi

echo "> Retrieving current git.hash from ${DEST}/git.hash, please wait ..."
remote_hash=$(ssh "${DEST_HOST}" "cat $PIPELINE_DEST/git.hash")

local_hash=$(git rev-parse HEAD)

echo "> Local hash is ${local_hash}, remote hash is ${remote_hash}"
echo "> Writing local hash to ${DIR}/git.hash"
echo "${local_hash}" > "${DIR}/git.hash"

echo "> Showing diff to ${remote_hash}"
git diff ${remote_hash} --stat

echo "> Do you want to view the full diff? (y/n)"
read -r view_full

if [[ "${view_full}" =~ ^[Yy]$ ]]; then
    git diff ${remote_hash}
fi

echo "> Do you want to proceed with deploying? (y/n)"
read -r response

if [[ "${response}" =~ ^[Yy]$ ]]; then

    echo "Deploying to ${DEST} ..."

    # Copy pipeline script
    scp "${DIR}/main.nf" "${DEST}"

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
