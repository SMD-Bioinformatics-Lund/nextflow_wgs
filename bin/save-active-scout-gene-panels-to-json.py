#!/usr/bin/env python3

"""
Generate a dump of gene panels from scout containing only the latest
versions of active/unarchived panels. Used to generate gene panel input
to create_yml.pl in the nextflow_wgs pipeline.

If the script breaks for any reason, then a JSON dump of scouts gene_panel
collection, generated using the python lib bson.json_util will work.

Run on lennart.
"""

import logging
import os
import subprocess
import sys
import tempfile
from datetime import date

import bson.json_util
import pymongo

LOG_FORMAT = "[%(asctime)s][%(levelname)s]: %(message)s"

logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
LOG = logging.getLogger(__name__)

SCOUT_MONGO_HOST = "mtcmdpgm01.lund.skane.se:27017"
REMOTE_HOST = "rs-fe1"
REMOTE_OUTPUT_DIR = "/fs1/resources/scout/active_gene_panels/"
REMOTE_LATEST_SYMLINK_PATH = "/fs1/resources/scout/scout_active_gene_panels_LATEST.json"


panels_db: pymongo.collection.Collection = pymongo.MongoClient(
    SCOUT_MONGO_HOST
).scout.gene_panel

# Fetch only active/unarchived panels, fetch only the latest versions, keep only name, version and institute
ONLY_ACTIVE_PANELS_PIPELINE = [
    # Step 1: Sort by 'panel_name' and 'version' (descending order to get the latest version first)
    {
        "$sort": {
            "panel_name": 1,  # Ascending by 'panel_name'
            "version": -1,  # Descending by 'version'
        }
    },
    # Step 2: Group by 'panel_name' to get the latest version of each panel
    {
        "$group": {
            "_id": "$panel_name",  # Group by 'panel_name'
            "latest_version": {
                "$first": "$version"
            },  # Take the first version (the latest)
            "institute": {"$first": "$institute"},  # Keep the corresponding institute
            "id": {"$first": "$_id"},  # Keep the ID of the latest version document
            "is_archived": {
                "$first": "$is_archived"
            },  # Capture the 'is_archived' field
            "hidden": {"$first": "$hidden"},  # Capture the 'hidden' field
        }
    },
    # Step 3: Filter out hidden or archived panels
    {
        "$match": {
            "$and": [
                {"$or": [{"hidden": None}, {"hidden": False}]},  # Not hidden
                {
                    "$or": [{"is_archived": None}, {"is_archived": False}]
                },  # Not archived
            ]
        }
    },
    # Step 4: Project the required fields
    {
        "$project": {
            "_id": "$id",  # Return the original _id of the latest version document
            "panel_name": "$_id",  # Rename _id to panel_name
            "version": "$latest_version",  # Keep the latest version
            "institute": 1,  # Keep the institute as it is
        }
    },
]


def main() -> None:
    """
    Generate gene panel dump, copy to remote and symlink copied file to
    REMOTE_LATEST_SYMLINK_PATH
    """
    LOG.info("Hello.")

    LOG.info("Fetching panel dump.")
    result: pymongo.command_cursor.CommandCursor = panels_db.aggregate(
        ONLY_ACTIVE_PANELS_PIPELINE
    )

    with tempfile.NamedTemporaryFile(mode="w", delete=True) as temp_file:

        LOG.debug("Dumping aggregated result to %s", temp_file.name)
        temp_file.write(bson.json_util.dumps(result))

        temp_file.flush()
        today = date.today().isoformat()
        remote_filepath = f"{REMOTE_OUTPUT_DIR}scout_gene_panels_active_{today}.json"
        LOG.debug("final remote_filepath: %s", remote_filepath)

        copy_panel_dump_to_remote(temp_file.name, remote_filepath)
        set_read_permissions(remote_filepath)

    symlink_to_latest(remote_filepath)
    LOG.info("All done.")
    LOG.info("Bye.")


def set_read_permissions(remote_filepath: str) -> None:
    """
    Make JSON readable by everyone. Otherwise the JSON is copied over
    with the restricted permissions of the NamedTemporaryFile
    """
    set_read_permissions_everyone = [
        "ssh",
        REMOTE_HOST,
        "chmod a+r",
        remote_filepath,
        "&&" "chmod a-w",
        remote_filepath,
    ]

    try:
        LOG.info("Setting permissions a+r, a-r for %s", remote_filepath)
        subprocess.run(set_read_permissions_everyone, check=True)
    except subprocess.CalledProcessError as e:
        LOG.error("Error occurred when setting permissions: %s", e)
        sys.exit(1)


def copy_panel_dump_to_remote(temp_filepath: str, remote_filepath: str) -> None:
    """
    Copy panel dump to the cluster with scp
    """

    LOG.debug(
        "Source file exists at %s: %s", temp_filepath, os.path.exists(temp_filepath)
    )
    LOG.info("Copying %s to %s:%s", temp_filepath, REMOTE_HOST, remote_filepath)

    scp_command = ["scp", temp_filepath, f"{REMOTE_HOST}:{remote_filepath}"]

    try:
        subprocess.run(scp_command, check=True)
        LOG.info("SCP transfer completed!")
        LOG.debug(
            "Remote file exists at: %s:%s: %s",
            REMOTE_HOST,
            remote_filepath,
            os.path.exists(remote_filepath),
        )
    except subprocess.CalledProcessError as e:
        LOG.error("Error occurred during SCP transfer: %s", e)


def backup_existing_symlink() -> None:
    """
    Back up the symlink that points to the latest panel dump
    before deleting and generating new
    """
    backup_ssh_command = [
        "ssh",
        REMOTE_HOST,
        "cp",
        "-P",
        REMOTE_LATEST_SYMLINK_PATH,
        f"{REMOTE_LATEST_SYMLINK_PATH}.bkp",
    ]

    try:
        subprocess.run(backup_ssh_command, check=True)
        LOG.info(
            "Successfully backed up old symlink to %s.bkp", REMOTE_LATEST_SYMLINK_PATH
        )
    except subprocess.CalledProcessError as e:
        LOG.error("Error occurred during symlink bkp transfer: %s", e)
        sys.exit(1)


def restore_symlink() -> None:
    """
    Restore old symlink in case of failure to generate new.
    """
    restore_symlink_command = [
        "ssh",
        REMOTE_HOST,
        "mv",
        f"{REMOTE_LATEST_SYMLINK_PATH}.bkp",
        REMOTE_LATEST_SYMLINK_PATH,
    ]

    try:
        subprocess.run(restore_symlink_command, check=True)
        LOG.info("Restored old symlink from %s.bkp", REMOTE_LATEST_SYMLINK_PATH)
    except subprocess.CalledProcessError as e:
        LOG.critical("Error occurred during symlink restore: %s", e)
        sys.exit(1)


def remove_existing_symlink() -> None:
    """
    Remove existing symlink pointing to LATEST panel dump on remote
    """
    ssh_command = ["ssh", REMOTE_HOST, "rm", REMOTE_LATEST_SYMLINK_PATH]

    try:
        subprocess.run(ssh_command, check=True)
        LOG.info("Removed current old symlink: %s", REMOTE_LATEST_SYMLINK_PATH)
    except subprocess.CalledProcessError as e:
        LOG.error("Error occurred during symlink bkp transfer: %s", e)
        sys.exit(1)


def symlink_to_latest(remote_filepath) -> None:
    """
    Generate symlink that points to the new panel dump
    """
    LOG.info("Symlinking %s", remote_filepath)

    # symlink_exists
    if os.path.exists(REMOTE_LATEST_SYMLINK_PATH):
        LOG.info(
            "Symlink already exists at %s. Attempting backup and remove.",
            REMOTE_LATEST_SYMLINK_PATH,
        )
        backup_existing_symlink()
        remove_existing_symlink()

    new_symlink_command = [
        "ssh",
        REMOTE_HOST,
        "ln",
        "-s",
        remote_filepath,
        REMOTE_LATEST_SYMLINK_PATH,
    ]

    try:
        subprocess.run(new_symlink_command, check=True)
        LOG.info(
            "%s:%s symlinked to %s:%s",
            REMOTE_HOST,
            remote_filepath,
            REMOTE_HOST,
            REMOTE_LATEST_SYMLINK_PATH,
        )
    except subprocess.CalledProcessError as e:
        LOG.error(
            "Error occurred during creation of new symlink bkp: %s. Restoring old symlink!",
            e,
        )
        restore_symlink()


if __name__ == "__main__":
    main()
