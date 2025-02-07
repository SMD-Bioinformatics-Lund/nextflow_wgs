#!/usr/bin/env python
import argparse
import logging

FILTER_COL_IDX = 6
INFO_COL_IDX = 7

parser = argparse.ArgumentParser(
    """
    Removes the following records from VCF:
        1. FILTER not marked as . or PASS.
        2  INFO-field missing CSQ field.
    """
)
parser.add_argument("vep_annotated_vcf", type=str, help="Path to VEP-annotated VCF")
args = parser.parse_args()

logging.basicConfig(
    level=logging.INFO, format="[%(asctime)s][%(levelname)s]: %(message)s"
)


def csq_missing(info_field):
    """
    Variants not passing VEP filters are kept in VCF
    but not annotated with a CSQ field as of VEP 109.

    Prior to VEP 109 they were removed.
    """
    return "CSQ=" not in info_field


def filter_pass(filter_field):
    return filter_field.lower() in ["pass", "."]


nbr_filtered = 0

with open(args.vep_annotated_vcf) as vcf:
    for idx, line in enumerate(vcf):
        line = line.strip()

        if line[0] == "#":
            print(line)
            continue

        content = line.split("\t")
        curr_var = " ".join(content[0:2] + content[5:7])

        if not filter_pass(content[FILTER_COL_IDX]):
            logging.warning("%s: FILTER not marked PASS/. Skipping!", curr_var)
            nbr_filtered += 1
            continue

        if csq_missing(content[INFO_COL_IDX]):
            logging.warning("%s: CSQ field missing. Skipping!", curr_var)
            nbr_filtered += 1
            continue

        print(line)

logging.info("Removed %s of %s variants", nbr_filtered, idx)
