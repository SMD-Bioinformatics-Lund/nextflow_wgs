#!/usr/bin/env python3

"""
This script ensures that the wgs trio-analysis SV `svdb_origin`
and `set` INFO values can be parsed by scout.

It converts the default nf_wgs output
`svdb_origin=manta1|manta2|tiddit1` -> `svdb_origin=manta|tiddit`
for both the set and svdb_origin fields.
"""

import argparse
import logging
import sys
import gzip
from collections import OrderedDict
from typing import Optional, List, Dict, Union

SVDB_ORIGIN_SEPARATOR = "|"
SVDB_ORIGIN_KEY = "svdb_origin"
SVDB_SET_KEY = "set"
SVDB_SET_SEPARATOR = "-"

SVDB_TARGET_FIELDS = [
    {"field_name": SVDB_SET_KEY, "separator": SVDB_SET_SEPARATOR},
    {"field_name": SVDB_ORIGIN_KEY, "separator": SVDB_ORIGIN_SEPARATOR},
]


def main() -> None:
    args = parse_args()

    vcf_infile = args.vcf_file
    callers = args.callers

    # Open the VCF file (supports gzip)
    if vcf_infile.endswith(".gz"):
        f = gzip.open(vcf_infile, "rt")
    else:
        f = open(vcf_infile, "r")

    for row_idx, line in enumerate(f):
        curr_record = str(line)

        # Retain header curr_records
        if curr_record.startswith("#"):
            print(curr_record, end="")
            continue

        fields = curr_record.rstrip("\n").split("\t")
        chrom, pos, id_, ref, alt, qual, filter_, info = fields[:8]
        others = fields[8:]  # If there are additional columns

        info_dict = parse_info_field(info_field=info)
        try:
            updated_info_dict = clean_up_caller_names_in_svdb_info_fields(
                info_dict, callers
            )
        except ValueError as e:
            raise ValueError(f"Line {row_idx}: {e}")

        # reconstruct the INFO field while retaining the order
        new_info = []
        for key, value in updated_info_dict.items():
            if value is True:
                new_info.append(f"{key}")
            else:
                new_info.append(f"{key}={value}")
        info_str = ";".join(new_info)

        # Reconstruct and print the modified VCF line
        new_fields = fields[:7] + [info_str] + others
        print("\t".join(new_fields))

    f.close()


def clean_up_caller_names_in_svdb_info_fields(
    info_row: Dict[str, Union[str, bool]], used_callers: List[str]
) -> Dict[str, Union[str, bool]]:
    """
    Iterate over all target SVDB fields in curr info row
    and clean up target caller names inside so they can
    be parsed by scout.
    """

    parsed_info_row = info_row.copy()

    for svdb_field_info in SVDB_TARGET_FIELDS:
        field_name: str = svdb_field_info["field_name"]
        separator: str = svdb_field_info["separator"]

        if field_name not in parsed_info_row:
            logging.warning(
                "%s not found in INFO.",
                field_name,
            )
            continue

        target_svdb_info_value = parsed_info_row[field_name]

        if not isinstance(target_svdb_info_value, str):
            raise ValueError(
                f"'{field_name}' parsed as {type(parsed_info_row[field_name])}, must be str."
            )

        parsed_info_row[field_name] = reduce_to_set_of_unique_callers(
            target_svdb_info_value, separator=separator, callers=used_callers
        )

    return parsed_info_row


def parse_info_field(info_field: str) -> Dict[str, Union[str, bool]]:
    """
    Split info string and parse into dict whilst retaining order
    """
    info_entries = info_field.split(";")
    info_dict: Dict[str, Union[str, bool]] = OrderedDict()
    for entry in info_entries:
        if "=" in entry:
            key, value = entry.split("=", 1)
            info_dict[key] = value
        else:
            info_dict[entry] = True  # Flag without a value

    return info_dict


def match_svdb_annotation_with_caller(
    unparsed_svdb_origin_annotation: str, callers: List[str]
) -> Optional[str]:
    """
    Convert caller annotations manta1|gatk1|gatk2 -> manta|gatk

    Returns:
     - caller_name (str) i
     - None if caller not supported
    """
    for caller in callers:
        if caller in unparsed_svdb_origin_annotation:
            return caller

    return None


def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(description="Modify SVDB merged VCF script.")
    parser.add_argument(
        "vcf_file", type=str, help="Input VCF file (can be .vcf or .vcf.gz)."
    )
    parser.add_argument(
        "--callers", nargs="+", type=str, required=True, help="List of SV caller names."
    )

    args = parser.parse_args()

    if not args.vcf_file.endswith((".vcf", ".vcf.gz")):
        print(
            "Error: Input file must be a VCF file with .vcf or .vcf.gz extension.",
            file=sys.stderr,
        )
        sys.exit(1)

    return args


def reduce_to_set_of_unique_callers(
    svdb_annotation: str, separator: str, callers: List[str]
) -> str:
    """
    Fetches callers from specified svdb INFO annotations
    and reduces annotation to list of detected unique callers

    Returns a string in the same format as the input annotation
    """
    callers_in_svdb_annotation = svdb_annotation.split(separator)
    parsed_callers = set()

    for caller in callers_in_svdb_annotation:
        matched_caller = match_svdb_annotation_with_caller(caller, callers)

        # Add caller as-is if not supported
        if matched_caller is None:
            parsed_callers.add(caller)
            continue

        parsed_callers.add(matched_caller)

    # Sort the diffs to avoid random/order diff between runs.
    sorted_parsed_callers = sorted(parsed_callers)
    return separator.join(sorted_parsed_callers)


if __name__ == "__main__":
    main()
