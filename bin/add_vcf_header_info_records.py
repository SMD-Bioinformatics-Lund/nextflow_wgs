#!/usr/bin/env python3

"""
Add INFO fields to VCF Header
"""

import argparse
import os
from typing import List, Dict, Optional

VALID_NUMBER_SPECIAL_CHARS = ("A", "R", "G", ".")
VALID_INFO_TYPES = ("Integer", "Float", "Flag", "Character", "String")


def main() -> None:
    """
    Add new INFO rows to VCF, output to new VCF file.
    """
    args = parse_arguments()

    if not os.path.exists(args.vcf):
        raise FileNotFoundError(f"Input VCF file '{args.vcf}' does not exist.")

    processed_header_data = process_info_args(args.info)
    info_lines_to_be_added = generate_info_lines(processed_header_data)
    add_info_to_vcf(args.vcf, args.output, info_lines_to_be_added)


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Add INFO fields to a VCF header.")

    parser.add_argument(
        "--vcf",
        required=True,
        help="Path to the input VCF file.",
    )

    parser.add_argument(
        "--info",
        action="append",
        nargs=6,  # Expect exactly six arguments
        metavar=("ID", "Number", "Type", "Description", "Source", "Version"),
        required=True,
        help=(
            "Define an INFO field in the header. Provide exactly six arguments: "
            "ID, Number, Type, Description, Source, and Version. Use an empty string "
            "(e.g., '') for Source and Version if not applicable. This argument can "
            "be repeated to add multiple INFO fields."
        ),
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output VCF file.",
    )

    return parser.parse_args()


def process_info_args(info_args: List[List[str]]) -> List[Dict[str, Optional[str]]]:
    """
    Process INFO arguments, ensuring valid Source and Version usage.

    Args:
        info_args (List[List[str]]): Parsed INFO arguments.

    Returns:
        List[Dict]: INFO fields as dictionaries.
    """

    id_: str
    number: str
    type_: str
    description: str
    source: Optional[str]
    version: Optional[str]

    processed_info = []
    for info in info_args:
        id_, number, type_, description, source, version = info

        if not source:
            source = None

        if not version:
            version = None

        validate_info_fields(number=number, info_type=type_)

        processed_info.append(
            {
                "ID": id_,
                "Number": number,
                "Type": type_,
                "Description": description,
                "Source": source,
                "Version": version,
            }
        )
    return processed_info


def validate_info_fields(*, number: str, info_type: str) -> None:
    if not number.isdigit() and number not in VALID_NUMBER_SPECIAL_CHARS:
        raise ValueError(
            f"Invalid number: '{number}'. Must be an integer or any of the following: {VALID_NUMBER_SPECIAL_CHARS}",
        )

    if info_type not in VALID_INFO_TYPES:
        raise ValueError(
            f"Invalid type: '{info_type}'. Must be one of the following: {VALID_INFO_TYPES}",
        )


def generate_info_lines(info_fields: List[Dict[str, Optional[str]]]) -> List[str]:
    """
    process info field data to header text rows
    """
    info_lines = []
    for info in info_fields:
        info_line = f'##INFO=<ID={info["ID"]},Number={info["Number"]},Type={info["Type"]},Description="{info["Description"]}"'

        if info.get("Source"):
            info_line += f',Source="{info["Source"]}"'

        if info.get("Version"):
            info_line += f',Version="{info["Version"]}"'

        info_line += ">"

        info_lines.append(info_line)
    return info_lines


def add_info_to_vcf(input_vcf: str, output_vcf: str, info_lines: List[str]) -> None:
    """
    add header data before VCF column header
    """

    found_column_header = False

    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#CHROM"):
                found_column_header = True
                # Insert INFO lines just before the column header line
                for info_line in info_lines:
                    outfile.write(info_line + "\n")
            outfile.write(line)

        if not found_column_header:
            raise ValueError(
                "Could not find column header (#CHROM and friends). Is this even a VCF file?"
            )


if __name__ == "__main__":
    main()
