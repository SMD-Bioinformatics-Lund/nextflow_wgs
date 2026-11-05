#!/usr/bin/env python3
"""Reject reference records that SVDB would skip or misread during annotation."""

import argparse
import gzip
import math


def validate_reference(path):
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    records = 0
    with opener(path, "rt") as source:
        for line_number, line in enumerate(source, 1):
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            location = f"{path}:{line_number}"
            if len(columns) < 8:
                raise ValueError(f"{location}: expected at least eight VCF columns")
            info = dict(field.split("=", 1) for field in columns[7].split(";") if "=" in field)
            try:
                count = int(info["N_HET"])
                frequency = float(info["AF"])
            except (KeyError, ValueError) as error:
                raise ValueError(
                    f"{location} ({columns[0]}:{columns[1]}): gnomAD SV annotation requires "
                    "scalar numeric N_HET and AF in every reference record; AC is not required"
                ) from error
            if count < 0 or not math.isfinite(frequency) or not 0 <= frequency <= 1:
                raise ValueError(f"{location}: expected N_HET >= 0 and finite AF between 0 and 1")
            records += 1
    if not records:
        raise ValueError(f"{path}: gnomAD SV reference contains no variant records")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", help="gnomAD slim sites VCF, optionally gzip compressed")
    args = parser.parse_args()
    try:
        validate_reference(args.reference)
    except (ValueError, OSError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
