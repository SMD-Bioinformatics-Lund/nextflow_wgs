#!/usr/bin/env python3

import argparse
from collections import defaultdict


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Copy caller names from SVDB INFO.set to INFO.SCOUT_CUSTOM in SVDB merged VCFs"
    )
    parser.add_argument(
        "--merged_vcf", required=True, help="Path to the merged VCF file."
    )
    parser.add_argument(
        "--callers", required=True, help="Comma-separated list of used callers."
    )
    return parser.parse_args()


def process_vcf(merged_vcf, used_callers):
    """Process the VCF file and copy caller information to SCOUT_CUSTOM variant by variant."""
    used_callers_set = set(used_callers.split(","))

    with open(merged_vcf, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                print(line.strip())
                continue

            columns = line.strip().split("\t")
            chrom, pos, id_, ref, alt, qual, filter_, info, fmt = columns[:9]
            samples = columns[9:]

            info_dict = parse_info_field(info)

            callers_flag = defaultdict(
                bool, {"manta": False, "cnvkit": False, "gatk": False}
            )

            # Parse caller information from the "set" field in INFO
            set_field = info_dict.get("set", "")
            callers_in_set = set(set_field.split("-"))
            for caller in callers_in_set:
                if caller in callers_flag:
                    callers_flag[caller] = True
                elif "Intersection" in caller:
                    for used_caller in used_callers_set:
                        if used_caller in callers_flag:
                            callers_flag[used_caller] = True

            found_callers = [caller for caller, flag in callers_flag.items() if flag]

            # Do not touch SCOUT_CUSTOM if no callers found in set
            if not any(found_callers):
                continue

            new_scout_custom_values = {"Caller": "&".join(found_callers)}

            if "SCOUT_CUSTOM" in info_dict:
                existing_scout_custom = parse_scout_custom(info_dict["SCOUT_CUSTOM"])
                existing_scout_custom.update(new_scout_custom_values)
                info_dict["SCOUT_CUSTOM"] = build_scout_custom(existing_scout_custom)
            else:
                info_dict["SCOUT_CUSTOM"] = build_scout_custom(new_scout_custom_values)

            new_info = reconstruct_info_field(info_dict)
            print("\t".join(columns[:7] + [new_info, fmt] + samples))


def parse_info_field(info_field):
    """Parse the INFO field into a dictionary, handling key-value pairs and flags."""
    info_dict = {}
    for item in info_field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info_dict[key] = value
        else:
            # Treat key-only items as flags with a value of True
            info_dict[item] = True
    return info_dict


def reconstruct_info_field(info_dict):
    """Reconstruct the INFO field from the dictionary."""
    info_items = []
    for key, value in info_dict.items():
        if value is True:
            # Handle flags (key-only)
            info_items.append(key)
        else:
            # Handle key-value pairs
            info_items.append(f"{key}={value}")
    return ";".join(info_items)


def parse_scout_custom(scout_custom):
    """Parse the SCOUT_CUSTOM field into a dictionary."""
    custom_dict = {}
    for item in scout_custom.split(","):
        key, value = item.split("|", 1)
        custom_dict[key] = value
    return custom_dict


def build_scout_custom(custom_dict):
    """Rebuild the SCOUT_CUSTOM field from a dictionary."""
    return ",".join(f"{key}|{value}" for key, value in custom_dict.items())


if __name__ == "__main__":
    args = parse_arguments()
    process_vcf(args.merged_vcf, args.callers)