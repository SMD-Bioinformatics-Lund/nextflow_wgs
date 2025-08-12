#!/usr/bin/env python3

import json
import sys

# Contact:     Alexander Koc <alexander.koc@skane.se>
# Date:        2023-06-21
# Description: Merge JSON files into one.


def merge_json_files(file_paths):
    merged_data = {}

    for file_path in file_paths:
        with open(file_path, "r") as file:
            try:
                json_data = json.load(file)
                merged_data.update(json_data)
            except json.JSONDecodeError:
                print(f"Error: Invalid JSON file '{file_path}'", file=sys.stderr)
                continue

    return merged_data


if __name__ == "__main__":
    if len(sys.argv) < 1:
        print("Usage: python3 merge_json_files.py <file1.json> <file2.json> ...")
        sys.exit(1)

    file_paths = sys.argv[1:]
    merged_json = merge_json_files(file_paths)
    print(json.dumps(merged_json, indent=4, sort_keys=True))
