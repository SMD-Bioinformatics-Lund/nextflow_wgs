#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path

try:
    from bin.panel_coverage import read_gene_list, summarize_coverage
except ModuleNotFoundError:
    from panel_coverage import read_gene_list, summarize_coverage


DEFAULT_THRESHOLD = 500.0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compatibility wrapper for panel_coverage.py summary output. "
            "A gene is counted as covered when its length-weighted mean CDS coverage "
            "is at least the selected threshold."
        )
    )
    parser.add_argument("coverage_json", type=Path, help="panel_coverage.py JSON output")
    parser.add_argument(
        "-g",
        "--genes",
        required=True,
        type=Path,
        help="Text file with one gene symbol per line",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        default=DEFAULT_THRESHOLD,
        type=float,
        help=f"Mean CDS coverage threshold. Default: {DEFAULT_THRESHOLD:g}",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write summary JSON to this file instead of stdout",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    with args.coverage_json.open("r", encoding="utf-8") as infile:
        coverage_data = json.load(infile)

    genes_to_include = read_gene_list(args.genes)
    summary = summarize_coverage(coverage_data, genes_to_include, args.threshold)

    if args.output:
        args.output.write_text(json.dumps(summary, indent=4) + "\n", encoding="utf-8")
    else:
        json.dump(summary, sys.stdout, indent=4)
        print()


if __name__ == "__main__":
    main()
