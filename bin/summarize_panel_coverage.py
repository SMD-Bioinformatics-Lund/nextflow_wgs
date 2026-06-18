#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path
from typing import Any


DEFAULT_THRESHOLD = 500.0


def read_gene_list(gene_file: Path) -> list[str]:
    genes: list[str] = []
    seen: set[str] = set()

    with gene_file.open("r", encoding="utf-8") as infile:
        for line in infile:
            gene = line.strip()
            if not gene or gene.startswith("#"):
                continue
            if gene not in seen:
                genes.append(gene)
                seen.add(gene)

    return genes


def region_length(region: dict[str, Any]) -> int:
    start = int(region["start"])
    end = int(region["end"])
    return max(end - start, 0)


def region_coverage(region: dict[str, Any]) -> float:
    if "cov" not in region:
        return 0.0
    return float(region["cov"])


def weighted_mean_cds_coverage(cds_regions: dict[str, dict[str, Any]]) -> float | None:
    total_bases = 0
    coverage_sum = 0.0

    for region in cds_regions.values():
        length = region_length(region)
        total_bases += length
        coverage_sum += region_coverage(region) * length

    if total_bases == 0:
        return None

    return coverage_sum / total_bases


def summarize_coverage(
    coverage_data: dict[str, Any], genes_to_include: list[str], threshold: float
) -> dict[str, Any]:
    genes = coverage_data.get("genes", {})

    gene_results = []
    missing_genes = []
    genes_below_threshold = []

    for gene in genes_to_include:
        gene_record = genes.get(gene)
        if gene_record is None:
            missing_genes.append(gene)
            genes_below_threshold.append(gene)
            gene_results.append(
                {
                    "gene": gene,
                    "found": False,
                    "mean_cds_coverage": None,
                    "covered_by_mean_cds_coverage_threshold": False,
                }
            )
            continue

        cds_regions = gene_record.get("CDS", {})
        mean_coverage = weighted_mean_cds_coverage(cds_regions)
        is_covered = mean_coverage is not None and mean_coverage >= threshold

        if not is_covered:
            genes_below_threshold.append(gene)

        gene_results.append(
            {
                "gene": gene,
                "found": True,
                "mean_cds_coverage": round(mean_coverage, 2) if mean_coverage is not None else None,
                "covered_by_mean_cds_coverage_threshold": is_covered,
                "cds_regions": len(cds_regions),
                "cds_regions_without_coverage": sum(
                    1 for region in cds_regions.values() if "cov" not in region
                ),
                "cds_regions_below_threshold": [
                    region_name
                    for region_name, region in cds_regions.items()
                    if region_coverage(region) < threshold
                ],
            }
        )

    total_genes = len(genes_to_include)
    covered_genes = total_genes - len(genes_below_threshold)
    percent_covered = (covered_genes / total_genes * 100) if total_genes else 0.0

    return {
        "threshold": threshold,
        "genes_requested": total_genes,
        "genes_found": total_genes - len(missing_genes),
        "genes_covered_by_mean_cds_coverage_threshold": covered_genes,
        "percent_genes_covered_by_at_least_threshold_mean_cds_coverage": round(
            percent_covered, 2
        ),
        "genes_not_covered_by_at_least_threshold_mean_cds_coverage": genes_below_threshold,
        "missing_genes": missing_genes,
        "gene_results": gene_results,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize panel_coverage.py JSON for a selected gene list. "
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
