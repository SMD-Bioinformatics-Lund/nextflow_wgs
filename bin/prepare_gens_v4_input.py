#!/usr/bin/env python3

import argparse
from collections import defaultdict
import gzip
from pathlib import Path
from typing import Dict, List, Literal, TextIO
import logging

logging.basicConfig(level=logging.INFO, format="%(message)s")
LOG = logging.getLogger(__name__)

description = """
Generates inputs for Gens v4+

* Parses ROH output and UPD output into a single bed file
* Calculates overal ROH %
* Calculates average chromosome coverage
* Summarizes various UPD / ROH metrics
"""

ROH_COLOR = "rgb(255,186,60)"
UPD_COLOR = "rgb(255,75,75)"

AUTO_CHROMS = [str(i) for i in range(1, 23)]
SEX_CHROMS = ["X", "Y"]
CHROMS = AUTO_CHROMS + SEX_CHROMS


def main(
    roh_path: Path,
    upd_path: Path,
    upd_sites_path: Path,
    cov_path: Path,
    chrom_length_path: Path,
    sample: str,
    sex: str,
    roh_quality_threshold: int,
    cov_diff_threshold: int,
    output_upd_roh: Path,
    output_chr_cov: Path,
    output_meta: Path,
):
    LOG.info("Starting up")

    LOG.info("Parsing chromosomes")
    chrom_lengths: dict[str, int] = parse_chrom_lengths(chrom_length_path)
    total_chrom_length = sum(chrom_lengths.values())
    LOG.info("Parsing ROH")
    roh_entries: List[RohEntry] = parse_roh(roh_path, sample, roh_quality_threshold)
    LOG.info("Parsing UPD")
    upd_entries: List[UPDEntry] = parse_upd(upd_path)
    LOG.info("Parsing UPD sites")
    upd_site_info: Dict[str, dict[str, str]] = parse_upd_sites(upd_sites_path)
    LOG.info("Parsing coverage")
    avg_cov_entries: List[ChromCovEntry] = parse_cov(cov_path, cov_diff_threshold, sex)

    tot_roh_length = sum(
        [entry.get_length() for entry in roh_entries if entry.chrom in AUTO_CHROMS]
    )
    roh_perc = float(tot_roh_length) / float(total_chrom_length) * 100

    LOG.info("Writing global meta")
    with open_file(output_upd_roh, "w") as out_fh:
        for entry in roh_entries:
            print("\t".join(entry.get_bed_fields()), file=out_fh)
        for entry in upd_entries:
            print("\t".join(entry.get_bed_fields()), file=out_fh)

    LOG.info("Writing per-chromosome meta")
    with open_file(output_chr_cov, "w") as out_fh:

        print("\t".join(["Chromosome", "type", "value", "color"]), file=out_fh)
        for chrom_cov in avg_cov_entries:
            label = "Estimated chromosomal copy numbers"
            print(f"{chrom_cov.chrom}\t{label}\t{chrom_cov.cov}\t{chrom_cov.color}", file=out_fh)

        upd_labels = [
            # "Chromosome",
            "Total SNPs",
            "Non-informative",
            "Mismatch father",
            "Mismatch mother",
            "Anti-UPD",
        ]
        for chrom in CHROMS:
            for upd_label in upd_labels:
                field_value = upd_site_info[chrom][upd_label]
                print(f"{chrom}\t{upd_label}\t{field_value}\trgb(0,0,0)", file=out_fh)

    with open_file(output_meta, "w") as out_fh:
        print("\t".join(["type", "value"]), file=out_fh)
        print("\t".join(["%ROH", str(roh_perc)]), file=out_fh)


class RohEntry:
    def __init__(self, line):
        line = line.rstrip()
        fields = line.split("\t")
        self.sample = fields[1]
        self.chrom = fields[2]
        self.start = int(fields[3])
        self.end = int(fields[4])
        self.length = float(fields[5])
        self.nbr_markers = float(fields[6])
        self.quality = float(fields[7])

    def get_length(self) -> int:
        return self.end - self.start

    def get_bed_fields(self) -> List[str]:
        return [self.chrom, str(self.start), str(self.end), "ROH", ".", ".", ".", ".", ROH_COLOR]


class UPDEntry:
    def __init__(self, line):
        line = line.rstrip()
        fields = line.split("\t")
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        details_str = fields[3]

        details = {}
        for entry in details_str.split(";"):
            key, value = entry.split("=")
            details[key] = value
        self.details = details

    def get_length(self) -> int:
        return self.end - self.start

    def get_bed_fields(self) -> List[str]:
        return [self.chrom, str(self.start), str(self.end), "UPD", ".", ".", ".", ".", UPD_COLOR]


class CovEntry:
    def __init__(self, line):
        line = line.rstrip()
        fields = line.split("\t")
        self.chrom = fields[0]
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.cn = float(fields[3])

    def get_length(self) -> int:
        return self.end - self.start


class ChromCovEntry:
    def __init__(self, chrom: str, cov: float, color: str):
        self.chrom = chrom
        self.cov = cov
        self.color = color

    def get_fields(self) -> List[str]:
        return [self.chrom, str(self.cov), self.color]


def parse_upd_sites(upd_sites: Path) -> dict[str, dict[str, str]]:

    sum_per_chrom: dict[str, int] = {}
    sum_chrom_type: dict[str, dict[str, int]] = {}

    for chrom in CHROMS:
        sum_per_chrom[chrom] = 0
        sum_chrom_type[chrom] = {}

    with open_file(upd_sites, "r") as in_fh:
        for line in in_fh:
            line = line.rstrip()
            fields = line.split("\t")
            chrom = fields[0]
            upd_type = fields[3]

            sum_per_chrom[chrom] += 1

            if sum_chrom_type[chrom].get(upd_type) is None:
                sum_chrom_type[chrom][upd_type] = 0

            sum_chrom_type[chrom][upd_type] += 1

    out_fields: dict[str, dict[str, str]] = defaultdict(dict)
    for chrom in CHROMS:
        chrom_tot = sum_per_chrom[chrom]
        chrom_types: dict[str, int] = sum_chrom_type[chrom]

        paternal = chrom_types.get("UPD_PATERNAL_ORIGIN") or 0
        paternal_perc = round(100 * paternal / chrom_tot, 1) if chrom_tot > 0 else 0
        maternal = chrom_types.get("UPD_MATERNAL_ORIGIN") or 0
        maternal_perc = round(100 * maternal / chrom_tot, 1) if chrom_tot > 0 else 0
        anti = chrom_types.get("ANTI_UPD") or 0
        anti_perc = round(100 * anti / chrom_tot, 1) if chrom_tot > 0 else 0
        non_informative = paternal + maternal + anti
        non_informative_perc = round(100 * non_informative / chrom_tot, 1) if chrom_tot > 0 else 0

        chrom_info: dict[str, str] = {
            "Total SNPs": str(chrom_tot),
            "Non-informative": f"{non_informative} ({non_informative_perc}%)",
            "Mismatch father": f"{paternal} ({paternal_perc}%)",
            "Mismatch mother": f"{maternal} ({maternal_perc}%)",
            "Anti-UPD": f"{anti} ({anti_perc}%)",
        }

        out_fields[chrom] = chrom_info

    return out_fields


def parse_chrom_lengths(chrom_lengths_path: Path) -> dict[str, int]:
    chrom_lenghts = {}
    with open(chrom_lengths_path, "r") as in_fh:
        for line in in_fh:
            line = line.rstrip()
            if not line.startswith("@SQ"):
                continue
            fields = line.split("\t")
            chrom = fields[1].split(":")[1]
            length = fields[2].split(":")[1]
            if chrom in AUTO_CHROMS:
                chrom_lenghts[chrom] = int(length)
    return chrom_lenghts


def parse_cov(cov: Path, cov_diff_thres: float, sex: str) -> List[ChromCovEntry]:
    cov_sums = {}
    with open_file(cov, "r") as cov_fh:
        for line in cov_fh:
            if line.startswith("@") or line.startswith("CONTIG"):
                continue
            cov_entry = CovEntry(line)
            if cov_sums.get(cov_entry.chrom) is None:
                cov_sums[cov_entry.chrom] = {"nval": 0, "acc": 0, "sum_cov": 0, "sum_length": 0}

            cov_sums[cov_entry.chrom]["nval"] += 1
            cov_sums[cov_entry.chrom]["acc"] += cov_entry.cn

    avg_covs = {}
    for chrom, sums in cov_sums.items():
        avg_covs[chrom] = sums["acc"] / sums["nval"]

    scaled_covs = []
    for chrom in CHROMS:

        color = "rgb(0,0,0)"
        chrom_count = 1 if sex == "M" and chrom in SEX_CHROMS else 2
        # Calculate the full value
        cn_val = chrom_count * 2 ** avg_covs[chrom]
        if abs(cn_val - chrom_count) > cov_diff_thres:
            color = "rgb(255,0,0)"

        scaled_covs.append(ChromCovEntry(chrom, round(cn_val, 2), color))

    return scaled_covs


def parse_roh(roh_path: Path, sample: str, qual_thres: float) -> List[RohEntry]:

    roh_entries_per_chrom: List[RohEntry] = []
    with open_file(roh_path, "r") as roh_fh:
        for line in roh_fh:
            if line.startswith("#") or line.startswith("ST"):
                continue
            line = line.rstrip()
            fields = line.split("\t")
            # There is no quality information available
            if len(fields) < 7:
                continue
            roh_entry = RohEntry(line)
            if roh_entry.sample == sample and roh_entry.quality > qual_thres:
                roh_entries_per_chrom.append(roh_entry)
    return roh_entries_per_chrom


def parse_upd(upd: Path) -> List[UPDEntry]:

    upd_entries: List[UPDEntry] = []
    with open_file(upd, "r") as upd_fh:
        for line in upd_fh:
            upd_entry = UPDEntry(line)
            upd_entries.append(upd_entry)
    return upd_entries


def open_file(path: Path, read_or_write: Literal["r", "w"]) -> TextIO:
    if path.suffix == ".gz":
        if read_or_write == "r":
            return gzip.open(path, "rt")
        else:
            return gzip.open(path, "wt")
    else:
        return path.open(read_or_write)


def parse_arguments():
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--roh", required=True, type=Path, help="ROH output as provided by the bcftools command roh"
    )
    parser.add_argument(
        "--upd_regions",
        required=True,
        type=Path,
        help="UPD regions output from: https://github.com/bjhall/upd",
    )
    parser.add_argument(
        "--upd_sites",
        required=True,
        type=Path,
        help="UPD sites output from: https://github.com/bjhall/upd",
    )
    parser.add_argument(
        "--cov",
        required=True,
        type=Path,
        help="Copy ratio output from GATK's CollectReadCounts+DenoiseReadCounts",
    )
    parser.add_argument(
        "--chrom_lengths",
        required=True,
        type=Path,
        help=".dict file generated by calculating FASTA sequences lengths using samtools faidx and Picard's CreateSequenceDictionary",
    )

    parser.add_argument(
        "--sample",
        required=True,
        type=str,
        help="Sample ID used to filter out proband-only rows from ROH input",
    )
    parser.add_argument(
        "--sex",
        required=True,
        type=str,
        help="Sex information (M / F) used to scale copy ratio to chromosome level for sex chromosomes",
    )

    parser.add_argument(
        "--roh_quality_threshold",
        default=85,
        type=int,
        help="Only entries with quality threshold above this in --roh are considered",
    )
    parser.add_argument(
        "--cov_diff_threshold",
        default=0.1,
        type=int,
        help="Chromosome coverage values differing more than this are colored in red in output data",
    )

    parser.add_argument(
        "--out_gens_track", required=True, type=Path, help="Bed ranges with UPD / ROH information"
    )
    parser.add_argument(
        "--out_meta", required=True, type=Path, help="Output global meta info (i.e. ROH%)"
    )
    parser.add_argument(
        "--out_chrom_meta",
        required=True,
        type=Path,
        help="Writes per-chromosome information about coverage and UPD details",
    )

    args = parser.parse_args()

    if not args.sex in {"M", "F"}:
        raise ValueError("--sex must be either 'M' or 'F'")

    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.roh,
        args.upd,
        args.upd_sites,
        args.cov,
        args.chrom_lengths,
        args.sample,
        args.sex,
        args.roh_quality_threshold,
        args.cov_diff_threshold,
        args.output_upd_roh,
        args.output_chr_cov,
        args.output_meta,
    )
