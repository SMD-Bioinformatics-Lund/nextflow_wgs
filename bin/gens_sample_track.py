#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path
from typing import Literal, TextIO


description = """
Parse ROH output and UPD output into a single bed file
This bed will then be loaded into Gens as a sample annotation track
"""

# FIXME: Is this the UPD?
TOP_COLOR = "rgba(255,186,60)"
BOTTOM_COLOR = "rgba(255,75,75)"

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
SEX_CHROMS = ["X", "Y"]


def open_file(path: Path, read_or_write: Literal["r", "w"]) -> TextIO:
    if path.suffix == ".gz":
        if read_or_write == "r":
            return gzip.open(path, "rt")
        else:
            return gzip.open(path, "wt")
    else:
        return path.open(read_or_write)


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

    def get_bed_fields(self) -> list[str]:
        return [self.chrom, str(self.start), str(self.end), "ROH", ".", ".", ".", ".", TOP_COLOR]


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

    def get_bed_fields(self) -> list[str]:
        return [self.chrom, str(self.start), str(self.end), "UPD", ".", ".", ".", ".", BOTTOM_COLOR]


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

    def get_fields(self) -> list[str]:
        return [self.chrom, str(self.cov), self.color]


def main(
    roh_path: Path,
    upd_path: Path,
    cov_path: Path,
    chrom_length_path: Path,
    sample: str,
    sex: str,
    roh_quality_threshold: int,
    min_length: int,
    cov_diff_threshold: int,
    output_upd_roh: Path,
    output_chr_cov: Path,
    output_meta: Path,
):
    chrom_lengths: dict[str, int] = parse_chrom_lengths(chrom_length_path)
    total_chrom_length = sum(chrom_lengths.values())
    roh_entries: list[RohEntry] = parse_roh(roh_path, sample, roh_quality_threshold)
    upd_entries: list[UPDEntry] = parse_upd(upd_path)
    avg_cov_entries: list[ChromCovEntry] = parse_cov(cov_path, cov_diff_threshold, sex)

    print("Number ROH entries", len(roh_entries))

    tot_roh_length = sum(
        [
            entry.get_length()
            for entry in roh_entries
            if entry.chrom not in SEX_CHROMS and entry.chrom in CHROMS
        ]
    )
    print("ROH tot length", tot_roh_length)
    print("Tot chrom length", total_chrom_length)
    roh_perc = float(tot_roh_length) / float(total_chrom_length) * 100

    with open_file(output_upd_roh, "w") as out_fh:
        for entry in roh_entries:
            if entry.get_length() > min_length:
                print("\t".join(entry.get_bed_fields()), file=out_fh)
        for entry in upd_entries:
            if entry.get_length() > min_length:
                print("\t".join(entry.get_bed_fields()), file=out_fh)

    with open_file(output_chr_cov, "w") as out_fh:
        print("\t".join(["row_name", "type", "value", "color"]), file=out_fh)
        for chrom_cov in avg_cov_entries:
            label = "Estimated chromosomal copy numbers"
            print(f"{chrom_cov.chrom}\t{label}\t{chrom_cov.cov}\t{chrom_cov.color}", file=out_fh)

    with open_file(output_meta, "w") as out_fh:
        print("\t".join(["type", "value"]), file=out_fh)
        print(roh_perc)


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
            if chrom not in SEX_CHROMS and chrom in CHROMS:
                chrom_lenghts[chrom] = int(length)
    return chrom_lenghts


def parse_cov(cov: Path, cov_diff_thres: float, sex: str) -> list[ChromCovEntry]:
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

    sex_chroms = ["X", "Y"]

    scaled_covs = []
    for chrom in CHROMS:

        color = "rgb(0,0,0)"
        chrom_count = 1 if sex == "M" and chrom in sex_chroms else 2
        # Calculate the full value
        cn_val = chrom_count * 2 ** avg_covs[chrom]
        if abs(cn_val - chrom_count) > cov_diff_thres:
            color = "rgb(255,0,0)"

        scaled_covs.append(ChromCovEntry(chrom, cn_val, color))

    return scaled_covs


def parse_roh(roh_path: Path, sample: str, qual_thres: float) -> list[RohEntry]:

    roh_entries_per_chrom: list[RohEntry] = []
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


def parse_upd(upd: Path) -> list[UPDEntry]:

    upd_entries: list[UPDEntry] = []
    with open_file(upd, "r") as upd_fh:
        for line in upd_fh:
            upd_entry = UPDEntry(line)
            upd_entries.append(upd_entry)
    return upd_entries


def parse_arguments():
    parser = argparse.ArgumentParser()

    # FIXME: Understand the difference roh / upd

    parser.add_argument("--roh", required=True, type=Path)
    parser.add_argument("--upd", required=True, type=Path)
    parser.add_argument("--cov", required=True, type=Path)
    parser.add_argument("--chrom_lengths", required=True, type=Path)

    parser.add_argument("--sample", required=True, type=str)
    parser.add_argument("--sex", required=True, type=str)

    parser.add_argument("--roh_quality_threshold", default=85, type=int)
    parser.add_argument("--min_length", default=100000, type=int)
    parser.add_argument("--cov_diff_threshold", default=10, type=int)

    parser.add_argument("--output_upd_roh", type=Path)
    parser.add_argument("--output_chr_cov", type=Path)
    parser.add_argument("--output_meta", type=Path)

    args = parser.parse_args()

    if not args.sex in {"M", "F"}:
        raise ValueError("--sex must be either 'M' or 'F'")

    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.roh,
        args.upd,
        args.cov,
        args.chrom_lengths,
        args.sample,
        args.sex,
        args.roh_quality_threshold,
        args.min_length,
        args.cov_diff_threshold,
        args.output_upd_roh,
        args.output_chr_cov,
        args.output_meta,
    )
