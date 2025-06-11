#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path
from typing import TextIO


description = """
Parse ROH output and UPD output into a single bed file
This bed will then be loaded into Gens as a sample annotation track
"""


def open_file(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    else:
        return path.open("r")


# # This file was produced by: bcftools roh(1.9+htslib-1.9)
# # The command line was: bcftools roh --rec-rate 1e-9 --AF-tag GNOMADAF -o roh.txt giab-trio-viktortestrun_dev-vh-split_prescore_annotsvbump-bam.SNPs.vcf
# #
# # RG    [2]Sample       [3]Chromosome   [4]Start        [5]End  [6]Length (bp)  [7]Number of markers    [8]Quality (average fwd-bwd phred score)
# # ST    [2]Sample       [3]Chromosome   [4]Position     [5]State (0:HW, 1:AZ)   [6]Quality (fwd-bwd phred score)
# ST      hg002   1       13868   0       3.0
# ST      hg002   1       14210   0       99.0
# ST      hg002   1       14354   0       99.0
# ST      hg002   1       14464   0       99.0


class RohEntry:
    def __init__(self, line):
        fields = line.split("\t")
        self.sample = fields[1]
        self.chrom = fields[2]
        self.start = int(fields[3])
        self.end = int(fields[4])
        self.qual = float(fields[5])

    def get_bed_fields(self) -> list[str]:
        return [self.chrom, str(self.start), str(self.end), "ROH", ".", ".", ".", ".", "255,0,0"]


class UPDEntry:
    def __init__(self, line):
        fields = line.split("\t")
        self.chrom = fields[0]
        self.start = fields[1]
        self.end = fields[2]
        details_str = fields[3]

        details = {}
        for entry in details_str.split(";"):
            key, value = entry.split("=")
            details[key] = value
        self.details = details

    def get_bed_fields(self) -> list[str]:
        return [self.chrom, str(self.start), str(self.end), "UPD", ".", ".", ".", ".", "0,0,255"]


def main(roh: Path, upd: Path, roh_quality_threshold: int, output: Path):

    roh_entries: list[RohEntry] = parse_roh(roh, roh_quality_threshold)
    upd_entries: list[UPDEntry] = parse_upd(upd)

    with output.open("w") as out_fh:
        for entry in roh_entries:
            print("\t".join(entry.get_bed_fields()), file=out_fh)
        for entry in upd_entries:
            print("\t".join(entry.get_bed_fields()), file=out_fh)

    print(f"ROH {len(roh_entries)} UPD {len(upd_entries)}")


def parse_roh(roh_path: Path, qual_thres: float) -> list[RohEntry]:

    roh_entries: list[RohEntry] = []
    with open_file(roh_path) as roh_fh:
        for line in roh_fh:
            if line.startswith("#") or line.startswith("ST"):
                continue
            line = line.rstrip()
            roh_entry = RohEntry(line)
            if roh_entry.qual > qual_thres:
                roh_entries.append(roh_entry)
    return roh_entries


def parse_upd(upd: Path) -> list[UPDEntry]:

    upd_entries: list[UPDEntry] = []
    with open_file(upd) as upd_fh:
        for line in upd_fh:
            upd_entry = UPDEntry(line)
            upd_entries.append(upd_entry)
    return upd_entries


def write_bed():
    pass


def parse_arguments():
    parser = argparse.ArgumentParser()

    # FIXME: Understand the difference roh / upd

    parser.add_argument("--roh", required=True, type=Path)
    parser.add_argument("--upd", required=True, type=Path)

    parser.add_argument("--roh_quality_threshold", default=85, type=int)

    parser.add_argument("--output", required=True, type=Path)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    main(args.roh, args.upd, args.roh_quality_threshold, args.output)
