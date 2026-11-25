#!/usr/bin/env python3
"""Build a singleton VCF and pedigree for an isolated variant score."""

import argparse
import gzip
from pathlib import Path

SAMPLE = "query"
BASES = {"A", "C", "G", "T"}


def read_record(path):
    with path.open("rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rt") as handle:
        records = (line.rstrip("\r\n") for line in handle if line.strip() and not line.startswith("#"))
        record = next(records, None)
        if record is None or next(records, None) is not None:
            raise ValueError("expected exactly one VCF record")
    fields = record.split("\t")
    if len(fields) < 8:
        raise ValueError("expected the eight tab-separated VCF columns")
    return fields[:8]


def info_fields(info):
    return dict(item.partition("=")[::2] for item in info.split(";") if item != ".")


def reference_base(fasta, chrom, pos):
    index = {}
    with open(f"{fasta}.fai") as handle:
        for line in handle:
            name, length, offset, bases, width, *_ = line.split("\t")
            index[name] = (int(length), int(offset), int(bases), int(width))

    bare = chrom.removeprefix("chr")
    aliases = (chrom, bare, f"chr{bare}")
    if bare in {"M", "MT"}:
        aliases += ("M", "MT", "chrM", "chrMT")
    contig = next((name for name in aliases if name in index), None)
    if contig is None:
        raise ValueError(f"contig {chrom!r} is absent from the reference FASTA")

    length, offset, bases, width = index[contig]
    if not 1 <= pos <= length:
        raise ValueError(f"POS must be between 1 and {length} for {contig}")
    row, column = divmod(pos - 1, bases)
    with fasta.open("rb") as handle:
        handle.seek(offset + row * width + column)
        base = handle.read(1).decode("ascii").upper()
    if base not in BASES:
        raise ValueError(f"expected an A/C/G/T reference base at {contig}:{pos}")
    return contig, length, base


def prepare(args):
    chrom, position, identifier, ref, alt, *_ = read_record(args.variant)
    pos = int(position)
    contig, contig_length, base = reference_base(args.reference, chrom, pos)
    ref, alt = ref.upper(), alt.upper()
    if ref not in BASES or alt not in BASES or ref == alt:
        raise ValueError("expected a biallelic SNV with different A/C/G/T alleles")
    if contig in {"M", "MT", "chrM", "chrMT"}:
        contig = "MT"

    if ref != base:
        raise ValueError(f"REF {ref!r} does not match reference base {base} at {contig}:{pos}")
    if not identifier or any(char.isspace() for char in identifier):
        raise ValueError("VCF ID must be nonempty and contain no whitespace")

    headers = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={contig},length={contig_length}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}",
    ]
    record = [contig, str(pos), identifier, ref, alt, ".", "PASS", ".", "GT", args.genotype]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "snv.vcf").write_text("\n".join([*headers, "\t".join(record)]) + "\n")
    (args.output_dir / "query.ped").write_text(f"{SAMPLE}\t{SAMPLE}\t0\t0\t{args.sex}\t2\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", type=Path, required=True, help="VCF or single tab-separated VCF row")
    parser.add_argument("--reference", type=Path, required=True, help="Uncompressed GRCh38 FASTA with a .fai index")
    parser.add_argument("--genotype", choices=["0/1", "1/1"], default="0/1")
    parser.add_argument("--sex", choices=["0", "1", "2"], default="0", help="PED sex code (0 unknown, 1 male, 2 female)")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        prepare(args)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
