#!/usr/bin/env python3
"""Extract the query family's numeric RankScore from a scored singleton VCF."""

import argparse
from decimal import Decimal, InvalidOperation
from pathlib import Path

from prepare_variant_query import SAMPLE, info_fields, read_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", type=Path)
    args = parser.parse_args()
    try:
        chrom, pos, _, ref, alt, _, _, info = read_record(args.vcf)
        fields = info_fields(info)
        if fields.get("CSQ") in {None, ".", ""}:
            raise ValueError("the scored variant is missing its VEP consequence annotation")
        scores = dict(entry.split(":", 1) for entry in fields["RankScore"].split(","))
        score = scores[SAMPLE]
        if not Decimal(score).is_finite():
            raise ValueError("RankScore must be finite")
    except (OSError, KeyError, ValueError, InvalidOperation) as error:
        parser.error(f"could not extract a RankScore for {SAMPLE}: {error}")

    print("CHROM\tPOS\tEND\tREF\tALT\tRankScore")
    print(chrom, pos, fields.get("END", pos), ref, alt, score, sep="\t")


if __name__ == "__main__":
    main()
