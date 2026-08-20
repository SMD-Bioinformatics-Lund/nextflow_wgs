#!/usr/bin/env python3

"""
Estimate duplication depth ratios from a VCF and BAM using mosdepth.
"""

import argparse
import gzip
import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, TextIO, Tuple


TARGET_CALLERS = {"manta", "tiddit"}


@dataclass
class Duplication:
    idx: int
    chrom: str
    start: int
    end: int
    variant_id: str
    callers: str
    pr_ref: Optional[int]
    pr_alt: Optional[int]
    sr_ref: Optional[int]
    sr_alt: Optional[int]


@dataclass
class RegionDepths:
    variant: Optional[float] = None
    upstream: Optional[float] = None
    downstream: Optional[float] = None


def main() -> None:
    args = parse_args()

    duplications = list(read_duplications(args.vcf, args.proband_id, args.sample_index))
    if not duplications:
        write_results([], {}, args.output, args.ploidy)
        return

    with tempfile.TemporaryDirectory(prefix="dup_depth_ratio.") as tmpdir:
        bed_path = Path(tmpdir) / "regions.bed"
        prefix = Path(tmpdir) / "mosdepth"

        write_regions_bed(
            duplications,
            bed_path,
            args.flank_size,
            args.variant_window_size,
        )
        run_mosdepth(args.mosdepth, args.bam, bed_path, prefix)
        depths = read_mosdepth_regions(f"{prefix}.regions.bed.gz")

    write_results(duplications, depths, args.output, args.ploidy)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "For DUP records in a CNV VCF, calculate variant depth, flanking "
            "depth, depth ratio, rounded ratio, and estimated copy number."
        )
    )
    parser.add_argument("--vcf", required=True, help="Input CNV VCF, optionally gzipped.")
    parser.add_argument("--bam", required=True, help="Input indexed BAM/CRAM for mosdepth.")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument(
        "--flank-size",
        type=int,
        default=None,
        help="Flank size in bp. Defaults to the duplication length for each record.",
    )
    parser.add_argument(
        "--variant-window-size",
        type=int,
        default=None,
        help=(
            "Use at most this many centered bp from each duplication for variant "
            "coverage. Defaults to the full duplication."
        ),
    )
    parser.add_argument(
        "--ploidy",
        type=int,
        default=2,
        help="Baseline ploidy used for estimated copy number. Default: 2.",
    )
    parser.add_argument(
        "--mosdepth",
        default="mosdepth",
        help="Path to mosdepth executable. Default: mosdepth.",
    )
    parser.add_argument(
        "--sample-index",
        type=int,
        default=1,
        help="Fallback 1-based sample column index to parse PR/SR from. Default: 1.",
    )
    parser.add_argument(
        "--proband-id",
        help="VCF sample ID for the proband. When set, PR/SR are parsed from this sample.",
    )
    return parser.parse_args()


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_duplications(
    vcf_path: str,
    proband_id: Optional[str] = None,
    sample_index: int = 1,
) -> Iterable[Duplication]:
    if sample_index < 1:
        raise ValueError("--sample-index must be 1 or greater")

    selected_sample_index = sample_index
    found_column_header = False

    with open_text(vcf_path) as handle:
        for idx, line in enumerate(handle, start=1):
            if line.startswith("#CHROM"):
                found_column_header = True
                selected_sample_index = resolve_sample_index(line, proband_id, sample_index)
                continue

            if line.startswith("#"):
                continue

            if not found_column_header:
                raise ValueError("Could not find VCF #CHROM header before variant records")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record at input line {idx}: fewer than 8 columns")

            chrom, pos, variant_id, _ref, alt, _qual, _filter, info = fields[:8]
            info_dict = parse_info(info)
            svtype = str(info_dict.get("SVTYPE", ""))

            if svtype != "DUP" and "<DUP>" not in alt:
                continue

            callers = parse_callers(str(info_dict.get("set", "")))
            if not callers & TARGET_CALLERS:
                continue

            pr_ref, pr_alt, sr_ref, sr_alt = parse_pr_sr(fields, selected_sample_index, idx)
            end = parse_end(info_dict, pos, idx)
            start = int(pos)
            if end < start:
                raise ValueError(f"Malformed DUP at input line {idx}: END is smaller than POS")

            yield Duplication(
                idx=idx,
                chrom=chrom,
                start=start,
                end=end,
                variant_id=variant_id,
                callers="-".join(sorted(callers)),
                pr_ref=pr_ref,
                pr_alt=pr_alt,
                sr_ref=sr_ref,
                sr_alt=sr_alt,
            )

    if not found_column_header:
        raise ValueError("Could not find VCF #CHROM header")


def parse_info(info: str) -> Dict[str, object]:
    parsed: Dict[str, object] = {}
    if info == ".":
        return parsed

    for entry in info.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            parsed[key] = value
        else:
            parsed[entry] = True
    return parsed


def parse_end(info: Dict[str, object], pos: str, line_number: int) -> int:
    if "END" in info:
        return int(str(info["END"]))
    if "SVLEN" in info:
        svlen = abs(int(str(info["SVLEN"]).split(",", 1)[0]))
        return int(pos) + svlen - 1
    raise ValueError(f"DUP at input line {line_number} has neither END nor SVLEN in INFO")


def parse_callers(set_value: str) -> Set[str]:
    if not set_value or set_value == ".":
        return set()
    return {caller.strip().lower() for caller in set_value.split("-") if caller.strip()}


def resolve_sample_index(
    column_header_line: str,
    proband_id: Optional[str],
    fallback_sample_index: int,
) -> int:
    fields = column_header_line.rstrip("\n").split("\t")
    sample_ids = fields[9:]

    if proband_id is None:
        if sample_ids and fallback_sample_index > len(sample_ids):
            raise ValueError(
                f"--sample-index {fallback_sample_index} was requested, but VCF has {len(sample_ids)} samples"
            )
        return fallback_sample_index

    if proband_id not in sample_ids:
        raise ValueError(
            f"Proband ID '{proband_id}' was not found in VCF samples: {', '.join(sample_ids)}"
        )

    return sample_ids.index(proband_id) + 1


def parse_pr_sr(
    fields: List[str],
    sample_index: int,
    line_number: int,
) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[int]]:
    if len(fields) < 10:
        return None, None, None, None

    sample_column = 8 + sample_index
    if len(fields) <= sample_column:
        raise ValueError(f"Input line {line_number} does not have sample column {sample_index}")

    format_keys = fields[8].split(":")
    sample_values = fields[sample_column].split(":")
    sample_data = dict(zip(format_keys, sample_values))
    pr_ref, pr_alt = parse_ref_alt_count(sample_data.get("PR"))
    sr_ref, sr_alt = parse_ref_alt_count(sample_data.get("SR"))
    return pr_ref, pr_alt, sr_ref, sr_alt


def parse_ref_alt_count(value: Optional[str]) -> Tuple[Optional[int], Optional[int]]:
    if value is None or value in (".", ""):
        return None, None

    counts = value.split(",")
    if len(counts) < 2 or "." in counts[:2]:
        return None, None

    return int(counts[0]), int(counts[1])


def region_name(dup: Duplication, region_type: str) -> str:
    return f"dup{dup.idx}|{region_type}"


def write_regions_bed(
    duplications: List[Duplication],
    bed_path: Path,
    flank_size: Optional[int],
    variant_window_size: Optional[int] = None,
) -> None:
    if flank_size is not None and flank_size <= 0:
        raise ValueError("--flank-size must be positive")
    if variant_window_size is not None and variant_window_size <= 0:
        raise ValueError("--variant-window-size must be positive")

    with open(bed_path, "w") as bed:
        for dup in duplications:
            event_start = dup.start - 1
            event_end = dup.end
            event_length = event_end - event_start
            flank = flank_size if flank_size is not None else event_length
            variant_start, variant_end = centered_variant_window(
                event_start,
                event_end,
                variant_window_size,
            )

            upstream_start = max(0, event_start - flank)
            upstream_end = event_start
            downstream_start = event_end
            downstream_end = event_end + flank

            if upstream_start < upstream_end:
                bed.write(
                    f"{dup.chrom}\t{upstream_start}\t{upstream_end}\t{region_name(dup, 'upstream')}\n"
                )

            bed.write(
                f"{dup.chrom}\t{variant_start}\t{variant_end}\t{region_name(dup, 'variant')}\n"
            )
            bed.write(
                f"{dup.chrom}\t{downstream_start}\t{downstream_end}\t{region_name(dup, 'downstream')}\n"
            )


def centered_variant_window(
    event_start: int,
    event_end: int,
    variant_window_size: Optional[int],
) -> Tuple[int, int]:
    event_length = event_end - event_start
    if event_length <= 0:
        raise ValueError("Duplication interval length must be positive")
    if variant_window_size is None or event_length <= variant_window_size:
        return event_start, event_end

    center = event_start + event_length // 2
    window_start = center - variant_window_size // 2
    window_end = window_start + variant_window_size
    return window_start, window_end


def run_mosdepth(mosdepth: str, bam: str, bed_path: Path, prefix: Path) -> None:
    command = [
        mosdepth,
        "--no-per-base",
        "--fast-mode",
        "--by",
        str(bed_path),
        str(prefix),
        bam,
    ]
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as error:
        raise FileNotFoundError(f"Could not find mosdepth executable: {mosdepth}") from error


def read_mosdepth_regions(regions_path: str) -> Dict[int, RegionDepths]:
    depths: Dict[int, RegionDepths] = {}
    with gzip.open(regions_path, "rt") as handle:
        for line in handle:
            chrom, start, end, name, mean_depth = line.rstrip("\n").split("\t")[:5]
            del chrom, start, end
            dup_label, region_type = name.split("|", 1)
            dup_idx = int(dup_label[3:])
            record_depths = depths.setdefault(dup_idx, RegionDepths())
            setattr(record_depths, region_type, float(mean_depth))
    return depths


def write_results(
    duplications: List[Duplication],
    depths: Dict[int, RegionDepths],
    output_path: str,
    ploidy: int,
) -> None:
    with open(output_path, "w") as out:
        out.write(
            "\t".join(
                [
                    "chrom",
                    "start",
                    "end",
                    "id",
                    "callers",
                    "variant_mean_depth",
                    "upstream_mean_depth",
                    "downstream_mean_depth",
                    "surrounding_mean_depth",
                    "depth_ratio",
                    "rounded_depth_ratio",
                    "estimated_copy_number",
                    "pr_ref_count",
                    "pr_alt_count",
                    "pr_alt_ref_ratio",
                    "sr_ref_count",
                    "sr_alt_count",
                    "sr_alt_ref_ratio",
                ]
            )
            + "\n"
        )

        for dup in duplications:
            record_depths = depths.get(dup.idx, RegionDepths())
            surrounding = mean_ignore_none([record_depths.upstream, record_depths.downstream])
            ratio = safe_divide(record_depths.variant, surrounding)
            rounded_ratio = nearest_int(ratio) if ratio is not None else None
            copy_number = nearest_int(ratio * ploidy) if ratio is not None else None
            pr_ratio = safe_divide_ints(dup.pr_alt, dup.pr_ref)
            sr_ratio = safe_divide_ints(dup.sr_alt, dup.sr_ref)

            out.write(
                "\t".join(
                    [
                        dup.chrom,
                        str(dup.start),
                        str(dup.end),
                        dup.variant_id,
                        dup.callers,
                        format_optional_float(record_depths.variant),
                        format_optional_float(record_depths.upstream),
                        format_optional_float(record_depths.downstream),
                        format_optional_float(surrounding),
                        format_optional_float(ratio),
                        format_optional_int(rounded_ratio),
                        format_optional_int(copy_number),
                        format_optional_int(dup.pr_ref),
                        format_optional_int(dup.pr_alt),
                        format_optional_float(pr_ratio),
                        format_optional_int(dup.sr_ref),
                        format_optional_int(dup.sr_alt),
                        format_optional_float(sr_ratio),
                    ]
                )
                + "\n"
            )


def mean_ignore_none(values: List[Optional[float]]) -> Optional[float]:
    real_values = [value for value in values if value is not None]
    if not real_values:
        return None
    return sum(real_values) / len(real_values)


def safe_divide(numerator: Optional[float], denominator: Optional[float]) -> Optional[float]:
    if numerator is None or denominator in (None, 0):
        return None
    return numerator / denominator


def safe_divide_ints(numerator: Optional[int], denominator: Optional[int]) -> Optional[float]:
    if numerator is None or denominator in (None, 0):
        return None
    return numerator / denominator


def nearest_int(value: float) -> int:
    return math.floor(value + 0.5)


def format_optional_float(value: Optional[float]) -> str:
    if value is None:
        return "NA"
    return f"{value:.4f}"


def format_optional_int(value: Optional[int]) -> str:
    if value is None:
        return "NA"
    return str(value)


if __name__ == "__main__":
    main()
