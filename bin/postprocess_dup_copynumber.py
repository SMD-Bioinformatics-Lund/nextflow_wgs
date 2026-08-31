#!/usr/bin/env python3

"""
This script is needed to make geneticists aware of CN 4 and above
which can be missed in single-analysis without parents.

Add copy-number estimates to duplication records in a VCF.
Uses different logic depending on what callers are merged.

1. GATK-only should not be adjusted, since they should have CN
2. GATK merged with manta and/or tiddit gets gatkCN -> format CN
3. manta tiddit combinations gets new CN estimated from mosdepth
   by calculating varcov/mean flank cov -> closest integer

The script can take any variant callers by providing what cov-based
caller with INFO-CN-flag and readpair based caller(s). 

This script need to be run after SVDB merging of trio calls
(postprocess_vep_sv). 

Alternative to this method could be:
only estimate manta calls. Let SVDB merge CN in format for all
and then re-prioritize CN calls from GATK where possible. This
would retain CN calls for parents where GATK is not the only caller
"""

import argparse
import gzip
import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, TextIO, Tuple


@dataclass
class Duplication:
    idx: int
    chrom: str
    start: int
    end: int


@dataclass
class RegionDepths:
    variant: Optional[float] = None
    upstream: Optional[float] = None
    downstream: Optional[float] = None


def main() -> None:
    args = parse_args()
    copy_number_callers = parse_caller_list(args.copy_number_callers)
    readpair_callers = parse_caller_list(args.readpair_callers)

    duplications = list(
        read_readpair_duplications(
            args.vcf,
            copy_number_callers,
            readpair_callers,
        )
    )

    depths: Dict[int, RegionDepths] = {}
    if duplications:
        if args.bam is None:
            raise ValueError("--bam is required when readpair-aware duplications need CN estimation")

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

    write_vcf(
        args.vcf,
        args.output,
        depths,
        args.ploidy,
        copy_number_callers,
        readpair_callers,
        args.copy_number_info_field,
        args.proband_id,
        args.sample_index,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "For DUP/TDUP records in a VCF, add FORMAT/CN from either a "
            "copy-number-aware caller INFO field or a mosdepth depth estimate."
        )
    )
    parser.add_argument("--vcf", required=True, help="Input VCF, optionally gzipped.")
    parser.add_argument(
        "--bam",
        help="Input indexed BAM/CRAM for mosdepth. Required for readpair-only duplications.",
    )
    parser.add_argument("--output", required=True, help="Output VCF path.")
    parser.add_argument(
        "--copy-number-callers",
        default="gatk",
        help="Comma-separated callers that already provide copy number. Default: gatk.",
    )
    parser.add_argument(
        "--readpair-callers",
        default="manta,tiddit",
        help="Comma-separated readpair-aware callers to estimate CN for. Default: manta,tiddit.",
    )
    parser.add_argument(
        "--copy-number-info-field",
        default="gatkCN",
        help="INFO field containing copy number for CN-aware merged calls. Default: gatkCN.",
    )
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
        help="Fallback 1-based sample column index to update with FORMAT/CN. Default: 1.",
    )
    parser.add_argument(
        "--proband-id",
        help="VCF sample ID to update with FORMAT/CN. Overrides --sample-index.",
    )
    return parser.parse_args()


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def open_output_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def read_readpair_duplications(
    vcf_path: str,
    copy_number_callers: Set[str],
    readpair_callers: Set[str],
) -> Iterable[Duplication]:
    found_column_header = False

    with open_text(vcf_path) as handle:
        for idx, line in enumerate(handle, start=1):
            if line.startswith("#CHROM"):
                found_column_header = True
                continue

            if line.startswith("#"):
                continue

            if not found_column_header:
                raise ValueError("Could not find VCF #CHROM header before variant records")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record at input line {idx}: fewer than 8 columns")

            chrom, pos, _variant_id, _ref, alt, _qual, _filter, info = fields[:8]
            info_dict = parse_info(info)
            callers = parse_callers(str(info_dict.get("set", "")))

            if not is_duplication_record(info_dict, alt):
                continue
            if callers & copy_number_callers:
                continue
            if not callers & readpair_callers:
                continue

            end = parse_end(info_dict, pos, idx)
            start = int(pos)
            if end < start:
                raise ValueError(f"Malformed DUP/TDUP at input line {idx}: END is smaller than POS")

            yield Duplication(idx=idx, chrom=chrom, start=start, end=end)

    if not found_column_header:
        raise ValueError("Could not find VCF #CHROM header")


def write_vcf(
    vcf_path: str,
    output_path: str,
    depths: Dict[int, RegionDepths],
    ploidy: int,
    copy_number_callers: Set[str],
    readpair_callers: Set[str],
    copy_number_info_field: str,
    proband_id: Optional[str],
    sample_index: int,
) -> None:
    if sample_index < 1:
        raise ValueError("--sample-index must be 1 or greater")

    selected_sample_index = sample_index
    found_column_header = False
    has_cn_header = False
    has_estimated_cn_header = False
    needs_cn_header = vcf_needs_cn_header(vcf_path, copy_number_callers, readpair_callers)
    needs_estimated_cn_header = vcf_needs_estimated_cn_header(
        vcf_path,
        copy_number_callers,
        readpair_callers,
    )

    with open_text(vcf_path) as inp, open_output_text(output_path) as out:
        for idx, line in enumerate(inp, start=1):
            if line.startswith("##FORMAT=<ID=CN,"):
                has_cn_header = True
            if line.startswith("##INFO=<ID=ESTIMATED_CN,"):
                has_estimated_cn_header = True

            if line.startswith("#CHROM"):
                found_column_header = True
                selected_sample_index = resolve_sample_index(line, proband_id, sample_index)
                if needs_cn_header and not has_cn_header:
                    out.write(
                        '##FORMAT=<ID=CN,Number=.,Type=Integer,Description="Copy number">\n'
                    )
                if needs_estimated_cn_header and not has_estimated_cn_header:
                    out.write(
                        "##INFO=<ID=ESTIMATED_CN,Number=5,Type=Float,"
                        'Description="Mosdepth-estimated copy number and mean coverage values: '
                        'CN,variant,upstream_flank,downstream_flank,mean_flanks">\n'
                    )
                out.write(line)
                continue

            if line.startswith("#"):
                out.write(line)
                continue

            if not found_column_header:
                raise ValueError("Could not find VCF #CHROM header before variant records")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record at input line {idx}: fewer than 8 columns")

            info_dict = parse_info(fields[7])
            callers = parse_callers(str(info_dict.get("set", "")))
            cn_value = copy_number_for_record(
                idx,
                fields[4],
                info_dict,
                callers,
                depths,
                ploidy,
                copy_number_callers,
                readpair_callers,
                copy_number_info_field,
            )

            if cn_value is not None:
                estimated_cn_value = estimated_cn_info_for_record(
                    idx,
                    fields[4],
                    info_dict,
                    callers,
                    depths,
                    ploidy,
                    copy_number_callers,
                    readpair_callers,
                )
                if estimated_cn_value is not None:
                    fields[7] = add_info_value(fields[7], "ESTIMATED_CN", estimated_cn_value)
                fields = add_sample_format_value(fields, selected_sample_index, "CN", cn_value, idx)
                out.write("\t".join(fields) + "\n")
            else:
                out.write(line)

    if not found_column_header:
        raise ValueError("Could not find VCF #CHROM header")


def vcf_needs_cn_header(
    vcf_path: str,
    copy_number_callers: Set[str],
    readpair_callers: Set[str],
) -> bool:
    found_column_header = False

    with open_text(vcf_path) as handle:
        for idx, line in enumerate(handle, start=1):
            if line.startswith("#CHROM"):
                found_column_header = True
                continue

            if line.startswith("#"):
                continue

            if not found_column_header:
                raise ValueError("Could not find VCF #CHROM header before variant records")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record at input line {idx}: fewer than 8 columns")

            info_dict = parse_info(fields[7])
            callers = parse_callers(str(info_dict.get("set", "")))
            if not is_duplication_record(info_dict, fields[4]):
                continue

            has_copy_number_caller = bool(callers & copy_number_callers)
            has_readpair_caller = bool(callers & readpair_callers)
            if has_copy_number_caller and callers - copy_number_callers:
                return True
            if not has_copy_number_caller and has_readpair_caller:
                return True

    if not found_column_header:
        raise ValueError("Could not find VCF #CHROM header")

    return False


def vcf_needs_estimated_cn_header(
    vcf_path: str,
    copy_number_callers: Set[str],
    readpair_callers: Set[str],
) -> bool:
    found_column_header = False

    with open_text(vcf_path) as handle:
        for idx, line in enumerate(handle, start=1):
            if line.startswith("#CHROM"):
                found_column_header = True
                continue

            if line.startswith("#"):
                continue

            if not found_column_header:
                raise ValueError("Could not find VCF #CHROM header before variant records")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record at input line {idx}: fewer than 8 columns")

            info_dict = parse_info(fields[7])
            callers = parse_callers(str(info_dict.get("set", "")))
            if not is_duplication_record(info_dict, fields[4]):
                continue

            if not callers & copy_number_callers and callers & readpair_callers:
                return True

    if not found_column_header:
        raise ValueError("Could not find VCF #CHROM header")

    return False


def copy_number_for_record(
    idx: int,
    alt: str,
    info_dict: Dict[str, object],
    callers: Set[str],
    depths: Dict[int, RegionDepths],
    ploidy: int,
    copy_number_callers: Set[str],
    readpair_callers: Set[str],
    copy_number_info_field: str,
) -> Optional[str]:
    if not is_duplication_record(info_dict, alt):
        return None

    has_copy_number_caller = bool(callers & copy_number_callers)
    has_readpair_caller = bool(callers & readpair_callers)
    is_merged_copy_number_call = bool(callers - copy_number_callers)

    if has_copy_number_caller and is_merged_copy_number_call:
        if copy_number_info_field not in info_dict:
            raise ValueError(
                f"DUP/TDUP at input line {idx} is merged with a CN-aware caller but "
                f"does not have INFO/{copy_number_info_field}"
            )
        return first_list_value(str(info_dict[copy_number_info_field]))

    if has_copy_number_caller:
        return None

    if has_readpair_caller:
        copy_number = estimate_copy_number(depths.get(idx, RegionDepths()), ploidy)
        return "." if copy_number is None else str(copy_number)

    return None


def estimated_cn_info_for_record(
    idx: int,
    alt: str,
    info_dict: Dict[str, object],
    callers: Set[str],
    depths: Dict[int, RegionDepths],
    ploidy: int,
    copy_number_callers: Set[str],
    readpair_callers: Set[str],
) -> Optional[str]:
    if not is_duplication_record(info_dict, alt):
        return None

    has_copy_number_caller = bool(callers & copy_number_callers)
    has_readpair_caller = bool(callers & readpair_callers)
    if has_copy_number_caller or not has_readpair_caller:
        return None

    record_depths = depths.get(idx, RegionDepths())
    return format_estimated_cn_info(record_depths, ploidy)


def format_estimated_cn_info(record_depths: RegionDepths, ploidy: int) -> str:
    surrounding = mean_ignore_none([record_depths.upstream, record_depths.downstream])
    copy_number = estimate_copy_number(record_depths, ploidy)
    return ",".join(
        [
            format_optional_int(copy_number),
            format_optional_float(record_depths.variant),
            format_optional_float(record_depths.upstream),
            format_optional_float(record_depths.downstream),
            format_optional_float(surrounding),
        ]
    )


def add_info_value(info: str, key: str, value: str) -> str:
    entries = [] if info in ("", ".") else info.split(";")
    replaced = False
    for entry_index, entry in enumerate(entries):
        if entry == key or entry.startswith(f"{key}="):
            entries[entry_index] = f"{key}={value}"
            replaced = True
            break
    if not replaced:
        entries.append(f"{key}={value}")
    return ";".join(entries) if entries else "."


def add_sample_format_value(
    fields: List[str],
    sample_index: int,
    key: str,
    value: str,
    line_number: int,
) -> List[str]:
    if len(fields) < 10:
        raise ValueError(
            f"Cannot add FORMAT/{key} at input line {line_number}: VCF record has no sample columns"
        )

    sample_column = 8 + sample_index
    if len(fields) <= sample_column:
        raise ValueError(f"Input line {line_number} does not have sample column {sample_index}")

    format_keys = [] if fields[8] in (".", "") else fields[8].split(":")
    if key in format_keys:
        key_index = format_keys.index(key)
    else:
        format_keys.append(key)
        key_index = len(format_keys) - 1
        fields[8] = ":".join(format_keys)

    for column in range(9, len(fields)):
        sample_values = fields[column].split(":")
        while len(sample_values) < len(format_keys):
            sample_values.append(".")
        if column == sample_column:
            sample_values[key_index] = value
        fields[column] = ":".join(sample_values)

    return fields


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
    raise ValueError(f"DUP/TDUP at input line {line_number} has neither END nor SVLEN in INFO")


def parse_callers(set_value: str) -> Set[str]:
    if not set_value or set_value == ".":
        return set()
    return {caller.strip().lower() for caller in set_value.split("-") if caller.strip()}


def parse_caller_list(callers: str) -> Set[str]:
    return {caller.strip().lower() for caller in callers.split(",") if caller.strip()}


def first_list_value(value: str) -> str:
    return value.split(",", 1)[0]


def is_duplication_record(info: Dict[str, object], alt: str) -> bool:
    svtype = str(info.get("SVTYPE", ""))
    return svtype in {"DUP", "TDUP"} or "<DUP>" in alt or "<TDUP>" in alt


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


def estimate_copy_number(record_depths: RegionDepths, ploidy: int) -> Optional[int]:
    surrounding = mean_ignore_none([record_depths.upstream, record_depths.downstream])
    ratio = safe_divide(record_depths.variant, surrounding)
    if ratio is None:
        return None
    return nearest_int(ratio * ploidy)


def mean_ignore_none(values: List[Optional[float]]) -> Optional[float]:
    real_values = [value for value in values if value is not None]
    if not real_values:
        return None
    return sum(real_values) / len(real_values)


def safe_divide(numerator: Optional[float], denominator: Optional[float]) -> Optional[float]:
    if numerator is None or denominator in (None, 0):
        return None
    return numerator / denominator


def nearest_int(value: float) -> int:
    return math.floor(value + 0.5)


def format_optional_float(value: Optional[float]) -> str:
    if value is None:
        return "."
    return f"{value:.4f}"


def format_optional_int(value: Optional[int]) -> str:
    if value is None:
        return "."
    return str(value)


if __name__ == "__main__":
    main()
