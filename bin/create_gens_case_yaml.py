#!/usr/bin/env python3
"""Create Gens case YAML for `gens load case` from per-sample inputs."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any


TYPE_ORDER = {"proband": 0, "mother": 1, "father": 2}
PLACE_LAST = math.inf
EMPTY_VALUES = {"", "none", "null", "."}


def main(
    case_id: str,
    gens_accessdir: str,
    sample_ids: list[str],
    sample_types: list[str],
    sexes: list[str],
    roh_tracks: list[str],
    upd_tracks: list[str],
    meta_files: list[str],
    chrom_meta_files: list[str],
    output: Path,
    is_trio: bool,
) -> None:
    samples = load_samples_from_cli_args(
        sample_ids=sample_ids,
        sample_types=sample_types,
        sexes=sexes,
        roh_tracks=roh_tracks,
        upd_tracks=upd_tracks,
        meta_files=meta_files,
        chrom_meta_files=chrom_meta_files,
    )

    yaml_lines = build_yaml_lines(
        case_id=case_id,
        gens_accessdir=gens_accessdir,
        samples=samples,
        trio=is_trio,
    )
    output.write_text("\n".join(yaml_lines) + "\n", encoding="utf-8")


def load_samples_from_cli_args(
    sample_ids: list[str],
    sample_types: list[str],
    sexes: list[str],
    roh_tracks: list[str],
    upd_tracks: list[str],
    meta_files: list[str],
    chrom_meta_files: list[str],
) -> list[Sample]:
    """Generate per-sample objects from grouped CLI arguments."""
    fields = {
        "sample_id": sample_ids,
        "sample_type": sample_types,
        "sex": sexes,
        "roh_track": roh_tracks,
        "upd_track": upd_tracks,
        "meta_file": meta_files,
        "chrom_meta_file": chrom_meta_files,
    }
    lengths = {name: len(values or []) for name, values in fields.items()}
    unique_lengths = set(lengths.values())
    if len(unique_lengths) != 1:
        raise ValueError(
            "All sample list options must be provided the same number of values: "
            + ", ".join([f"{key}={value}" for key, value in sorted(lengths.items())])
        )

    sample_count = unique_lengths.pop()
    if sample_count == 0:
        raise ValueError("At least one sample must be provided")

    samples: list[Sample] = []
    for idx in range(sample_count):
        sample_id = normalize(sample_ids[idx])
        if sample_id is None:
            continue
        samples.append(
            Sample(
                sample_id=sample_id,
                sample_type=normalize(sample_types[idx]),
                sex=normalize(sexes[idx]),
                roh_track=normalize(roh_tracks[idx]),
                upd_track=normalize(upd_tracks[idx]),
                meta_file=normalize(meta_files[idx]),
                chrom_meta_file=normalize(chrom_meta_files[idx]),
            )
        )

    if not samples:
        raise ValueError("No sample rows found in CLI sample options")

    samples.sort(
        key=lambda sample: (
            TYPE_ORDER.get(sample.sample_type or "", 99),
            sample.sample_id,
        )
    )
    return samples


def build_yaml_lines(
    case_id: str,
    gens_accessdir: str,
    samples: list[Sample],
    trio: bool,
) -> list[str]:
    lines: list[str] = [
        f"case_id: '{case_id}'",
        "genome_build: 38",
        "samples:",
    ]

    for sample in samples:
        sample_id = sample.sample_id
        lines.append(f"  - sample_id: '{sample_id}'")
        lines.append(
            f"    baf: '{path_in_accessdir(gens_accessdir, f'{sample_id}.baf.bed.gz')}'"
        )
        lines.append(
            f"    coverage: '{path_in_accessdir(gens_accessdir, f'{sample_id}.cov.bed.gz')}'"
        )

        if sample.sample_type is not None:
            lines.append(f"    sample_type: '{sample.sample_type}'")
        if sample.sex is not None:
            lines.append(f"    sex: '{sample.sex}'")

        if sample.sample_type == "proband":
            meta_paths = [
                path_in_accessdir(gens_accessdir, value)
                for value in (sample.meta_file, sample.chrom_meta_file)
                if value is not None
            ]
            if meta_paths:
                lines.append("    meta_files:")
                for path in meta_paths:
                    lines.append(f"      - '{path}'")

            annotation_rows: list[dict[str, Any]] = []
            if sample.roh_track is not None:
                annotation_rows.append(
                    {
                        "file": path_in_accessdir(gens_accessdir, sample.roh_track),
                        "name": "LOH",
                    }
                )
            if trio and sample.upd_track is not None:
                annotation_rows.append(
                    {
                        "file": path_in_accessdir(gens_accessdir, sample.upd_track),
                        "name": "UPS",
                    }
                )
            if annotation_rows:
                lines.append("    sample_annotations:")
                for annotation in annotation_rows:
                    lines.append(f"      - file: '{annotation['file']}'")
                    lines.append(f"        name: '{annotation['name']}'")

    return lines


def normalize(value: str | None) -> str | None:
    """Remove whitespaces and return None for empty string values."""
    if value is None:
        return None
    stripped = value.strip()
    if stripped.lower() in EMPTY_VALUES:
        return None
    return stripped


def path_in_accessdir(gens_accessdir: str, filename_or_path: str) -> str:
    return f"{gens_accessdir.rstrip('/')}/{Path(filename_or_path).name}"


@dataclass
class Sample:
    sample_id: str
    sample_type: str | None
    sex: str | None
    roh_track: str | None
    upd_track: str | None
    meta_file: str | None
    chrom_meta_file: str | None


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case_id", required=True)
    parser.add_argument("--gens_accessdir", required=True)
    parser.add_argument("--sample_ids", nargs="+", required=True)
    parser.add_argument("--sample_types", nargs="+")
    parser.add_argument("--sexes", nargs="+")
    parser.add_argument("--roh_tracks", nargs="+")
    parser.add_argument("--upd_tracks", nargs="+")
    parser.add_argument("--meta_files", nargs="+")
    parser.add_argument("--chrom_meta_files", nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--trio", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    main(
        args.case_id,
        args.gens_accessdir,
        args.sample_ids,
        args.sample_types or [],
        args.sexes or [],
        args.roh_tracks or [],
        args.upd_tracks or [],
        args.meta_files or [],
        args.chrom_meta_files or [],
        args.output,
        args.trio,
    )
