#!/usr/bin/env python3
"""Create Gens case YAML for `gens load case` from per-sample inputs."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any


TYPE_ORDER = {"proband": 0, "mother": 1, "father": 2}
PLACE_LAST = math.inf
EMPTY_VALUES = {"", "none", "null", "."}


def normalize(value: str | None) -> str | None:
    """Remove whitespaces and return None for empty string values"""
    if value is None:
        return None
    stripped = value.strip()
    if stripped.lower() in EMPTY_VALUES:
        return None
    return stripped


def sort_samples(samples: list[dict[str, str | None]]) -> list[dict[str, str | None]]:
    samples.sort(
        key=lambda sample: (
            TYPE_ORDER.get(sample["sample_type"] or "", PLACE_LAST),
            sample["sample_id"] or "",
        )
    )
    return samples


def load_ped_roles(ped: Path) -> dict[str, str]:
    rows: list[tuple[str, str, str, str]] = []
    with ped.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            sample_id = normalize(fields[1])
            father = normalize(fields[2]) or "0"
            mother = normalize(fields[3]) or "0"
            phenotype = normalize(fields[5]) or "0"
            if sample_id is None:
                continue
            rows.append((sample_id, father, mother, phenotype))

    if not rows:
        return {}

    role_by_id: dict[str, str] = {}
    father_ids = {father for _, father, _, _ in rows if father not in EMPTY_VALUES and father != "0"}
    mother_ids = {mother for _, _, mother, _ in rows if mother not in EMPTY_VALUES and mother != "0"}
    for sample_id in father_ids:
        role_by_id[sample_id] = "father"
    for sample_id in mother_ids:
        role_by_id[sample_id] = "mother"

    # Prefer an affected sample with known parents. Fall back to any affected sample.
    proband_id: str | None = None
    for sample_id, father, mother, phenotype in rows:
        if phenotype == "2" and (father != "0" or mother != "0"):
            proband_id = sample_id
            break
    if proband_id is None:
        for sample_id, _father, _mother, phenotype in rows:
            if phenotype == "2":
                proband_id = sample_id
                break
    if proband_id is None and len(rows) == 1:
        proband_id = rows[0][0]
    if proband_id is not None:
        role_by_id[proband_id] = "proband"

    return role_by_id


def apply_ped_roles(
    samples: list[dict[str, str | None]], ped_roles: dict[str, str] | None
) -> list[dict[str, str | None]]:
    if not ped_roles:
        return sort_samples(samples)

    updated_samples: list[dict[str, str | None]] = []
    for sample in samples:
        sample_id = sample["sample_id"]
        if sample_id is None or sample_id not in ped_roles:
            updated_samples.append(sample)
            continue
        updated = dict(sample)
        updated["sample_type"] = ped_roles[sample_id]
        updated_samples.append(updated)

    return sort_samples(updated_samples)


def load_samples_from_cli_args(
    sample_ids: list[str] | None,
    sample_types: list[str] | None,
    sexes: list[str] | None,
    roh_tracks: list[str] | None,
    upd_tracks: list[str] | None,
    meta_files: list[str] | None,
    chrom_meta_files: list[str] | None,
) -> list[dict[str, str | None]]:
    """Generate a per-sample dict"""
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

    samples: list[dict[str, str | None]] = []
    for idx in range(sample_count):
        assert sample_ids is not None
        assert sample_types is not None
        assert sexes is not None
        assert roh_tracks is not None
        assert upd_tracks is not None
        assert meta_files is not None
        assert chrom_meta_files is not None
        sample_id = normalize(sample_ids[idx])
        if sample_id is None:
            continue
        samples.append(
            {
                "sample_id": sample_id,
                "sample_type": normalize(sample_types[idx]),
                "sex": normalize(sexes[idx]),
                "roh_track": normalize(roh_tracks[idx]),
                "upd_track": normalize(upd_tracks[idx]),
                "meta_file": normalize(meta_files[idx]),
                "chrom_meta_file": normalize(chrom_meta_files[idx]),
            }
        )

    if not samples:
        raise ValueError("No sample rows found in CLI sample options")

    return sort_samples(samples)


def select_samples(samples: list[dict[str, str | None]], trio: bool) -> list[dict[str, str | None]]:
    if trio:
        by_type: dict[str, dict[str, str | None]] = {}
        for sample in samples:
            sample_type = sample["sample_type"]
            if sample_type in TYPE_ORDER and sample_type not in by_type:
                by_type[sample_type] = sample
        selected = [by_type[key] for key in ("proband", "mother", "father") if key in by_type]
        if selected:
            return selected[:3]
        return samples[:3]

    for sample in samples:
        if sample["sample_type"] == "proband":
            return [sample]
    return samples[:1]


def path_in_accessdir(gens_accessdir: str, filename_or_path: str) -> str:
    return f"{gens_accessdir.rstrip('/')}/{Path(filename_or_path).name}"


def build_yaml_lines(
    case_id: str,
    gens_accessdir: str,
    selected_samples: list[dict[str, str | None]],
    trio: bool,
) -> list[str]:
    lines: list[str] = [
        f"case_id: '{case_id}'",
        "genome_build: 38",
        "samples:",
    ]

    for sample in selected_samples:
        sample_id = sample["sample_id"]
        if sample_id is None:
            continue
        lines.append(f"  - sample_id: '{sample_id}'")
        lines.append(
            f"    baf: '{path_in_accessdir(gens_accessdir, f'{sample_id}.baf.bed.gz')}'"
        )
        lines.append(
            f"    coverage: '{path_in_accessdir(gens_accessdir, f'{sample_id}.cov.bed.gz')}'"
        )

        if sample["sample_type"] is not None:
            lines.append(f"    sample_type: '{sample['sample_type']}'")
        if sample["sex"] is not None:
            lines.append(f"    sex: '{sample['sex']}'")

        if sample["sample_type"] == "proband":
            meta_paths = [
                path_in_accessdir(gens_accessdir, value)
                for value in (sample["meta_file"], sample["chrom_meta_file"])
                if value is not None
            ]
            if meta_paths:
                lines.append("    meta_files:")
                for path in meta_paths:
                    lines.append(f"      - '{path}'")

            annotation_rows: list[dict[str, Any]] = []
            if sample["roh_track"] is not None:
                annotation_rows.append(
                    {"file": path_in_accessdir(gens_accessdir, sample["roh_track"]), "name": "LOH"}
                )
            if trio and sample["upd_track"] is not None:
                annotation_rows.append(
                    {"file": path_in_accessdir(gens_accessdir, sample["upd_track"]), "name": "UPS"}
                )
            if annotation_rows:
                lines.append("    sample_annotations:")
                for annotation in annotation_rows:
                    lines.append(f"      - file: '{annotation['file']}'")
                    lines.append(f"        name: '{annotation['name']}'")

    return lines


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--gens-accessdir", required=True)
    parser.add_argument("--sample-id", nargs="+", dest="sample_ids", required=True)
    parser.add_argument("--sample-type", nargs="+", dest="sample_types")
    parser.add_argument("--sex", nargs="+", dest="sexes")
    parser.add_argument("--roh-track", nargs="+", dest="roh_tracks")
    parser.add_argument("--upd-track", nargs="+", dest="upd_tracks")
    parser.add_argument("--meta-file", nargs="+", dest="meta_files")
    parser.add_argument(
        "--chrom-meta-file", nargs="+", dest="chrom_meta_files"
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ped", type=Path)
    parser.add_argument("--trio", action="store_true")
    args = parser.parse_args()

    sample_ids = args.sample_ids
    sample_types = args.sample_types or []
    sexes = args.sexes or []
    roh_tracks = args.roh_tracks or []
    upd_tracks = args.upd_tracks or []
    meta_files = args.meta_files or []
    chrom_meta_files = args.chrom_meta_files or []

    samples = load_samples_from_cli_args(
        sample_ids=sample_ids,
        sample_types=sample_types,
        sexes=sexes,
        roh_tracks=roh_tracks,
        upd_tracks=upd_tracks,
        meta_files=meta_files,
        chrom_meta_files=chrom_meta_files,
    )

    ped_roles = load_ped_roles(args.ped) if args.ped is not None else None
    samples_with_ped_roles = apply_ped_roles(samples, ped_roles)

    selected = select_samples(samples_with_ped_roles, args.trio)
    yaml_lines = build_yaml_lines(
        case_id=args.case_id,
        gens_accessdir=args.gens_accessdir,
        selected_samples=selected,
        trio=args.trio,
    )
    args.output.write_text("\n".join(yaml_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
