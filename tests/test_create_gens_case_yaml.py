from pathlib import Path
import subprocess
import sys


def run_create_case_yaml(
    tmp_path: Path,
    *,
    case_id: str,
    sample_ids: list[str],
    sample_types: list[str],
    sexes: list[str],
    roh_tracks: list[str],
    upd_tracks: list[str],
    meta_files: list[str],
    chrom_meta_files: list[str],
    trio: bool,
) -> str:
    output_yaml = tmp_path / f"{case_id}.gens_const.yaml"

    cmd = [
        sys.executable,
        str(Path("bin/create_gens_case_yaml.py")),
        "--case_id",
        case_id,
        "--gens_accessdir",
        "/gens/data",
        "--sample_ids",
        *sample_ids,
        "--sample_types",
        *sample_types,
        "--sexes",
        *sexes,
        "--roh_tracks",
        *roh_tracks,
        "--upd_tracks",
        *upd_tracks,
        "--meta_files",
        *meta_files,
        "--chrom_meta_files",
        *chrom_meta_files,
        "--output",
        str(output_yaml),
    ]
    if trio:
        cmd.append("--trio")

    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return output_yaml.read_text(encoding="utf-8")


def test_create_case_yaml_trio_has_proband_mother_father_and_annotations(
    tmp_path: Path,
) -> None:
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_TRIO",
        sample_ids=["dad", "kid", "mom"],
        sample_types=["father", "proband", "mother"],
        sexes=["M", "F", "F"],
        roh_tracks=["dad.roh.bed", "kid.roh.bed", "mom.roh.bed"],
        upd_tracks=["dad.upd.bed", "kid.upd.bed", "mom.upd.bed"],
        meta_files=["dad.meta.tsv", "kid.meta.tsv", "mom.meta.tsv"],
        chrom_meta_files=["dad.chrom.tsv", "kid.chrom.tsv", "mom.chrom.tsv"],
        trio=True,
    )

    sample_lines = [
        line.strip()
        for line in result.splitlines()
        if line.strip().startswith("- sample_id:")
    ]
    assert sample_lines == [
        "- sample_id: 'kid'",
        "- sample_id: 'mom'",
        "- sample_id: 'dad'",
    ]
    assert "name: 'LOH'" in result
    assert "name: 'UPS'" in result
    assert "meta_files:" in result


def test_create_case_yaml_single_only_keeps_proband_and_no_ups(tmp_path: Path) -> None:
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_SINGLE",
        sample_ids=["mom", "kid", "dad"],
        sample_types=["mother", "proband", "father"],
        sexes=["F", "M", "M"],
        roh_tracks=["mom.roh.bed", "kid.roh.bed", "dad.roh.bed"],
        upd_tracks=["mom.upd.bed", "kid.upd.bed", "dad.upd.bed"],
        meta_files=["mom.meta.tsv", "kid.meta.tsv", "dad.meta.tsv"],
        chrom_meta_files=["mom.chrom.tsv", "kid.chrom.tsv", "dad.chrom.tsv"],
        trio=False,
    )

    sample_lines = [
        line.strip()
        for line in result.splitlines()
        if line.strip().startswith("- sample_id:")
    ]
    assert sample_lines == [
        "- sample_id: 'kid'",
        "- sample_id: 'mom'",
        "- sample_id: 'dad'",
    ]
    assert "name: 'LOH'" in result
    assert "name: 'UPS'" not in result


def test_create_case_yaml_single_falls_back_to_ordered_first_sample_without_proband(
    tmp_path: Path,
) -> None:
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_FALLBACK",
        sample_ids=["alpha", "beta"],
        sample_types=["relative", "mother"],
        sexes=["F", "F"],
        roh_tracks=["alpha.roh.bed", "beta.roh.bed"],
        upd_tracks=["alpha.upd.bed", "beta.upd.bed"],
        meta_files=["alpha.meta.tsv", "beta.meta.tsv"],
        chrom_meta_files=["alpha.chrom.tsv", "beta.chrom.tsv"],
        trio=False,
    )

    sample_lines = [
        line.strip()
        for line in result.splitlines()
        if line.strip().startswith("- sample_id:")
    ]
    assert sample_lines == [
        "- sample_id: 'beta'",
        "- sample_id: 'alpha'",
    ]
    assert "sample_annotations:" not in result


def test_create_case_yaml_trio_from_grouped_cli_args(tmp_path: Path) -> None:
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_GROUPED",
        sample_ids=["kid", "mom", "dad"],
        sample_types=["proband", "mother", "father"],
        sexes=["F", "F", "M"],
        roh_tracks=["kid.roh.bed", "mom.roh.bed", "dad.roh.bed"],
        upd_tracks=["kid.upd.bed", "mom.upd.bed", "dad.upd.bed"],
        meta_files=["kid.meta.tsv", "mom.meta.tsv", "dad.meta.tsv"],
        chrom_meta_files=["kid.chrom.tsv", "mom.chrom.tsv", "dad.chrom.tsv"],
        trio=True,
    )

    sample_lines = [
        line.strip()
        for line in result.splitlines()
        if line.strip().startswith("- sample_id:")
    ]
    assert sample_lines == [
        "- sample_id: 'kid'",
        "- sample_id: 'mom'",
        "- sample_id: 'dad'",
    ]
    assert "name: 'LOH'" in result
    assert "name: 'UPS'" in result
