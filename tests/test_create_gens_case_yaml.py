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


def test_create_case_yaml_trio(
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


def test_create_case_yaml_single(tmp_path: Path) -> None:
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_SINGLE",
        sample_ids=["kid"],
        sample_types=["proband"],
        sexes=["F"],
        roh_tracks=["kid.roh.bed"],
        upd_tracks=[],
        meta_files=["kid.meta.tsv"],
        chrom_meta_files=["kid.chrom.tsv"],
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


def test_create_case_yaml_trio_no_non_mandatory(tmp_path: Path) -> None:
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_TRIO_NO_OPTIONALS",
        sample_ids=["kid", "mom", "dad"],
        sample_types=[],
        sexes=[],
        roh_tracks=[],
        upd_tracks=[],
        meta_files=[],
        chrom_meta_files=[],
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
    assert "meta_files:" not in result
    assert "sample_annotations:" not in result
    assert "name: 'LOH'" not in result
    assert "name: 'UPS'" not in result
