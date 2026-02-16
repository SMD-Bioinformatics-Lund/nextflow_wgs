from pathlib import Path
import subprocess
import sys


def run_create_case_yaml(
    tmp_path: Path,
    *,
    case_id: str,
    trio: bool,
    sample_ids: list[str],
    sample_types: list[str],
    sexes: list[str],
    roh_tracks: list[str],
    upd_tracks: list[str],
    meta_files: list[str],
    chrom_meta_files: list[str],
    ped_content: str | None = None,
) -> str:
    output_yaml = tmp_path / f"{case_id}.gens_const.yaml"

    cmd = [
        sys.executable,
        str(Path("bin/create_gens_case_yaml.py")),
        "--case-id",
        case_id,
        "--gens-accessdir",
        "/gens/data",
        "--sample-id",
        *sample_ids,
        "--sample-type",
        *sample_types,
        "--sex",
        *sexes,
        "--roh-track",
        *roh_tracks,
        "--upd-track",
        *upd_tracks,
        "--meta-file",
        *meta_files,
        "--chrom-meta-file",
        *chrom_meta_files,
        "--output",
        str(output_yaml),
    ]
    if trio:
        cmd.append("--trio")
    if ped_content is not None:
        ped = tmp_path / f"{case_id}_base.ped"
        ped.write_text(ped_content, encoding="utf-8")
        cmd.extend(["--ped", str(ped)])

    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return output_yaml.read_text(encoding="utf-8")


def test_create_case_yaml_trio_has_proband_mother_father_and_annotations(
    tmp_path: Path,
) -> None:
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_TRIO",
        trio=True,
        sample_ids=["dad", "kid", "mom"],
        sample_types=["father", "proband", "mother"],
        sexes=["M", "F", "F"],
        roh_tracks=["dad.roh.bed", "kid.roh.bed", "mom.roh.bed"],
        upd_tracks=["dad.upd.bed", "kid.upd.bed", "mom.upd.bed"],
        meta_files=["dad.meta.tsv", "kid.meta.tsv", "mom.meta.tsv"],
        chrom_meta_files=["dad.chrom.tsv", "kid.chrom.tsv", "mom.chrom.tsv"],
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
        trio=False,
        sample_ids=["mom", "kid", "dad"],
        sample_types=["mother", "proband", "father"],
        sexes=["F", "M", "M"],
        roh_tracks=["mom.roh.bed", "kid.roh.bed", "dad.roh.bed"],
        upd_tracks=["mom.upd.bed", "kid.upd.bed", "dad.upd.bed"],
        meta_files=["mom.meta.tsv", "kid.meta.tsv", "dad.meta.tsv"],
        chrom_meta_files=["mom.chrom.tsv", "kid.chrom.tsv", "dad.chrom.tsv"],
    )

    assert result.count("- sample_id:") == 1
    assert "- sample_id: 'kid'" in result
    assert "name: 'LOH'" in result
    assert "name: 'UPS'" not in result


def test_create_case_yaml_trio_uses_ped_roles_when_sample_types_are_relative(
    tmp_path: Path,
) -> None:
    ped = (
        "CASE_REL\tkid\tdad\tmom\t2\t2\n"
        "CASE_REL\tdad\t0\t0\t1\t1\n"
        "CASE_REL\tmom\t0\t0\t2\t1\n"
    )
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_REL",
        trio=True,
        sample_ids=["dad", "kid", "mom"],
        sample_types=["relative", "relative", "relative"],
        sexes=["M", "F", "F"],
        roh_tracks=["dad.roh.bed", "kid.roh.bed", "mom.roh.bed"],
        upd_tracks=["dad.upd.bed", "kid.upd.bed", "mom.upd.bed"],
        meta_files=["dad.meta.tsv", "kid.meta.tsv", "mom.meta.tsv"],
        chrom_meta_files=["dad.chrom.tsv", "kid.chrom.tsv", "mom.chrom.tsv"],
        ped_content=ped,
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
    assert "sample_type: 'proband'" in result
    assert "sample_type: 'mother'" in result
    assert "sample_type: 'father'" in result


def test_create_case_yaml_single_uses_ped_to_mark_proband(tmp_path: Path) -> None:
    ped = "CASE_ONE\tsolo\t0\t0\t1\t2\n"
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_ONE",
        trio=False,
        sample_ids=["solo"],
        sample_types=["relative"],
        sexes=["M"],
        roh_tracks=["solo.roh.bed"],
        upd_tracks=["solo.upd.bed"],
        meta_files=["solo.meta.tsv"],
        chrom_meta_files=["solo.chrom.tsv"],
        ped_content=ped,
    )

    assert "- sample_id: 'solo'" in result
    assert "sample_type: 'proband'" in result
    assert "name: 'LOH'" in result
    assert "name: 'UPS'" not in result
