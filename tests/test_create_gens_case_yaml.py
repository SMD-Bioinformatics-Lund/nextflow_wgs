from pathlib import Path
import subprocess
import sys


def run_create_case_yaml(
    tmp_path: Path, sample_table_content: str, *, case_id: str, trio: bool
) -> str:
    sample_table = tmp_path / "samples.tsv"
    output_yaml = tmp_path / f"{case_id}.gens_const.yaml"
    sample_table.write_text(sample_table_content, encoding="utf-8")

    cmd = [
        sys.executable,
        str(Path("bin/create_gens_case_yaml.py")),
        "--case-id",
        case_id,
        "--gens-accessdir",
        "/gens/data",
        "--sample-table",
        str(sample_table),
        "--output",
        str(output_yaml),
    ]
    if trio:
        cmd.append("--trio")

    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return output_yaml.read_text(encoding="utf-8")


def run_create_case_yaml_grouped_args(tmp_path: Path, *, case_id: str, trio: bool) -> str:
    output_yaml = tmp_path / f"{case_id}.gens_const.yaml"
    cmd = [
        sys.executable,
        str(Path("bin/create_gens_case_yaml.py")),
        "--case-id",
        case_id,
        "--gens-accessdir",
        "/gens/data",
        "--sample-id",
        "kid",
        "mom",
        "dad",
        "--sample-type",
        "proband",
        "mother",
        "father",
        "--sex",
        "F",
        "F",
        "M",
        "--roh-track",
        "kid.roh.bed",
        "mom.roh.bed",
        "dad.roh.bed",
        "--upd-track",
        "kid.upd.bed",
        "mom.upd.bed",
        "dad.upd.bed",
        "--meta-file",
        "kid.meta.tsv",
        "mom.meta.tsv",
        "dad.meta.tsv",
        "--chrom-meta-file",
        "kid.chrom.tsv",
        "mom.chrom.tsv",
        "dad.chrom.tsv",
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
    sample_table = (
        "sample_id\tsample_type\tsex\troh_track\tupd_track\tmeta_file\tchrom_meta_file\n"
        "dad\tfather\tM\tdad.roh.bed\tdad.upd.bed\tdad.meta.tsv\tdad.chrom.tsv\n"
        "kid\tproband\tF\tkid.roh.bed\tkid.upd.bed\tkid.meta.tsv\tkid.chrom.tsv\n"
        "mom\tmother\tF\tmom.roh.bed\tmom.upd.bed\tmom.meta.tsv\tmom.chrom.tsv\n"
    )

    result = run_create_case_yaml(
        tmp_path, sample_table, case_id="CASE_TRIO", trio=True
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
    sample_table = (
        "sample_id\tsample_type\tsex\troh_track\tupd_track\tmeta_file\tchrom_meta_file\n"
        "mom\tmother\tF\tmom.roh.bed\tmom.upd.bed\tmom.meta.tsv\tmom.chrom.tsv\n"
        "kid\tproband\tM\tkid.roh.bed\tkid.upd.bed\tkid.meta.tsv\tkid.chrom.tsv\n"
        "dad\tfather\tM\tdad.roh.bed\tdad.upd.bed\tdad.meta.tsv\tdad.chrom.tsv\n"
    )

    result = run_create_case_yaml(
        tmp_path, sample_table, case_id="CASE_SINGLE", trio=False
    )

    assert result.count("- sample_id:") == 1
    assert "- sample_id: 'kid'" in result
    assert "name: 'LOH'" in result
    assert "name: 'UPS'" not in result


def test_create_case_yaml_single_falls_back_to_ordered_first_sample_without_proband(
    tmp_path: Path,
) -> None:
    sample_table = (
        "sample_id\tsample_type\tsex\troh_track\tupd_track\tmeta_file\tchrom_meta_file\n"
        "alpha\trelative\tF\talpha.roh.bed\talpha.upd.bed\talpha.meta.tsv\talpha.chrom.tsv\n"
        "beta\tmother\tF\tbeta.roh.bed\tbeta.upd.bed\tbeta.meta.tsv\tbeta.chrom.tsv\n"
    )

    result = run_create_case_yaml(
        tmp_path, sample_table, case_id="CASE_FALLBACK", trio=False
    )

    assert result.count("- sample_id:") == 1
    assert "- sample_id: 'beta'" in result
    assert "sample_annotations:" not in result


def test_create_case_yaml_trio_from_grouped_cli_args(tmp_path: Path) -> None:
    result = run_create_case_yaml_grouped_args(
        tmp_path, case_id="CASE_GROUPED", trio=True
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
