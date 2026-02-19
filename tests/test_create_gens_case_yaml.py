from pathlib import Path
import subprocess
import sys

import yaml

from bin.create_gens_case_yaml import main as create_case_yaml_main


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

    create_case_yaml_main(
        case_id=case_id,
        gens_accessdir="/gens/data",
        sample_ids=sample_ids,
        sample_types=sample_types,
        sexes=sexes,
        roh_tracks=roh_tracks,
        upd_tracks=upd_tracks,
        meta_files=meta_files,
        chrom_meta_files=chrom_meta_files,
        output=output_yaml,
        is_trio=trio,
    )
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

    case_data = yaml.safe_load(result)
    sample_ids = [sample["sample_id"] for sample in case_data["samples"]]
    assert sample_ids == ["kid", "mom", "dad"]

    proband_row = case_data["samples"][0]
    annotation_names = [annotation["name"] for annotation in proband_row["sample_annotations"]]
    assert "LOH" in annotation_names
    assert "UPS" in annotation_names
    assert "meta_files" in proband_row


def test_create_case_yaml_single(tmp_path: Path) -> None:
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

    case_data = yaml.safe_load(result)
    sample_ids = [sample["sample_id"] for sample in case_data["samples"]]
    assert sample_ids == ["kid", "mom", "dad"]

    proband_row = case_data["samples"][0]
    annotation_names = [annotation["name"] for annotation in proband_row["sample_annotations"]]
    assert "LOH" in annotation_names
    assert "UPS" not in annotation_names


def test_create_case_yaml_trio_no_non_mandatory(tmp_path: Path) -> None:
    result = run_create_case_yaml(
        tmp_path,
        case_id="CASE_TRIO_NO_OPTIONALS",
        sample_ids=["kid", "mom", "dad"],
        sample_types=["proband", "mother", "father"],
        sexes=["F", "F", "M"],
        roh_tracks=["none", "mom.roh.bed", "dad.roh.bed"],
        upd_tracks=["none", "mom.upd.bed", "dad.upd.bed"],
        meta_files=["none", "mom.meta.tsv", "dad.meta.tsv"],
        chrom_meta_files=["none", "mom.chrom.tsv", "dad.chrom.tsv"],
        trio=True,
    )

    case_data = yaml.safe_load(result)
    sample_ids = [sample["sample_id"] for sample in case_data["samples"]]
    assert sample_ids == ["kid", "mom", "dad"]

    proband_row = case_data["samples"][0]
    assert "meta_files" not in proband_row
    assert "sample_annotations" not in proband_row


def test_create_case_yaml_cli_smoke(tmp_path: Path) -> None:
    output_yaml = tmp_path / "CASE_CLI_SMOKE.gens_const.yaml"
    subprocess.run(
        [
            sys.executable,
            str(Path("bin/create_gens_case_yaml.py")),
            "--case_id",
            "CASE_CLI_SMOKE",
            "--gens_accessdir",
            "/gens/data",
            "--sample_ids",
            "kid",
            "--sample_types",
            "proband",
            "--sexes",
            "F",
            "--roh_tracks",
            "kid.roh.bed",
            "--upd_tracks",
            "kid.upd.bed",
            "--meta_files",
            "kid.meta.tsv",
            "--chrom_meta_files",
            "kid.chrom.tsv",
            "--output",
            str(output_yaml),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    result = output_yaml.read_text(encoding="utf-8")
    case_data = yaml.safe_load(result)
    assert case_data["case_id"] == "CASE_CLI_SMOKE"
    assert case_data["samples"][0]["sample_id"] == "kid"
