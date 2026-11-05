"""Optional native SVDB 2.12.2/Nextflow and rank model regression checks."""

import gzip
import os
from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
HEADER = """##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##INFO=<ID=AF,Number=1,Type=Float,Description="Frequency">
##INFO=<ID=N_HET,Number=1,Type=Integer,Description="Heterozygote count">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def run_annotation(tmp_path, source, query):
    if not shutil.which("svdb") or not shutil.which("nextflow"):
        pytest.skip("Native annotation check requires SVDB 2.12.2 and Nextflow")
    (tmp_path / "query.vcf").write_text(HEADER + "\n".join(query) + "\n")
    with gzip.open(tmp_path / "source.vcf.gz", "wt") as output:
        output.write(HEADER + "\n".join(source) + "\n")
    (tmp_path / "test.nf").write_text(
        f"include {{ gnomad_sv }} from '{ROOT}/modules/gnomad_sv.nf'\n"
        "workflow { gnomad_sv(Channel.of(tuple('test', 'sample', file('query.vcf'))), file('source.vcf.gz')) }\n"
    )
    (tmp_path / "test.config").write_text("""
params.container_gnomad_sv = ''
docker.enabled = false
singularity.enabled = false
process.executor = 'local'
process.cpus = 1
process.memory = '1 GB'
process.publishDir = [path: 'out', mode: 'copy']
""")
    return subprocess.run(
        ["nextflow", "run", "test.nf", "-c", "test.config", "-ansi-log", "false"],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
        env={**os.environ, "NXF_OFFLINE": "true"},
    )


def info_fields(record):
    return dict(field.split("=", 1) for field in record[7].split(";") if "=" in field)


def test_native_matching_keeps_af_with_zero_heterozygotes(tmp_path):
    query = []
    source = []
    for pos, count in [(10000, 5), (50000, 0)]:
        for kind, ref, alt in [("literal", "A" * 5001, "A"), ("symbolic", "A", "<DEL>")]:
            query.append(f"1\t{pos}\t{kind}_{pos}\t{ref}\t{alt}\t50\tPASS\tSVTYPE=DEL;SVLEN=-5000")
        source.extend([
            f"1\t{pos}\tbest\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END={pos + 4900};AF=0.02;N_HET={count}",
            f"1\t{pos}\tother\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END={pos + 4100};AF=0.5;N_HET=100",
        ])
    query.extend([
        "1\t10000\twrong_type\tN\t<DUP>\t50\tPASS\tSVTYPE=DUP;END=15000",
        "1\t90000\tunmatched\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=95000",
        "1\t130000\tonly_small\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=135000",
        "1\t170000\tins_match\tN\t<INS>\t50\tPASS\tSVTYPE=INS;SVLEN=100",
        "1\t210000\tins_wrong_length\tN\t<INS>\t50\tPASS\tSVTYPE=INS;SVLEN=1000",
    ])
    source.extend([
        "1\t132000\tsmall\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=132126;AF=0.5;N_HET=100",
        "1\t170020\tins\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=90;AF=0.03;N_HET=0",
        "1\t210020\tins\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=90;AF=0.03;N_HET=10",
    ])
    result = run_annotation(tmp_path, source, query)
    assert result.returncode == 0, result.stdout + result.stderr
    annotated = [
        line.split("\t") for line in (tmp_path / "out/test.gnomad.vcf").read_text().splitlines()
        if not line.startswith("#")
    ]
    assert len(annotated) == len(query)
    for original, record in zip(query, annotated):
        assert record[:7] == original.split("\t")[:7]
        info = info_fields(record)
        assert "gnomad_svAC" not in info
        if record[2].endswith(("_10000", "_50000")):
            assert float(info["gnomad_svAF"]) == 0.02
            assert int(info["gnomad_svN_HET"]) == (0 if record[2].endswith("_50000") else 5)
        elif record[2] == "ins_match":
            assert float(info["gnomad_svAF"]) == 0.03
            assert int(info["gnomad_svN_HET"]) == 0
        else:
            assert "gnomad_svAF" not in info
    assert "2.12.2" in (tmp_path / "out/gnomad_sv_versions.yml").read_text()


@pytest.mark.parametrize("model", ["sv_rank_model_v5.3.ini", "sv_rank_model_panel_v1.1.ini"])
@pytest.mark.parametrize("frequency,expected", [(None, 6), (0.0001, 4), (0.001, 2), (0.01, 1), (0.02, -12)])
def test_rank_models_consume_info_frequency(model, frequency, expected):
    pytest.importorskip("genmod")
    from genmod.score_variants.config_parser import ConfigParser
    from genmod.score_variants.score_variant import get_plugin_score

    config = ConfigParser(str(ROOT / "rank_models/live/sv" / model))
    variant = {"info_dict": {} if frequency is None else {"gnomad_svAF": str(frequency)}}
    assert get_plugin_score(variant, "gnomad_sv", config)[0] == expected
