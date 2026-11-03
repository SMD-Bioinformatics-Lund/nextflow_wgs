"""Optional native regression checks: install SVDB 2.8.2, Nextflow and GENMOD 3.9."""
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
##INFO=<ID=AC,Number=1,Type=Integer,Description="Count">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def run_annotation(tmp_path, source_info):
    if not shutil.which("svdb") or not shutil.which("nextflow"):
        pytest.skip("Native annotation check requires SVDB and Nextflow")
    # The only match for the second pair is a tiny deletion contained in the query.
    queries = []
    for pos in (10000, 50000):
        for kind, ref, alt in (("literal", "A" * 5001, "A"), ("symbolic", "A", "<DEL>")):
            queries.append(f"1\t{pos}\t{kind}_{pos}\t{ref}\t{alt}\t50\tPASS\tSVTYPE=DEL;SVLEN=-5000")
    queries.append("1\t10000\twrong_type\tN\t<DUP>\t50\tPASS\tSVTYPE=DUP;END=15000")
    (tmp_path / "query.vcf").write_text(HEADER + "\n".join(queries) + "\n")
    (tmp_path / "source.vcf").write_text(HEADER + "\n".join([
        f"1\t10000\tbest\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=14900;{source_info}",
        "1\t10000\tless_similar_common\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=14100;AF=0.5;AC=50",
        "1\t12000\tsmall\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=12126;AF=0.5;AC=50",
        "1\t52000\tonly_small\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=52126;AF=0.5;AC=50",
    ]) + "\n")
    (tmp_path / "test.nf").write_text(
        f"include {{ gnomad_sv }} from '{ROOT}/modules/gnomad_sv.nf'\n"
        "workflow { gnomad_sv(Channel.of(tuple('test', 'sample', file('query.vcf'))), file('source.vcf')) }\n"
    )
    (tmp_path / "test.config").write_text("""
params.container_svdb = ''
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


def test_literal_and_symbolic_queries_use_the_same_matching(tmp_path):
    result = run_annotation(tmp_path, "AF=0.001;AC=1")
    assert result.returncode == 0, result.stdout + result.stderr
    original = [line.split("\t") for line in (tmp_path / "query.vcf").read_text().splitlines() if not line.startswith("#")]
    annotated = [line.split("\t") for line in (tmp_path / "out/test.gnomad.vcf").read_text().splitlines() if not line.startswith("#")]
    assert len(annotated) == len(original)
    for before, after in zip(original, annotated):
        assert before[:7] == after[:7]
        assert after[7].startswith(before[7])
        if after[2].endswith("_10000"):
            # Best geometry wins, rather than the highest AF among passing matches.
            assert ";gnomad_svAF=0.001" in after[7]
            assert ";gnomad_svAC=1" in after[7]
        else:
            assert "gnomad_svAF=" not in after[7]


def test_old_slim_resource_missing_ac_fails_the_process(tmp_path):
    result = run_annotation(tmp_path, "AF=0.001")
    assert result.returncode != 0
    assert "Check source AC/AF" in result.stdout + result.stderr


@pytest.mark.parametrize("model", ["sv_rank_model_v5.3.ini", "sv_rank_model_panel_v1.1.ini"])
@pytest.mark.parametrize("frequency,expected", [(None, 6), (0.0001, 4), (0.001, 2), (0.01, 1), (0.02, -12)])
def test_rank_models_consume_info_frequency(model, frequency, expected):
    pytest.importorskip("genmod")
    from genmod.score_variants.config_parser import ConfigParser
    from genmod.score_variants.score_variant import get_plugin_score

    config = ConfigParser(str(ROOT / "rank_models/live/sv" / model))
    variant = {"info_dict": {} if frequency is None else {"gnomad_svAF": str(frequency)}}
    assert get_plugin_score(variant, "gnomad_sv", config)[0] == expected
