import gzip
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "bin/gnomad_del.py"
FIELDS = ["gnomad", "gnomad_AF", "gnomad_FREQ_HOMALT", "gnomad_SVTYPE"]
SCHEMA = ["Allele", "Consequence", "SYMBOL", "VARIANT_CLASS"] + FIELDS


def variant(pos=1000, length=500, kind="DEL", symbolic=False, name="q", csq=True):
    ref, alt = ("A" + "T" * length, "A") if kind == "DEL" else ("A", "A" + "T" * length)
    if symbolic:
        ref, alt = "N", "<{}>".format(kind)
    info = "END={};SVLEN={};SVTYPE={};FLAG".format(
        pos + length if kind == "DEL" else pos, -length if kind == "DEL" else length, kind)
    if csq:
        cls = "deletion" if kind == "DEL" else "insertion"
        info += ";CSQ=-|intron_variant|GENE1|{}|old|0.9|0.8|{},-|splice_region_variant|GENE2|{}|old|0.9|0.8|{}".format(cls, kind, cls, kind)
    return ["1", str(pos), name, ref, alt, "100", "PASS", info, "GT:DP", "0/1:24"]


def reference(pos=1000, length=500, kind="DEL", name="r", af="0.03", hom="0.001"):
    return ["chr1", str(pos), name, "N", "<{}>".format(kind), ".", "PASS",
            "END={};SVLEN={};SVTYPE={};AF={};FREQ_HOMALT={}".format(pos + length, length, kind, af, hom)]


def write_vcf(path, rows, csq=False):
    header = "##fileformat=VCFv4.2\n"
    if csq:
        header += '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {}">\n'.format("|".join(SCHEMA))
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    with (gzip.open(path, "wt") if str(path).endswith(".gz") else path.open("w")) as handle:
        handle.write(header)
        handle.writelines("\t".join(row) + "\n" for row in rows)


@pytest.fixture
def annotate(tmp_path):
    def run(queries, references, csq=True):
        vcf, source, bed, intervals, matches, output = (
            tmp_path / name for name in ("input.vcf", "reference.vcf.gz", "reference.bed.gz",
                                        "intervals.vcf", "matches.vcf", "output.vcf"))
        write_vcf(vcf, queries, csq)
        write_vcf(source, references)

        def helper(mode, input_path, output_path, *extra):
            subprocess.run([sys.executable, str(SCRIPT), mode, "--vcf", str(input_path),
                            "--output", str(output_path), *extra], check=True, capture_output=True, text=True)

        helper("reference", source, bed)
        helper("query", vcf, intervals)
        subprocess.run(["bcftools", "annotate", "--no-version", "-a", str(bed),
                        "-c", "CHROM,FROM,TO,GNOMAD_DEL_MATCHES", "-l", "GNOMAD_DEL_MATCHES:append",
                        "--min-overlap", "0.699999:0.699999", "-Ov", "-o", str(matches), str(intervals)],
                       check=True, capture_output=True, text=True)
        helper("restore", vcf, output, "--matches", str(matches))
        lines = output.read_text().splitlines()
        rows = [line.split("\t") for line in lines if not line.startswith("#")]
        annotations = []
        for row in rows:
            info = dict(entry.split("=", 1) for entry in row[7].split(";") if "=" in entry)
            annotations.append([dict(zip(SCHEMA, entry.split("|"))) for entry in info["CSQ"].split(",")]
                               if csq and "CSQ" in info else [info])
        return rows, annotations, lines
    return run


@pytest.mark.parametrize("symbolic", [False, True])
@pytest.mark.parametrize("length", [79, 500, 999, 1000, 1001])
def test_symbolic_and_sequence_deletions_match_and_keep_original_data(annotate, symbolic, length):
    original = variant(length=length, symbolic=symbolic)
    rows, annotations, lines = annotate([original], [reference(length=length)])
    assert rows[0][:7] == original[:7]
    assert rows[0][8:] == original[8:]
    assert ";FLAG;" in rows[0][7]
    assert len([line for line in lines if line.startswith("##INFO=<ID=CSQ,")]) == 1
    for entry in annotations[0]:
        assert [entry[field] for field in FIELDS] == ["r", "0.03", "0.001", "DEL"]


@pytest.mark.parametrize("pos, length, expected", [
    (1150, 500, "r"), (1151, 500, ""), (1000, 800, ""), (1000, 300, ""),
])
def test_reciprocal_overlap_inclusive_at_seventy_percent(annotate, pos, length, expected):
    _, annotations, _ = annotate([variant()], [reference(pos, length)])
    assert annotations[0][0]["gnomad"] == expected


@pytest.mark.parametrize("length, shift, expected", [
    (50000000, 15000000, "r"), (50000000, 15000001, ""), (49999990, 14999997, "r"),
])
def test_exact_boundary_despite_bcftools_float_rounding(annotate, length, shift, expected):
    _, annotations, _ = annotate([variant(length=length, symbolic=True)], [reference(1000 + shift, length)])
    assert annotations[0][0]["gnomad"] == expected


def test_only_ordinary_deletion_class_is_matched_and_all_frequencies_are_kept(annotate):
    _, annotations, _ = annotate([variant()], [
        reference(name="exact", af="0.0001"), reference(1150, name="shifted", af="0.03"),
        reference(kind="DUP", name="dup", af="0.8"), reference(kind="DEL:ME:LINE1", name="mobile", af="0.9")])
    assert annotations[0][0]["gnomad"] == "exact&shifted"
    assert annotations[0][0]["gnomad_AF"] == "0.0001&0.03"


@pytest.mark.parametrize("symbolic", [False, True])
def test_insertions_are_unchanged_including_preexisting_vep_matches(annotate, symbolic):
    original = variant(kind="INS", symbolic=symbolic)
    rows, _, _ = annotate([original], [reference()])
    assert rows == [original]


def test_manta_mixed_deletion_with_unknown_reference_uses_svlen(annotate):
    original = variant(length=201, symbolic=True)
    original[3:5] = ["A" + "N" * 201, "AGA"]
    original[7] += ";CIGAR=1M2I201D"
    rows, annotations, _ = annotate([original], [reference(length=201)])
    assert rows[0][3:5] == original[3:5]
    assert annotations[0][0]["gnomad"] == "r"


def test_first_fifty_matches_follow_reference_source_order(annotate):
    refs = [reference(1001, name="r{}".format(i), af="0.0001") for i in range(50)]
    refs.append(reference(1000, name="excluded", af="0.9"))
    _, annotations, _ = annotate([variant()], refs)
    assert annotations[0][0]["gnomad"].split("&") == ["r{}".format(i) for i in range(50)] + ["..."]
    assert max(map(float, annotations[0][0]["gnomad_AF"].split("&"))) == 0.0001


def test_missing_values_and_vep_field_formatting(annotate):
    refs = [reference(name="first;second", af="0", hom="."), reference(name="r2", af="0.97,0.03")]
    refs[1][7] = refs[1][7].replace(";FREQ_HOMALT=0.001", "")
    _, annotations, _ = annotate([variant()], refs)
    assert [annotations[0][0][field] for field in FIELDS] == ["first&r2", "0&0.03", ".", "DEL&DEL"]


def test_duplicate_query_ids_and_records_without_csq_are_not_confused(annotate):
    queries = [variant(name="."), variant(10000, name="."), variant(1000, name="skipped", csq=False)]
    rows, annotations, _ = annotate(queries, [reference(10000)])
    assert annotations[0][0]["gnomad"] == ""
    assert annotations[1][0]["gnomad"] == "r"
    assert rows[2] == queries[2]


def test_padding_changes_only_the_matching_interval(annotate):
    original = variant(length=3)
    original[3:5] = ["AATTT", "AA"]
    rows, annotations, _ = annotate([original], [reference(1001, length=3)])
    assert rows[0][1:5] == original[1:5]
    assert annotations[0][0]["gnomad"] == "r"


def test_unannotated_benchmark_vcf_gets_info_fields(annotate):
    original = variant(csq=False)
    rows, annotations, _ = annotate([original], [reference()], csq=False)
    assert rows[0][:7] == original[:7]
    assert annotations[0][0]["gnomad_AF"] == "0.03"


def test_raw_benchmark_allele_without_shared_anchor_requires_vep_classification(annotate):
    original = variant(length=50, csq=False)
    original[4] = "C"
    rows, _, _ = annotate([original], [reference(length=50)], csq=False)
    assert rows == [original]


def test_empty_reference_and_empty_query(annotate):
    _, annotations, _ = annotate([variant()], [])
    assert annotations[0][0]["gnomad"] == ""
    rows, annotations, _ = annotate([], [reference()])
    assert rows == annotations == []
