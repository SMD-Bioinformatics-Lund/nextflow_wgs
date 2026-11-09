import gzip
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "bin/annotate_gnomad_sv_bedtools.py"
SCHEMA = ["Allele", "Consequence", "SYMBOL", "VARIANT_CLASS"]
GNOMAD_FIELDS = ["gnomad", "gnomad_AF", "gnomad_FREQ_HOMALT", "gnomad_SVTYPE"]


def query(pos=1000, length=500, *, symbolic=False, kind="DEL", name="q", csq=True):
    variant_class = {"DEL": "deletion", "DUP": "duplication", "INS": "insertion",
                     "INV": "inversion", "DUP:TANDEM": "tandem_duplication"}[kind]
    ref, alt = "N", "<{}>".format(kind)
    if not symbolic and kind in ("DEL", "INS"):
        ref, alt = ("A" + "T" * length, "A") if kind == "DEL" else ("A", "A" + "T" * length)
    end = pos if "INS" in kind else pos + length
    info = "END={};SVLEN={};SVTYPE={};FLAG".format(end, -length if kind == "DEL" else length, kind)
    if csq:
        info += ";CSQ=-|intron_variant|GENE1|{},-|splice_region_variant|GENE2|{}".format(variant_class, variant_class)
    return ["1", str(pos), name, ref, alt, "100", "PASS", info, "GT:DP", "0/1:24"]


def reference(pos=1000, length=500, *, kind="DEL", af="0.03", hom="0.001", name="r", chrom="chr1", filter_value="PASS"):
    return [chrom, str(pos), name, "N", "<{}>".format(kind), ".", filter_value,
            "END={};SVLEN={};SVTYPE={};AF={};FREQ_HOMALT={}".format(pos + length, length, kind, af, hom)]


def write_vcf(path, variants, schema=None, extra_headers=()):
    text = "##fileformat=VCFv4.2\n"
    if schema:
        text += '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {}">\n'.format("|".join(schema))
    text += "".join(line + "\n" for line in extra_headers)
    text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    text += "".join("\t".join(variant) + "\n" for variant in variants)
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


@pytest.fixture
def annotate(tmp_path):
    def run(queries, references, schema=None, extra_headers=()):
        source, database, synonyms, output = (tmp_path / name for name in (
            "query.vcf", "reference.vcf.gz", "synonyms", "annotated.vcf"))
        write_vcf(source, queries, schema or SCHEMA, extra_headers)
        write_vcf(database, references)
        synonyms.write_text("1\tchr1\tNC_000001.11\nMT\tchrM\n")
        subprocess.run([sys.executable, str(SCRIPT), "--vcf", str(source),
                        "--database", str(database), "--synonyms", str(synonyms),
                        "--output", str(output)], check=True, capture_output=True, text=True)
        lines = output.read_text().splitlines()
        header = next(line for line in lines if line.startswith("##INFO=<ID=CSQ,"))
        fields = header.split("Format: ")[1].split('"')[0].split("|")
        rows = [line.split("\t") for line in lines if not line.startswith("#")]
        annotations = []
        for row in rows:
            csq = next((field[4:] for field in row[7].split(";") if field.startswith("CSQ=")), None)
            annotations.append([dict(zip(fields, entry.split("|"))) for entry in csq.split(",")] if csq else [])
        return rows, annotations, fields
    return run


@pytest.mark.parametrize("symbolic", [False, True])
@pytest.mark.parametrize("length", [79, 500, 999, 1000, 1001])
def test_exact_deletions_and_original_vcf_are_retained(annotate, symbolic, length):
    original = query(length=length, symbolic=symbolic)
    rows, annotations, _ = annotate([original], [reference(length=length)])
    assert rows[0][:7] == original[:7]
    assert rows[0][8:] == original[8:]
    assert ";FLAG;" in rows[0][7]
    for entry in annotations[0]:
        assert [entry[field] for field in GNOMAD_FIELDS] == ["r", "0.03", "0.001", "DEL"]
    assert [(entry["Consequence"], entry["SYMBOL"]) for entry in annotations[0]] == [
        ("intron_variant", "GENE1"), ("splice_region_variant", "GENE2")]


@pytest.mark.parametrize("pos, length, matches", [
    (1100, 500, True), (1150, 500, True), (1151, 500, False),
    (1200, 500, False), (1000, 800, False), (2000, 500, False),
])
def test_true_reciprocal_overlap_is_inclusive_at_seventy_percent(annotate, pos, length, matches):
    _, annotations, _ = annotate([query()], [reference(pos, length)])
    assert annotations[0][0]["gnomad_AF"] == ("0.03" if matches else "")


def test_small_reference_inside_larger_query_is_also_rejected(annotate):
    _, annotations, _ = annotate([query(length=800)], [reference(length=500)])
    assert annotations[0][0]["gnomad_AF"] == ""


@pytest.mark.parametrize("length, shift, matches", [
    (50000000, 15000000, True), (50000000, 15000001, False),
    (49999990, 14999997, True),
])
def test_large_sv_threshold_is_exact_despite_bedtools_float_rounding(annotate, length, shift, matches):
    _, annotations, _ = annotate([query(length=length, symbolic=True)], [reference(1000 + shift, length)])
    assert annotations[0][0]["gnomad_AF"] == ("0.03" if matches else "")


def test_maximum_frequency_comes_from_all_reported_matches(annotate):
    refs = [reference(name="exact_rare", af="0.0001"),
            reference(1150, name="shifted_common", af="0.03", filter_value="LowQual"),
            reference(1200, name="invalid_overlap", af="0.9"),
            reference(kind="DUP", name="invalid_type", af="0.8")]
    _, annotations, _ = annotate([query()], refs)
    assert annotations[0][0]["gnomad"] == "exact_rare&shifted_common"
    assert max(map(float, annotations[0][0]["gnomad_AF"].split("&"))) == 0.03


def test_same_vep_class_including_tandem_duplications(annotate):
    _, annotations, _ = annotate([query(symbolic=True, kind="DUP:TANDEM")], [
        reference(kind="DUP", name="wrong_class"), reference(kind="DUP:TANDEM", name="right_class")])
    assert annotations[0][0]["gnomad"] == "right_class"


@pytest.mark.parametrize("symbolic", [False, True])
def test_sequence_and_symbolic_insertions_use_the_same_length_interval(annotate, symbolic):
    _, annotations, _ = annotate([query(length=79, symbolic=symbolic, kind="INS")], [
        reference(length=79, kind="INS", name="same"), reference(length=200, kind="INS", name="too_long")])
    assert annotations[0][0]["gnomad"] == "same"


def test_record_limit_and_source_order_match_vep_default(annotate):
    refs = [reference(name="r{}".format(index), af="0.0001") for index in range(50)]
    refs += [reference(name="over_limit", af="0.9")]
    _, annotations, _ = annotate([query()], refs)
    assert annotations[0][0]["gnomad"].split("&") == ["r{}".format(index) for index in range(50)] + ["..."]
    assert len(annotations[0][0]["gnomad_AF"].split("&")) == 50
    assert max(map(float, annotations[0][0]["gnomad_AF"].split("&"))) == 0.0001


def test_records_without_csq_are_not_rescued(annotate):
    original = query(csq=False)
    rows, annotations, _ = annotate([original], [reference()])
    assert rows == [original]
    assert annotations == [[]]


def test_duplicate_ids_do_not_mix_variants(annotate):
    queries = [query(name="."), query(10000, name=".")]
    _, annotations, _ = annotate(queries, [reference(10000)])
    assert annotations[0][0]["gnomad_AF"] == ""
    assert annotations[1][0]["gnomad_AF"] == "0.03"


# These output expectations were checked against the unmodified VEP 114.2
# AnnotationSource/File/VCF.pm and OutputFactory/VCF.pm methods.
def test_zero_af_and_literal_missing_values_are_retained(annotate):
    _, annotations, _ = annotate([query()], [reference(af="0", hom=".")])
    assert annotations[0][0]["gnomad_AF"] == "0"
    assert annotations[0][0]["gnomad_FREQ_HOMALT"] == "."


def test_absent_fields_are_omitted_from_lists_as_in_vep(annotate):
    refs = [reference(name="r1", hom="."), reference(name="r2", af="0.02")]
    refs[0][7] = refs[0][7].replace(";AF=0.03", "").replace(";SVTYPE=DEL", "")
    refs[1][7] = refs[1][7].replace(";FREQ_HOMALT=0.001", "")
    _, annotations, _ = annotate([query(symbolic=True)], refs)
    assert [annotations[0][0][field] for field in GNOMAD_FIELDS] == ["r1&r2", "0.02", ".", "DEL"]


@pytest.mark.parametrize("af, hom, expected_af, expected_hom", [
    ("0.97,0.03", "0.999,0.001", "0.03", "0.001"),
    ("0.01,0.02,0.03", ".", "0.01&0.02&0.03", "."),
])
def test_reference_field_lists_follow_vep_allele_handling(annotate, af, hom, expected_af, expected_hom):
    _, annotations, _ = annotate([query(symbolic=True)], [reference(af=af, hom=hom)])
    assert annotations[0][0]["gnomad_AF"] == expected_af
    assert annotations[0][0]["gnomad_FREQ_HOMALT"] == expected_hom


def test_reference_info_flags_are_reported_as_one(annotate):
    ref = reference()
    ref[7] = ref[7].replace("AF=0.03", "AF")
    _, annotations, _ = annotate([query(symbolic=True)], [ref])
    assert annotations[0][0]["gnomad_AF"] == "1"


@pytest.mark.parametrize("name, expected", [
    ("first;second", "first"), (".", "chr1:1001-1500"), ("0", "chr1:1001-1500"),
    ("id|one,two with space", "id&one&two_with_space"), ("-", ""),
])
def test_reference_names_follow_vep_output(annotate, name, expected):
    _, annotations, _ = annotate([query(symbolic=True)], [reference(name=name)])
    assert annotations[0][0]["gnomad"] == expected


def test_custom_headers_match_vep_and_other_output_is_retained(annotate, tmp_path):
    original = query(symbolic=True)
    command = "##VEP-command-line='vep --variant_class'"
    plugin = '##LoFtool="LoFtool score"'
    rows, _, fields = annotate([original], [reference()], extra_headers=[plugin, command])
    lines = (tmp_path / "annotated.vcf").read_text().splitlines()
    custom_headers = [line for line in lines if line.startswith("##INFO=<ID=gnomad")]
    assert custom_headers == [
        '##INFO=<ID=gnomad,Number=.,Type=String,Description="[PATH]/reference.vcf.gz">',
        '##INFO=<ID=gnomad_AF,Number=.,Type=String,Description="AF field from [PATH]/reference.vcf.gz">',
        '##INFO=<ID=gnomad_FREQ_HOMALT,Number=.,Type=String,Description="FREQ_HOMALT field from [PATH]/reference.vcf.gz">',
        '##INFO=<ID=gnomad_SVTYPE,Number=.,Type=String,Description="SVTYPE field from [PATH]/reference.vcf.gz">',
    ]
    assert fields == SCHEMA + GNOMAD_FIELDS
    start = lines.index(plugin)
    assert lines[start:start + 6] == [plugin] + custom_headers + [command]
    expected = original.copy()
    expected[7] += "|r|0.03|0.001|DEL"
    expected[7] = expected[7].replace("|deletion,", "|deletion|r|0.03|0.001|DEL,")
    assert rows == [expected]

    # Re-annotating a file with VEP custom headers must not duplicate them.
    annotate([original], [reference()], extra_headers=custom_headers + [command])
    updated = (tmp_path / "annotated.vcf").read_text().splitlines()
    assert [line for line in updated if line.startswith("##INFO=<ID=gnomad")] == custom_headers


def test_empty_reference_and_empty_query(annotate):
    _, annotations, _ = annotate([query()], [])
    assert annotations[0][0]["gnomad_AF"] == ""
    rows, annotations, fields = annotate([], [reference()])
    assert rows == annotations == []
    assert fields[-4:] == GNOMAD_FIELDS


def test_synonyms_and_no_breakpoint_distance_limit(annotate):
    q = query(length=1000000, symbolic=True)
    q[0] = "NC_000001.11"
    _, annotations, _ = annotate([q], [reference(101000, length=1000000)])
    assert annotations[0][0]["gnomad_AF"] == "0.03"


def test_existing_gnomad_fields_are_replaced_without_duplicate_columns(annotate):
    q = query()
    q[7] = q[7].replace("|deletion", "|deletion|old|0.9|0.8|DEL")
    _, annotations, fields = annotate([q], [reference()], schema=SCHEMA + GNOMAD_FIELDS)
    assert fields == SCHEMA + GNOMAD_FIELDS
    assert annotations[0][0]["gnomad_AF"] == "0.03"


def test_non_symbolic_padding_is_normalized_only_for_matching(annotate):
    q = query(length=3)
    q[3:5] = ["AATTT", "AA"]
    rows, annotations, _ = annotate([q], [reference(1001, length=3)])
    assert rows[0][1:5] == q[1:5]
    assert annotations[0][0]["gnomad_AF"] == "0.03"


def test_unhandled_sequence_sv_representation_fails_instead_of_silently_missing_af(annotate):
    q = query()
    q[4] = "A,AT"
    with pytest.raises(subprocess.CalledProcessError) as error:
        annotate([q], [reference()])
    assert "Cannot represent SV record" in error.value.stderr
