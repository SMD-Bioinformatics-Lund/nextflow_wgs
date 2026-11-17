import gzip
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "bin/gnomad_sv.py"
MODES = ("overlap", "ins_site", "ins_length")
FIELDS = ["gnomad", "gnomad_AF", "gnomad_FREQ_HOMALT", "gnomad_SVTYPE"]
SCHEMA = ["Allele", "Consequence", "SYMBOL", "VARIANT_CLASS"]


def variant(pos=1000, length=500, kind="DEL", symbolic=False, name="q", csq=True):
    ref, alt = ("A" + "T" * length, "A") if kind == "DEL" else ("A", "A" + "T" * length)
    if symbolic:
        ref, alt = "N", "<{}>".format(kind)
    info = "END={};SVLEN={};SVTYPE={};FLAG".format(
        pos + length if kind == "DEL" else pos, -length if kind == "DEL" else length, kind)
    if csq:
        cls = "deletion" if kind == "DEL" else "insertion"
        info += ";CSQ=-|intron_variant|GENE1|{},-|splice_region_variant|GENE2|{}".format(cls, cls)
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
        vcf, source, database, prefix, output = (
            tmp_path / name for name in ("input.vcf", "reference.vcf.gz", "reference", "query", "output.vcf"))
        write_vcf(vcf, queries, csq)
        write_vcf(source, references)

        def helper(mode, input_path, output_path, *extra):
            subprocess.run([sys.executable, str(SCRIPT), mode, "--vcf", str(input_path),
                            "--output", str(output_path), *extra], check=True, capture_output=True, text=True)

        helper("reference", source, database)
        helper("query", vcf, prefix)
        matches = []
        for mode in MODES:
            intervals = "{}.{}.intervals.vcf".format(prefix, mode)
            match = tmp_path / (mode + ".matches.vcf")
            subprocess.run(["bcftools", "annotate", "--no-version", "-a", str(database / (mode + ".bed.gz")),
                            "-c", "CHROM,FROM,TO,GNOMAD_SV_MATCHES", "-l", "GNOMAD_SV_MATCHES:append",
                            "--min-overlap", "0:0" if mode == "ins_site" else "0.699999:0.699999",
                            "-Ov", "-o", str(match), intervals], check=True, capture_output=True, text=True)
            matches.append(str(match))
        helper("restore", vcf, output, "--matches", *matches)
        lines = output.read_text().splitlines()
        rows = [line.split("\t") for line in lines if not line.startswith("#")]
        annotations = []
        for row in rows:
            info = dict(entry.split("=", 1) for entry in row[7].split(";") if "=" in entry)
            annotations.append([dict({field: "" for field in FIELDS}, **info)])
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
    assert rows[0][7].split("CSQ=", 1)[1].split(";", 1)[0] == original[7].split("CSQ=", 1)[1]
    assert len([line for line in lines if line.startswith("##INFO=<ID=CSQ,")]) == 1
    assert '##INFO=<ID=gnomad_AF,Number=.,Type=Float,' in "\n".join(lines)
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
    assert annotations[0][0]["gnomad"] == "exact,shifted"
    assert annotations[0][0]["gnomad_AF"] == "0.0001,0.03"


@pytest.mark.parametrize("symbolic", [False, True])
def test_insertions_match_the_ordinary_insertion_class(annotate, symbolic):
    original = variant(kind="INS", symbolic=symbolic)
    rows, annotations, _ = annotate([original], [reference(kind="INS"), reference(kind="DEL", af="0.9")])
    assert rows[0][:7] == original[:7]
    assert rows[0][8:] == original[8:]
    assert [annotations[0][0][field] for field in FIELDS] == ["r", "0.03", "0.001", "INS"]


@pytest.mark.parametrize("pos, length, expected", [(1000, 1000, "r"), (1001, 100, ""), (999, 100, "")])
def test_sequence_insertions_keep_master_site_matching_without_size_filter(annotate, pos, length, expected):
    _, annotations, _ = annotate([variant(kind="INS", length=100)], [reference(pos, length, kind="INS")])
    assert annotations[0][0]["gnomad"] == expected


@pytest.mark.parametrize("pos, length, expected", [(1030, 100, "r"), (1031, 100, ""), (1000, 200, "")])
def test_structural_insertions_keep_reciprocal_overlap(annotate, pos, length, expected):
    _, annotations, _ = annotate([variant(kind="INS", length=100, symbolic=True)], [reference(pos, length, kind="INS")])
    assert annotations[0][0]["gnomad"] == expected


def test_sequence_insertion_padding_determines_the_site(annotate):
    original = variant(kind="INS", length=3)
    original[3:5] = ["AA", "AATAA"]
    rows, annotations, _ = annotate([original], [reference(1001, kind="INS")])
    assert rows[0][1:5] == original[1:5]
    assert annotations[0][0]["gnomad"] == "r"


def test_ins_site_uses_raw_reference_span_not_inserted_length(annotate):
    ref = reference(999, kind="INS")
    ref[3] = "NN"
    _, annotations, _ = annotate([variant(kind="INS")], [ref])
    assert annotations[0][0]["gnomad"] == "r"


@pytest.mark.parametrize("index_end, expected", [(1000, ""), (1001, "r")])
def test_sequence_ins_preserves_vep_reference_index_boundary(annotate, index_end, expected):
    ref = reference(kind="INS")
    ref[7] = ref[7].replace("END=1500", "END={}".format(index_end))
    _, annotations, _ = annotate([variant(kind="INS")], [ref])
    assert annotations[0][0]["gnomad"] == expected


@pytest.mark.parametrize("index_end, expected", [(950, ""), (1400, "r")])
def test_structural_ins_preserves_vep_reference_index_boundary(annotate, index_end, expected):
    ref = reference(900, kind="INS")
    ref[7] = ref[7].replace("END=1400", "END={}".format(index_end))
    _, annotations, _ = annotate([variant(kind="INS", symbolic=True)], [ref])
    assert annotations[0][0]["gnomad"] == expected


def test_unknown_size_structural_ins_is_not_assigned_an_inserted_length(annotate):
    original = variant(kind="INS", symbolic=True)
    original[7] = original[7].replace(";SVLEN=500", "")
    _, annotations, _ = annotate([original], [reference(length=50, kind="INS")])
    assert annotations[0][0]["gnomad"] == ""


def test_manta_mixed_insertion_with_ns_takes_the_structural_path(annotate):
    original = variant(kind="INS", length=201, symbolic=True)
    original[3:5] = ["A" + "N" * 20, "A" + "T" * 201]
    original[7] = original[7].replace("END=1000", "END=1020") + ";CIGAR=1M201I20D"
    rows, annotations, _ = annotate([original], [reference(1060, length=201, kind="INS")])
    assert rows[0][3:5] == original[3:5]
    assert annotations[0][0]["gnomad"] == "r"


@pytest.mark.parametrize("kind, cls", [
    ("DUP", "duplication"), ("INV", "inversion"), ("CNV", "copy_number_variation"),
    ("DUP:TANDEM", "tandem_duplication"), ("INS:ME:ALU", "Alu_insertion"),
    ("INS:ME:LINE1", "LINE1_insertion"), ("DEL:ME:LINE1", "LINE1_deletion"),
    ("BND", "chromosome_breakpoint"),
])
def test_other_sv_classes_get_variant_info_without_vep_frequency_fields(annotate, kind, cls):
    original = variant(kind="INS", symbolic=True)
    original[4] = "<{}>".format(kind)
    original[7] = original[7].replace("|insertion", "|" + cls).replace("SVTYPE=INS", "SVTYPE=" + kind.split(":", 1)[0])
    rows, annotations, _ = annotate([original], [reference(kind=kind)])
    assert rows[0][:7] == original[:7]
    assert rows[0][8:] == original[8:]
    assert annotations[0][0]["gnomad"] == "r"
    assert annotations[0][0]["gnomad_AF"] == "0.03"
    assert annotations[0][0]["CSQ"] == original[7].split("CSQ=", 1)[1]


def test_tandem_duplication_does_not_match_ordinary_duplication(annotate):
    original = variant(symbolic=True)
    original[4] = "<DUP:TANDEM>"
    original[7] = original[7].replace("SVTYPE=DEL", "SVTYPE=DUP").replace("|deletion", "|tandem_duplication")
    _, annotations, _ = annotate([original], [reference(kind="DUP")])
    assert annotations[0][0]["gnomad"] == ""


def test_nonstandard_tdup_alt_uses_svtype_dup(annotate):
    original = variant(symbolic=True)
    original[4] = "<TDUP>"
    original[7] = original[7].replace("SVTYPE=DEL", "SVTYPE=DUP").replace("|deletion", "|duplication")
    _, annotations, _ = annotate([original], [reference(kind="DUP")])
    assert annotations[0][0]["gnomad"] == "r"


def test_breakend_without_svlen_ignores_mate_end(annotate):
    original = variant(kind="INS", symbolic=True)
    original[4] = "<BND>"
    original[7] = original[7].replace("SVTYPE=INS", "SVTYPE=BND").replace("|insertion", "|chromosome_breakpoint")
    original[7] = original[7].replace(";SVLEN=500", "").replace("END=1000", "END=2000000")
    _, annotations, _ = annotate([original], [reference(kind="BND", length=1)])
    assert annotations[0][0]["gnomad"] == "r"


def test_reannotation_replaces_info_and_retains_one_typed_header(annotate):
    original = variant()
    original[7] += ";gnomad=stale;gnomad_AF=0.9"
    rows, annotations, lines = annotate([original], [reference()])
    assert rows[0][7].count("gnomad_AF=") == 1
    assert annotations[0][0]["gnomad_AF"] == "0.03"
    assert len([line for line in lines if line.startswith("##INFO=<ID=gnomad_AF,")]) == 1


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
    assert annotations[0][0]["gnomad"].split(",") == ["r{}".format(i) for i in range(50)] + ["..."]
    assert max(map(float, annotations[0][0]["gnomad_AF"].split(","))) == 0.0001


def test_missing_values_keep_reference_fields_aligned_in_info(annotate):
    refs = [reference(name="first;second", af="0", hom="."), reference(name="r2", af="0.97,0.03")]
    refs[1][7] = refs[1][7].replace(";FREQ_HOMALT=0.001", "")
    _, annotations, _ = annotate([variant()], refs)
    assert [annotations[0][0][field] for field in FIELDS] == ["first,r2", "0,0.03", ".,.", "DEL,DEL"]


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
