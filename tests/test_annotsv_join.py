"""Regression controls for AnnotSV identity joins and ambiguous legacy tables."""
import csv
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotsv_join import AnnotSVTable


def record(identifier="call1", ref="A" + "C" * 100, alt="A", pos=100, end=None,
           chrom="chr1", kind="DEL"):
    return SimpleNamespace(id=identifier, ref=ref, alts=(alt,), contig=chrom, pos=pos,
                           stop=end if end is not None else pos + len(ref) - 1,
                           info={"SVTYPE": kind})


def row(**changes):
    value = dict(SV_chrom="1", SV_start=100, SV_end=199, SV_type="DEL", ID="call1",
                 REF="A" + "C" * 100, ALT="A", INFO="SVLEN=-100", Annotation_mode="full", ACMG_class="4")
    value.update(changes)
    return value


def table(tmp_path, rows):
    path = tmp_path / "annotsv.tsv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    return AnnotSVTable.read(path)


def test_literal_deletion_without_end_uses_original_identity(tmp_path):
    variant = record()
    before = vars(variant).copy()
    assert table(tmp_path, [row()]).lookup(variant) == {"ACMG_class": "4"}
    assert vars(variant) == before


def test_multi_base_alt_and_inconsistent_svlen_do_not_need_an_offset_guess(tmp_path):
    variant = record(alt="ATT")
    assert table(tmp_path, [row(ALT="ATT", SV_end=197, INFO="SVLEN=-1")]).lookup(variant) == {"ACMG_class": "4"}


@pytest.mark.parametrize("symbolic", [False, True])
def test_explicit_original_end_remains_part_of_identity(tmp_path, symbolic):
    ref, alt = ("A", "<DEL>") if symbolic else ("A" + "C" * 100, "A")
    values = table(tmp_path, [row(REF=ref, ALT=alt, INFO="END=400", SV_end=400)])
    assert values.lookup(record(ref=ref, alt=alt, end=400)) == {"ACMG_class": "4"}
    with pytest.raises(ValueError, match="different allele or coordinates"):
        values.lookup(record(ref=ref, alt=alt, end=401))


@pytest.mark.parametrize("change", [dict(pos=101), dict(chrom="chr2"), dict(ref="G" + "C" * 100),
                                    dict(alt="T"), dict(kind="DUP")])
def test_same_id_with_changed_event_fails_instead_of_borrowing_an_annotation(tmp_path, change):
    with pytest.raises(ValueError, match="different allele or coordinates"):
        table(tmp_path, [row()]).lookup(record(**change))


def test_merged_or_missing_id_requires_an_exact_allele_match(tmp_path):
    values = table(tmp_path, [row()])
    assert values.lookup(record(identifier="merged-id")) == {"ACMG_class": "4"}
    assert values.lookup(record(identifier=None)) == {"ACMG_class": "4"}
    assert values.lookup(record(identifier="unknown", pos=101)) == {}


def test_symbolic_type_and_chromosome_aliases(tmp_path):
    values = table(tmp_path, [row(REF="A", ALT="<DUP>", SV_type="DUP", SV_end=400)])
    assert values.lookup(record(ref="A", alt="<DUP>", kind="TDUP", end=400)) == {"ACMG_class": "4"}


def test_split_rows_do_not_overwrite_the_full_class(tmp_path):
    values = table(tmp_path, [row(), row(Annotation_mode="split", ACMG_class="1")])
    assert values.lookup(record()) == {"ACMG_class": "4"}


@pytest.mark.parametrize("missing", ["NA", ".", ""])
def test_missing_class_remains_missing(tmp_path, missing):
    assert table(tmp_path, [row(ACMG_class=missing)]).lookup(record()) == {}


def test_conflicting_full_rows_for_an_id_fail(tmp_path):
    with pytest.raises(ValueError, match="Conflicting full"):
        table(tmp_path, [row(), row(ACMG_class="5")])


def test_allele_fallback_rejects_conflicting_annotations(tmp_path):
    values = table(tmp_path, [row(), row(ID="other", ACMG_class="5")])
    assert values.lookup(record()) == {"ACMG_class": "4"}
    with pytest.raises(ValueError, match="Ambiguous"):
        values.lookup(record(identifier="merged"))


def test_coordinate_only_tables_keep_exact_matching_and_reject_conflicts(tmp_path):
    original = {key: value for key, value in row(SV_end=200).items()
                if key not in ("REF", "ALT", "ID", "INFO")}
    values = table(tmp_path, [original])
    assert values.lookup(record()) == {"ACMG_class": "4"}
    assert values.lookup(record(end=201)) == {}
    values = table(tmp_path, [original, {**original, "ACMG_class": "5"}])
    with pytest.raises(ValueError, match="Ambiguous"):
        values.lookup(record())


def test_large_info_fields_and_csv_limit_are_preserved(tmp_path):
    previous = csv.field_size_limit()
    assert table(tmp_path, [row(INFO="CSQ=" + "A" * 300000)]).lookup(record()) == {"ACMG_class": "4"}
    assert csv.field_size_limit() == previous


@pytest.mark.parametrize("header", ["1\n", "AnnotSV_ID\tACMG_class\n"])
def test_explicit_empty_outputs(tmp_path, header):
    path = tmp_path / "empty.tsv"; path.write_text(header)
    assert AnnotSVTable.read(path).lookup(record()) == {}


def test_real_pysam_records_keep_genotypes_and_phasing(tmp_path):
    pysam = pytest.importorskip("pysam")
    path = tmp_path / "input.vcf"
    path.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n'
                    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n'
                    'chr1\t100\tcall1\t' + 'A' + 'C' * 100 + '\tA\t50\tPASS\tSVTYPE=DEL\tGT\t0|1\n')
    with pysam.VariantFile(path) as vcf:
        variant = next(vcf); before = str(variant)
        assert variant.stop == 200
        assert table(tmp_path, [row()]).lookup(variant) == {"ACMG_class": "4"}
        assert str(variant) == before
