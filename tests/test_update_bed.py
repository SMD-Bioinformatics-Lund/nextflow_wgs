from pathlib import Path
import gzip

from bin.reference_tools.update_bed import (
    write_ensembl_bed,
    append_to_bed,
    compare_clinvar,
    read_clinvar,
    ClinVarVariant,
    main,
)

mock_gtf_content = """\
# Example GTF file
1\tENSEMBL\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gene1"; transcript_biotype "protein_coding";
1\tENSEMBL\texon\t1100\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
2\tENSEMBL\ttranscript\t3000\t4000\t.\t+\t.\tgene_id "gene2"; transcript_biotype "non_coding";
2\tENSEMBL\texon\t3100\t3200\t.\t+\t.\tgene_id "gene2"; transcript_id "transcript2";
3\tENSEMBL\ttranscript\t4000\t5000\t.\t+\t.\tgene_id "gene3"; transcript_biotype "protein_coding";
3\tENSEMBL\texon\t4400\t4500\t.\t+\t.\tgene_id "gene3"; transcript_id "transcript3";
"""

mock_chr_gtf_content = """\
# Example GTF file
chr1\tENSEMBL\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "gene1"; transcript_biotype "protein_coding";
chr1\tENSEMBL\texon\t1100\t1200\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";
chr2\tENSEMBL\ttranscript\t3000\t4000\t.\t+\t.\tgene_id "gene2"; transcript_biotype "non_coding";
chr2\tENSEMBL\texon\t3100\t3200\t.\t+\t.\tgene_id "gene2"; transcript_id "transcript2";
chr3\tENSEMBL\ttranscript\t4000\t5000\t.\t+\t.\tgene_id "gene3"; transcript_biotype "protein_coding";
chr3\tENSEMBL\texon\t4400\t4500\t.\t+\t.\tgene_id "gene3"; transcript_id "transcript3";
"""


def get_file_string(rows: list[list[str]]) -> str:
    """
    Convert
    [['chr1', '1', '2'], ['chr2', '2', '3']]
    to
    'chr1\t1\t2\nchr2\t2\t3'
    """
    row_strs: list[str] = []
    for row in rows:
        row_str = "\t".join(row)
        row_strs.append(row_str)

    return "\n".join(row_strs) + "\n"


def write_gz(path: Path, content: str):
    path.unlink(missing_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write(content)


def test_write_ensembl_bed(tmp_path: Path):
    """
    Checks that exons are extracted and padded correctly
    Here, a pre-downloaded GTF is provided

    Checks both GTF with chromosomes chr1, chr2, chr2 syntax and 1, 2, 3
    """

    release = "108"
    skip_download = True
    padding = 20
    out_bed = tmp_path / "output.bed"

    # Non-chr based reference
    write_gz(tmp_path / "tmp.gtf.gz", mock_gtf_content)

    write_ensembl_bed(tmp_path, out_bed, release, skip_download, padding)
    output_content = out_bed.read_text()
    row_fields = [["1", "1080", "1220"], ["3", "4380", "4520"]]
    expected_output = get_file_string(row_fields)
    assert output_content == expected_output, f"Output was: {output_content}"

    # Chr based reference
    out_bed2 = tmp_path / "output2.bed"
    write_gz(tmp_path / "tmp.gtf.gz", mock_chr_gtf_content)
    write_ensembl_bed(tmp_path, out_bed2, release, skip_download, padding)
    output_content = out_bed2.read_text()
    row_fields = [["chr1", "1080", "1220"], ["chr3", "4380", "4520"]]
    expected_output = get_file_string(row_fields)
    assert output_content == expected_output, f"Output was: {output_content}"


def test_append_to_bed(tmp_path: Path):
    """Checks that additional bed content is appended correctly"""

    base_bed_content = """\
chr1\t1\t5\told_annot
chr1\t10\t15\told_annot2
"""

    mock_bed_content = """\
chr1\t1000\t2000
chr1\t3000\t4000\tgene1
chr2\t5000\t6000
"""

    expected_output = """\
chr1\t1\t5\told_annot
chr1\t10\t15\told_annot2
chr1\t1000\t2000\tdefault_annot
chr1\t3000\t4000\tgene1
chr2\t5000\t6000\tdefault_annot
"""

    out_bed_path = tmp_path / "out_bed.bed"
    bed2add_path = tmp_path / "bed2add.bed"

    out_bed_path.write_text(base_bed_content)
    bed2add_path.write_text(mock_bed_content)

    fourth_col_default = "default_annot"
    append_to_bed(out_bed_path, bed2add_path, fourth_col_default)
    result = out_bed_path.read_text()

    assert result == expected_output, f"Output was: {result}"


def test_read_clinvar(tmp_path: Path):
    """
    Given a minimal VCF, check that the 'read_clinvar' function successfully separates
    the pathogenic and non-pathogenic variants
    """

    mock_vcf_content = """\
##fileformat=VCFv4.2
##source=ClinVar
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t12345\t.\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNSIGCONF=Conflicting_interpretations_of_pathogenicity;CLNSIGINCL=athogenic
chr1\t54321\t.\tC\tT\t.\t.\tCLNSIG=Benign
chr2\t12345\t.\tG\tA\t.\t.\tCLNSIG=Likely_pathogenic
"""

    clinvar_vcf_path = tmp_path / "tmp.vcf"
    write_gz(clinvar_vcf_path, mock_vcf_content)

    clinvar_variants, benign_clinvar_variants = read_clinvar(clinvar_vcf_path)

    expected_clinvar_reasons = {
        "chr1:12345_A_G": "Pathogenic",
        "chr2:12345_G_A": "Likely_pathogenic",
    }
    expected_benign_clinvar_reasons = {"chr1:54321_C_T": "Undefined"}

    for key, variant in clinvar_variants.items():
        assert key in expected_clinvar_reasons
        assert expected_clinvar_reasons[key] == variant.reason

    for key, variant in benign_clinvar_variants.items():
        assert key in expected_benign_clinvar_reasons
        assert expected_benign_clinvar_reasons[key] == variant.reason


def test_compare_clinvar(tmp_path: Path):
    """Given mock entries with different sets of clinvar variants, check what is added, removed and retained"""

    new_clinvar = {
        "chr1:1000_A_G": ClinVarVariant(
            "chr1", 1000, "A", "G", {"CLDN": "CLDN", "CLNACC": "CLNACC"}
        ),
        "chr2:3000_C_T": ClinVarVariant(
            "chr2", 3000, "C", "T", {"CLDN": "CLDN", "CLNACC": "CLNACC"}
        ),
        "chr3:5000_G_A": ClinVarVariant(
            "chr3", 5000, "G", "A", {"CLDN": "CLDN", "CLNACC": "CLNACC"}
        ),
    }

    old_clinvar = {
        "chr1:1000_A_G": ClinVarVariant(
            "chr1", 1000, "A", "G", {"CLDN": "CLDN", "CLNACC": "CLNACC"}
        ),
        "chr2:2000_A_G": ClinVarVariant(
            "chr2", 2000, "A", "G", {"CLDN": "CLDN", "CLNACC": "CLNACC"}
        ),
        "chr2:5000_C_T": ClinVarVariant(
            "chr2", 5000, "C", "T", {"CLDN": "CLDN", "CLNACC": "CLNACC"}
        ),
    }

    expected_new_to_add = [
        ["chr2", "2995", "3005", "CLINVAR-Undefined", "chr2:3000_C_T"],
        ["chr3", "4995", "5005", "CLINVAR-Undefined", "chr3:5000_G_A"],
    ]
    expected_old_to_remove = [
        ["chr2", "1995", "2005", "CLINVAR-Undefined", "chr2:2000_A_G"]
    ]

    final_bed_path = tmp_path / "final.bed"

    bed_test_content = """\
chr1\t1\t5\told_annot
chr1\t10\t15\told_annot2
chr1\t1000\t2000\tdefault_annot
chr1\t3000\t4000\tgene1
chr2\t5000\t6000\tdefault_annot
"""

    final_bed_path.write_text(bed_test_content)

    out_dir = tmp_path
    padding = 5

    current_new_to_add, current_old_to_remove = compare_clinvar(
        new_clinvar, old_clinvar, final_bed_path, padding, out_dir, keep_tmp=False
    )

    assert current_new_to_add == expected_new_to_add
    assert current_old_to_remove == expected_old_to_remove


def test_main(tmp_path: Path):
    """Check that the full script runs and check that outputs are present with correct information"""

    old_clinvar_content = """\
##fileformat=VCFv4.2
##source=ClinVar
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1000\t.\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNSIGCONF=Conflicting_interpretations_of_pathogenicity
2\t2000\t.\tA\tG\t.\t.\tCLNSIG=Benign
3\t3000\t.\tC\tT\t.\t.\tCLNSIG=Likely_benign
"""

    new_clinvar_content = """\
##fileformat=VCFv4.2
##source=ClinVar
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1000\t.\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNSIGCONF=Conflicting_interpretations_of_pathogenicity
2\t3000\t.\tC\tT\t.\t.\tCLNSIG=Likely_pathogenic
3\t5000\t.\tG\tA\t.\t.\tCLNSIG=Pathogenic
3\t5000\t.\tG\tA\t.\t.\tCLNSIG=Pathogenic
"""

    expected_results = """1\t995\t1005\tCLINVAR-Pathogenic
1\t1080\t1220\tEXONS-release
2\t2995\t3005\tCLINVAR-Likely_pathogenic
3\t4380\t4520\tEXONS-release
3\t4995\t5005\tCLINVAR-Pathogenic
"""

    old_clinvar = tmp_path / "clinvar_old.vcf"
    write_gz(old_clinvar, old_clinvar_content)
    new_clinvar = tmp_path / "clinvar_new.vcf"
    write_gz(new_clinvar, new_clinvar_content)
    clinvardate = "clinvardate"
    ensembl_release = "release"
    write_gz(tmp_path / "tmp.gtf.gz", mock_gtf_content)

    skip_download = True
    incl_bed: list[str] = []

    main(
        new_clinvar,
        old_clinvar,
        str(clinvardate),
        tmp_path,
        ensembl_release,
        skip_download,
        incl_bed,
        keep_tmp=False,
    )

    exon_pad = 20
    variant_pad = 5
    final_bed = (
        tmp_path
        / f"exons_{ensembl_release}padded{exon_pad}bp_clinvar-{clinvardate}padded{variant_pad}bp.bed"
    )

    assert final_bed is not None
    assert final_bed.read_text() == expected_results
