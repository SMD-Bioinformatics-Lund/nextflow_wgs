#!/usr/bin/env python
from collections import defaultdict
import json
import subprocess
import logging
import gzip
import argparse
"""
required
    mane.gtf
    bed-file containing probe placement of interest (genefocused, not backbone LoH SNPs), or genefilter for WGS samples
    bam-file
"""

def main():
    logging.debug(f"Hi!")
    args = parse_arguments()   
    sample_id = args.sample_id
    output_prefix = sample_id
    bam_file = args.bam
    mane_gtf = args.gtf
    exclude_y_regions = args.sex == "F"
    # Panel
    if args.design_bed:
        design_bed = args.design_bed
        cds_file, exon_file = intersect_gtf_with_design(
            mane_gtf,
            design_bed,
            output_prefix,
            exclude_y_regions,
        )
        gene_gtf,exon_to_gene,cds_to_gene,probe_mapper,_wgs_gtf = read_mane_gtf(mane_gtf)
        if exclude_y_regions:
            prune_y_chromosome_regions(gene_gtf, exon_to_gene, cds_to_gene, probe_mapper)
        if args.caveat_genes:
            caveat_genes = read_caveat_gene_list(args.caveat_genes)
            mark_partial_cds_genes(gene_gtf, caveat_genes)
    # WGS
    else:
        gene_filter = args.gene_filter
        gene_list, partial_cds_genes = read_coverage_gene_list(gene_filter)
        filtered_gtf = f"{output_prefix}.filtered.gtf"
        gene_gtf,exon_to_gene,cds_to_gene,probe_mapper,wgs_gtf = read_mane_gtf(
            mane_gtf,
            gene_list,
            filtered_gtf,
        )
        if exclude_y_regions:
            prune_y_chromosome_regions(gene_gtf, exon_to_gene, cds_to_gene, probe_mapper)
        mark_partial_cds_genes(gene_gtf, partial_cds_genes)
        cds_file,exon_file = generate_bed_regions(wgs_gtf, output_prefix, exclude_y_regions)

    # calculate coverage for sub-regions using mosdepth
    cds_cov   = query_cov_for_region(bam_file,sample_id,cds_file,"cds")
    exon_cov  = query_cov_for_region(bam_file,sample_id,exon_file,"exons")
    gene_gtf = open_case_cov(gene_gtf,cds_to_gene,cds_cov,"CDS")
    gene_gtf = open_case_cov(gene_gtf,exon_to_gene,exon_cov,"exons")
    if args.design_bed:
        probe_bed = filter_y_chromosome_bed(design_bed, output_prefix, "probes") if exclude_y_regions else design_bed
        probe_cov = query_cov_for_region(bam_file,sample_id,probe_bed,"probes")
        gene_gtf = assign_probes(gene_gtf,probe_mapper,probe_cov)

    # only keep genes with some kind of coverage, probe or gene filtered for wgs    
    filtered_gene_gtf = {key: value for key, value in gene_gtf.items() if value.get("covered_by_panel") is True}
    
    mongolike = {}
    mongolike['genes'] = filtered_gene_gtf
    if args.summary_genes:
        genes_to_include = read_gene_list(args.summary_genes)
        summary = summarize_coverage(mongolike, genes_to_include, args.threshold)
        mongolike["summary"] = summary
        if args.summary_output:
            write_json(args.summary_output, summary)
    json_output = json.dumps(mongolike, indent=4)
    output_json = f"{sample_id}.cov.json"
    if args.output:
        output_json = args.output
    logging.debug(f"printing output to {output_json}")
    with open(output_json, 'w') as outfile:
        outfile.write(json_output)
    logging.debug(f"Bye!")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Create coverage stats for coyote 3.0")
    parser.add_argument(
        '-b', '--bam',
        type=str,
        required=True,
        help="file path to bam-file"
    )
    parser.add_argument(
        '-g', '--gtf',
        type=str,
        required=True,
        help="file path to gtf-file"
    )
    parser.add_argument(
        '-d', '--design_bed',
        type=str,
        help="design probes, only genes overlapping these will be kept, mutually exclusive to gene_filter"
    )
    parser.add_argument(
        '--caveat_genes',
        type=str,
        help=(
            "gene list to skip from summary when using design_bed. "
            "One-column entries are treated as caveats; optional second column accepts "
            "true/false, yes/no, partial/full"
        )
    )
    parser.add_argument(
        '-s', '--sample_id',
        type=str,
        required=True,
        help="name of sample"
    )
    parser.add_argument(
        '-f', '--gene_filter',
        type=str,
        help=(
            "list of genes to keep, mutually exclusive to design_bed. "
            "Optional second column marks partial CDS coverage: true/false, yes/no, partial/full"
        )
    )
    parser.add_argument(
        '--summary_genes',
        type=str,
        help="list of genes to include in the coverage summary"
    )
    parser.add_argument(
        '-t', '--threshold',
        default=500.0,
        type=parse_thresholds,
        help="Mean CDS coverage threshold(s) for summary. Use comma-separated values for multiple thresholds. Default: 500"
    )
    parser.add_argument(
        '--summary_output',
        type=str,
        help="write summary JSON to this file in addition to embedding it in the main JSON output"
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        help="name of json output file. If not provided sample_id.cov.json"
    )
    parser.add_argument(
        '--sex',
        type=parse_sex,
        choices=["M", "F"],
        help="sample sex (M/F). Female samples ignore Y chromosome regions"
    )
    args = parser.parse_args()
    if args.design_bed and args.gene_filter:
        exit("Cannot use both gene filter and design bed")
    elif args.design_bed is None and args.gene_filter is None:
        exit("Please provide either design_bed or gene_filter")
    elif args.caveat_genes and not args.design_bed:
        exit("caveat_genes can only be used together with design_bed")
    return args

def parse_sex(value):
    sex = value.strip().upper()
    if sex not in {"M", "F"}:
        raise argparse.ArgumentTypeError("--sex must be either 'M' or 'F'")
    return sex

def write_json(output_file, data):
    with open(output_file, "w", encoding="utf-8") as outfile:
        outfile.write(json.dumps(data, indent=4) + "\n")

def parse_partial_cds_flag(value):
    normalized = value.strip().lower()
    if normalized in {"1", "true", "t", "yes", "y", "partial", "partially_covered", "partial_cds"}:
        return True
    if normalized in {"0", "false", "f", "no", "n", "full", "complete", "fully_covered", "full_cds"}:
        return False
    raise ValueError(
        f"Unknown partial CDS flag '{value}'. Use true/false, yes/no, partial/full, or 1/0."
    )

def parse_thresholds(value):
    thresholds = []
    for raw_threshold in str(value).split(","):
        raw_threshold = raw_threshold.strip()
        if not raw_threshold:
            continue
        try:
            thresholds.append(float(raw_threshold))
        except ValueError as error:
            raise argparse.ArgumentTypeError(
                f"Invalid threshold '{raw_threshold}'. Use a number or comma-separated numbers."
            ) from error

    if not thresholds:
        raise argparse.ArgumentTypeError("At least one threshold is required.")

    return thresholds[0] if len(thresholds) == 1 else thresholds

def read_coverage_gene_list(gene_file):
    """
    Reads the coverage gene list used to select genes for coverage calculation.
    Column 1 is gene symbol. Optional column 2 marks genes where only part of the
    CDS is covered by the design.
    """
    gene_list = []
    partial_cds_genes = set()
    seen = set()
    with open(gene_file, 'r', encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            gene = fields[0]
            if len(fields) > 2:
                raise ValueError(
                    f"{gene_file}:{line_number} has too many columns. Expected gene and optional partial CDS flag."
                )
            if gene not in seen:
                gene_list.append(gene)
                seen.add(gene)
            if len(fields) == 2 and parse_partial_cds_flag(fields[1]):
                partial_cds_genes.add(gene)
    return gene_list, partial_cds_genes

def read_caveat_gene_list(gene_file):
    """
    Reads genes that should be skipped from summary assessment in design-bed mode.
    A single-column file marks every listed gene as a caveat. If a second column is
    present, only truthy partial/caveat values are marked.
    """
    caveat_genes = set()
    with open(gene_file, 'r', encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) > 2:
                raise ValueError(
                    f"{gene_file}:{line_number} has too many columns. Expected gene and optional caveat flag."
                )
            if len(fields) == 1 or parse_partial_cds_flag(fields[1]):
                caveat_genes.add(fields[0])
    return caveat_genes

def read_gene_list(gene_file):
    genes = []
    seen = set()
    with open(gene_file, 'r', encoding="utf-8") as file:
        for line in file:
            gene = line.strip()
            if not gene or gene.startswith("#"):
                continue
            if gene not in seen:
                genes.append(gene)
                seen.add(gene)
    return genes

def mark_partial_cds_genes(gene_gtf, partial_cds_genes):
    for gene in partial_cds_genes:
        gene_gtf[gene]["partial_cds_coverage"] = True
    for gene in gene_gtf:
        gene_gtf[gene].setdefault("partial_cds_coverage", False)
    return gene_gtf

def region_length(region):
    start = int(region["start"])
    end = int(region["end"])
    return max(end - start, 0)

def region_coverage(region):
    if "cov" not in region:
        return 0.0
    return float(region["cov"])

def weighted_mean_cds_coverage(cds_regions):
    total_bases = 0
    coverage_sum = 0.0

    for region in cds_regions.values():
        length = region_length(region)
        total_bases += length
        coverage_sum += region_coverage(region) * length

    if total_bases == 0:
        return None

    return coverage_sum / total_bases

def normalize_thresholds(threshold):
    if isinstance(threshold, (list, tuple)):
        thresholds = [float(value) for value in threshold]
    else:
        thresholds = [float(threshold)]

    if not thresholds:
        raise ValueError("At least one threshold is required.")

    return thresholds

def summarize_threshold_results(gene_results, threshold):
    genes_below_threshold = []
    covered_genes = 0
    assessed_genes = 0

    for gene_result in gene_results:
        assessment = gene_result["coverage_assessment"]
        if assessment == "partial_cds_not_assessed":
            continue

        assessed_genes += 1
        mean_coverage = gene_result["mean_cds_coverage"]
        is_covered = mean_coverage is not None and mean_coverage >= threshold

        if is_covered:
            covered_genes += 1
        else:
            genes_below_threshold.append(gene_result["gene"])

    percent_covered = (covered_genes / assessed_genes * 100) if assessed_genes else 0.0

    return {
        "threshold": threshold,
        "genes_assessed_for_full_cds_mean_coverage": assessed_genes,
        "genes_covered_by_mean_cds_coverage_threshold": covered_genes,
        "percent_genes_covered_by_at_least_threshold_mean_cds_coverage": round(
            percent_covered, 2
        ),
        "genes_not_covered_by_at_least_threshold_mean_cds_coverage": genes_below_threshold,
    }

def summarize_coverage(coverage_data, genes_to_include, threshold):
    thresholds = normalize_thresholds(threshold)
    primary_threshold = thresholds[0]
    genes = coverage_data.get("genes", {})

    gene_results = []
    threshold_assessment_records = []
    missing_genes = []
    partial_cds_genes = []
    assessed_genes = 0

    for gene in genes_to_include:
        gene_record = genes.get(gene)
        if gene_record is None:
            missing_genes.append(gene)
            assessed_genes += 1
            threshold_assessment_records.append(
                {
                    "gene": gene,
                    "mean_cds_coverage": None,
                    "coverage_assessment": "missing_gene",
                }
            )
            gene_results.append(
                {
                    "gene": gene,
                    "found": False,
                    "partial_cds_coverage": False,
                    "mean_cds_coverage": None,
                    "covered_by_mean_cds_coverage_threshold": False,
                    "covered_by_mean_cds_coverage_thresholds": {
                        str(threshold): False for threshold in thresholds
                    },
                    "coverage_assessment": "missing_gene",
                }
            )
            continue

        cds_regions = gene_record.get("CDS", {})
        mean_coverage = weighted_mean_cds_coverage(cds_regions)
        is_partial_cds = bool(gene_record.get("partial_cds_coverage", False))

        if is_partial_cds:
            partial_cds_genes.append(gene)
            threshold_assessment_records.append(
                {
                    "gene": gene,
                    "mean_cds_coverage": mean_coverage,
                    "coverage_assessment": "partial_cds_not_assessed",
                }
            )
            gene_results.append(
                {
                    "gene": gene,
                    "found": True,
                    "partial_cds_coverage": True,
                    "mean_cds_coverage": None,
                    "available_cds_regions_mean_coverage": (
                        round(mean_coverage, 2) if mean_coverage is not None else None
                    ),
                    "covered_by_mean_cds_coverage_threshold": None,
                    "covered_by_mean_cds_coverage_thresholds": {
                        str(threshold): None for threshold in thresholds
                    },
                    "coverage_assessment": "partial_cds_not_assessed",
                    "cds_regions": len(cds_regions),
                    "cds_regions_without_coverage": sum(
                        1 for region in cds_regions.values() if "cov" not in region
                    ),
                    "cds_regions_below_threshold": [
                        region_name
                        for region_name, region in cds_regions.items()
                        if region_coverage(region) < primary_threshold
                    ],
                    "cds_regions_below_thresholds": {
                        str(threshold): [
                            region_name
                            for region_name, region in cds_regions.items()
                            if region_coverage(region) < threshold
                        ]
                        for threshold in thresholds
                    },
                }
            )
            continue

        assessed_genes += 1
        is_covered = mean_coverage is not None and mean_coverage >= primary_threshold
        threshold_assessment_records.append(
            {
                "gene": gene,
                "mean_cds_coverage": mean_coverage,
                "coverage_assessment": "assessed_full_cds",
            }
        )

        gene_results.append(
            {
                "gene": gene,
                "found": True,
                "partial_cds_coverage": False,
                "mean_cds_coverage": round(mean_coverage, 2) if mean_coverage is not None else None,
                "covered_by_mean_cds_coverage_threshold": is_covered,
                "covered_by_mean_cds_coverage_thresholds": {
                    str(threshold): mean_coverage is not None and mean_coverage >= threshold
                    for threshold in thresholds
                },
                "coverage_assessment": "assessed_full_cds",
                "cds_regions": len(cds_regions),
                "cds_regions_without_coverage": sum(
                    1 for region in cds_regions.values() if "cov" not in region
                ),
                "cds_regions_below_threshold": [
                    region_name
                    for region_name, region in cds_regions.items()
                    if region_coverage(region) < primary_threshold
                ],
                "cds_regions_below_thresholds": {
                    str(threshold): [
                        region_name
                        for region_name, region in cds_regions.items()
                        if region_coverage(region) < threshold
                    ]
                    for threshold in thresholds
                },
            }
        )

    total_genes = len(genes_to_include)
    threshold_summaries = [
        summarize_threshold_results(threshold_assessment_records, threshold)
        for threshold in thresholds
    ]
    primary_threshold_summary = threshold_summaries[0]

    summary = {
        "threshold": primary_threshold,
        "genes_requested": total_genes,
        "genes_found": total_genes - len(missing_genes),
        "genes_assessed_for_full_cds_mean_coverage": assessed_genes,
        "genes_skipped_due_to_partial_cds_coverage": partial_cds_genes,
        "genes_covered_by_mean_cds_coverage_threshold": primary_threshold_summary[
            "genes_covered_by_mean_cds_coverage_threshold"
        ],
        "percent_genes_covered_by_at_least_threshold_mean_cds_coverage": primary_threshold_summary[
            "percent_genes_covered_by_at_least_threshold_mean_cds_coverage"
        ],
        "genes_not_covered_by_at_least_threshold_mean_cds_coverage": primary_threshold_summary[
            "genes_not_covered_by_at_least_threshold_mean_cds_coverage"
        ],
        "missing_genes": missing_genes,
        "gene_results": gene_results,
    }

    if len(thresholds) > 1:
        summary["thresholds"] = thresholds
        summary["threshold_summaries"] = threshold_summaries

    return summary

def query_cov_for_region(bam_file, sample_id, region_bed, region):
    """
    Run mosdepth for a given BAM and BED file, then normalize the
    mosdepth regions output to the old mosdepth-stat-like format expected
    by open_case_cov() and assign_probes():

        chrom\tstart\tend\tmean_coverage

    mosdepth creates <prefix>.regions.bed.gz. The last column is the mean
    coverage for the region. If the BED contains a name column, mosdepth
    keeps it before the coverage column, so we always use fields[0:3]
    plus fields[-1].
    """
    prefix = f"{sample_id}.{region}.mosdepth"
    mosdepth_regions_gz = f"{prefix}.regions.bed.gz"
    normalized_cov = f"{sample_id}.{region}.cov"

    if not bed_has_regions(region_bed):
        logging.debug(f"No regions in {region_bed}; writing empty coverage file {normalized_cov}")
        with open(normalized_cov, "w", encoding="utf-8"):
            pass
        return normalized_cov

    logging.debug(f"Calculating {region} coverage for {region_bed} with mosdepth ..")
    command = [
        "mosdepth",
        "--no-per-base",
        "--by", region_bed,
        prefix,
        bam_file,
    ]
    subprocess.run(command, stderr=subprocess.PIPE, text=True, check=True)

    with gzip.open(mosdepth_regions_gz, "rt") as infile, open(normalized_cov, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 4:
                continue
            chrom, start, end = fields[0], fields[1], fields[2]
            mean_cov = fields[-1]
            outfile.write(f"{chrom}\t{start}\t{end}\t{mean_cov}\n")

    logging.debug(f"Done, wrote normalized coverage to {normalized_cov}")
    return normalized_cov

def bed_has_regions(region_bed):
    with open(region_bed, "r", encoding="utf-8") as file:
        for line in file:
            if line.strip() and not line.startswith("#"):
                return True
    return False

def gtf_info_field(annotation_info:list):
    anno_info_list = annotation_info.split(';')
    #trailing ;
    anno_info_list.pop()
    tags = []
    anno_info_dict = {}
    for anno in anno_info_list:
        anno = anno.lstrip()
        anno = anno.split(' ')
        # tag field is special
        if anno[0] == "tag":
            tags.append(" ".join(anno[1:]).replace('"', '').replace("'", ""))
        else:
            anno_info_dict[anno[0]] = " ".join(anno[1:]).replace('"', '').replace("'", "")
    anno_info_dict['tags'] = tags
    return anno_info_dict

def intervals_overlap(start, end, other_start, other_end):
    return start < other_end and end > other_start

def is_y_chromosome(chromosome):
    normalized = str(chromosome).lower()
    if normalized.startswith("chr"):
        normalized = normalized[3:]
    return normalized == "y"

def prune_y_chromosome_regions(gene_gtf, exon_to_gene, cds_to_gene, probe_mapper):
    for region_map in (exon_to_gene, cds_to_gene, probe_mapper):
        for region_name in list(region_map.keys()):
            chromosome = region_name.split("_", 1)[0]
            if is_y_chromosome(chromosome):
                del region_map[region_name]

    for gene in list(gene_gtf.keys()):
        gene_record = gene_gtf[gene]
        for region_type in ("exons", "CDS", "probes"):
            for region_name, region in list(gene_record.get(region_type, {}).items()):
                chromosome = region.get("chr", region_name.split("_", 1)[0])
                if is_y_chromosome(chromosome):
                    del gene_record[region_type][region_name]

        transcript = gene_record.get("transcript", {})
        if is_y_chromosome(transcript.get("chr", "")):
            del gene_gtf[gene]

    return gene_gtf

def filter_y_chromosome_bed(bed_file, output_prefix, label):
    filtered_bed = f"{output_prefix}.{label}.no_y.bed"
    with open(bed_file, "r", encoding="utf-8") as infile, open(
        filtered_bed,
        "w",
        encoding="utf-8",
    ) as outfile:
        for line in infile:
            if line.startswith("#") or not line.strip():
                outfile.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if fields and is_y_chromosome(fields[0]):
                continue
            outfile.write(line)
    return filtered_bed

def assign_probes(gene_gtf,probe_mapper,file_path):
    logging.debug(f"Reading probe cov and assigning to gene ..")
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            line = line.split('\t')
            chr = line[0]
            start = int(line[1])
            end = int(line[2])
            cov = line[3]
            gene = None
            probe_name = f"{line[0]}_{line[1]}_{line[2]}"
            for entry in probe_mapper:
                parts = entry.split("_")
                gchr = parts[0]
                gstart = int(parts[1])
                gend = int(parts[2])
                #print(f"{gend}:{end}  {gchr}={chr}  {gstart}:{start}")
                if gchr == chr and intervals_overlap(start, end, gstart, gend):
                    gene = probe_mapper[entry]["gene_name"]
            if gene:
                gene_gtf[gene]["probes"][probe_name] = { 'cov': cov, 'chr':chr, 'start': start, 'end': end}
                gene_gtf[gene]["covered_by_panel"] = True
    logging.debug(f"Done")
    return gene_gtf

def open_case_cov(gene_gtf,translate_dict,file_path,region):
    logging.debug(f"Reading {region} cov and assigning to gene ..")
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            line = line.split('\t')
            cov = line[3]
            probe_name = f"{line[0]}_{line[1]}_{line[2]}"
            gene = translate_dict.get(probe_name,None)
            if gene:
                gene_gtf[gene][region][probe_name]['cov'] = cov
                gene_gtf[gene]["covered_by_panel"] = True
    logging.debug(f"Done")
    return gene_gtf

def gtf_start_to_bed_start(start):
    return max(int(start) - 1, 0)

def read_mane_gtf(mane_gtf: str, gene_list=None, filtered_gtf=None):
    logging.debug(f"Reading gtf: {mane_gtf} ..")
    gene_gtf = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    exon_to_gene = defaultdict(str)
    cds_to_gene = defaultdict(str)
    probe_mapper = defaultdict(dict)
    selected_genes = set(gene_list) if gene_list is not None else None
    if selected_genes is not None and filtered_gtf is None:
        filtered_gtf = "filtered.gtf"
    wgs_gtf = filtered_gtf
    filtered_out = open(wgs_gtf, 'w') if selected_genes is not None else None
    with open(mane_gtf, 'r') as file:
        for line in file:
            original = line
            line = line.rstrip()
            line = line.split('\t')
            if len(line) < 9:
                continue
            chr = line[0]
            start = line[3]
            end = line[4]
            annotation_type = line[2]
            annotation_info = line[8]
            annotation = gtf_info_field(annotation_info)
            gene_name = annotation.get('gene_name',None)
            tags = annotation.get('tags',[])
            # only keep MANE Select transcripts
            if 'MANE_Select' not in tags:
                continue
            if gene_name is None:
                continue
            if selected_genes is not None and gene_name not in selected_genes:
                continue
            else:
                if selected_genes is not None:
                    filtered_out.write(original)
                    gene_gtf[gene_name]["covered_by_panel"] = True
                else:
                    gene_gtf[gene_name]["covered_by_panel"] = False
            bed_start = str(gtf_start_to_bed_start(start))
            if annotation_type == "exon":
                exon          = {}
                exon['chr']   = chr
                exon['start'] = bed_start
                exon['end']   = end
                exon['nbr']   = annotation['exon_number']
                gene_gtf[annotation['gene_name']]['exons'][f"{chr}_{bed_start}_{end}"] = exon
                exon_to_gene[f"{chr}_{bed_start}_{end}"] = annotation['gene_name']
            elif annotation_type == "CDS":
                cds          = {}
                cds['chr']   = chr
                cds['start'] = bed_start
                cds['end']   = end
                cds['nbr']   = annotation['exon_number']
                gene_gtf[annotation['gene_name']]['CDS'][f"{chr}_{bed_start}_{end}"] = cds
                cds_to_gene[f"{chr}_{bed_start}_{end}"] = annotation['gene_name']
            elif annotation_type == "transcript":
                ## catch probes in promoters
                start_slop     = max(gtf_start_to_bed_start(start) - 1000, 0)
                end_slop       = int(end) + 1000
                slopped_coords = f"{chr}_{start_slop}_{end_slop}"
                ## assign values
                probe_mapper[slopped_coords]['gene_name']     = gene_name
                probe_mapper[slopped_coords]['start']         = bed_start
                probe_mapper[slopped_coords]['end']           = end
                probe_mapper[slopped_coords]['chr']           = chr
                probe_mapper[slopped_coords]['transcript_id'] = annotation['transcript_id']
                
                gene_gtf[gene_name]['transcript'] = {"chr":chr,'start':bed_start,'transcript_id':annotation['transcript_id'], 'end': end}
            elif annotation_type == "three_prime_utr":
                pass
            elif annotation_type == "five_prime_utr":
                pass
    if filtered_out:
        filtered_out.close()
    logging.debug(f"Done")
    return gene_gtf,exon_to_gene,cds_to_gene,probe_mapper,wgs_gtf

def generate_bed_regions(panel_gtf, output_prefix, exclude_y_regions=False):
    """
    generates bed-regions from gtf as input for mosdepth stats
    This is used both in panels and WGS
    """
    logging.debug(f"generating exon and CDS input beds")
    cds_file = f"{output_prefix}.cds.bed"
    exon_file = f"{output_prefix}.exons.bed"
    with open(panel_gtf, 'r') as file:
        with open(cds_file, "w") as cdsfile, open(exon_file, "w") as exonfile:
            for line in file:
                line = line.rstrip()
                line = line.split('\t')
                if len(line) < 5:
                    continue
                chr = line[0]
                if exclude_y_regions and is_y_chromosome(chr):
                    continue
                start = gtf_start_to_bed_start(line[3])
                end = line[4]
                annotation_type = line[2]
                if annotation_type == "exon":
                    exonfile.write(f"{chr}\t{start}\t{end}\n")
                elif annotation_type == "CDS":
                    cdsfile.write(f"{chr}\t{start}\t{end}\n")
    logging.debug(f"Done, created input beds {cds_file}, {exon_file}")
    return cds_file,exon_file

def intersect_gtf_with_design(mane_gtf,design_bed,output_prefix,exclude_y_regions=False):
    """
    use bedtools to pickout regions overlapping provided 
    design file
    Input:
        bed-file with probeplacement for gene-centric probes
        gtf containing mane-transcripts from encode
    Returns:
        input file for generate_bed_regions
    """
    logging.debug(f"Generating design gtf via bedtools")
    output_file = f"{output_prefix}.design.gtf"
    command = ['bedtools', 'intersect', '-a', mane_gtf, '-b', design_bed, '-u']
    with open(output_file, "w") as file:
        subprocess.run(command, stdout=file, stderr=subprocess.PIPE, text=True, check=True)
    cds_file,exon_file = generate_bed_regions(output_file, output_prefix, exclude_y_regions)
    return cds_file,exon_file

if __name__ == "__main__":
    ## logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )
    main()
