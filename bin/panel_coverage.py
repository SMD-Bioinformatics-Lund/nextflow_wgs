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
    bam_file = args.bam
    mane_gtf = args.gtf
    # Panel
    if args.design_bed:
        design_bed = args.design_bed
        cds_file, exon_file = intersect_gtf_with_design(mane_gtf,design_bed)        
        gene_gtf,exon_to_gene,cds_to_gene,probe_mapper,wgs_gtf = read_mane_gtf(mane_gtf)
    # WGS
    else:
        gene_filter = args.gene_filter
        gene_list = read_gene_filter(gene_filter)
        gene_gtf,exon_to_gene,cds_to_gene,probe_mapper,wgs_gtf = read_mane_gtf(mane_gtf,gene_list)
        cds_file,exon_file = generate_bed_regions(wgs_gtf)

    # calculate coverage for sub-regions using mosdepth
    cds_cov   = query_cov_for_region(bam_file,sample_id,cds_file,"cds")
    exon_cov  = query_cov_for_region(bam_file,sample_id,exon_file,"exons")
    gene_gtf = open_case_cov(gene_gtf,cds_to_gene,cds_cov,"CDS")
    gene_gtf = open_case_cov(gene_gtf,exon_to_gene,exon_cov,"exons")
    if args.design_bed:
        probe_cov = query_cov_for_region(bam_file,sample_id,design_bed,"probes")
        gene_gtf = assign_probes(gene_gtf,probe_mapper,probe_cov)

    # only keep genes with some kind of coverage, probe or gene filtered for wgs    
    filtered_gene_gtf = {key: value for key, value in gene_gtf.items() if value.get("covered_by_panel") is True}
    
    mongolike = {}
    mongolike['genes'] = filtered_gene_gtf
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
        '-s', '--sample_id',
        type=str,
        required=True,
        help="name of sample"
    )
    parser.add_argument(
        '-f', '--gene_filter',
        type=str,
        help="list of genes to keep mutually exclusive to design_bed"
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        help="name of json output file. If not provided sample_id.cov.json"
    )
    args = parser.parse_args()
    if args.design_bed and args.gene_filter:
        exit("Cannot use both gene filter and design bed")
    elif args.design_bed is None and args.gene_filter is None:
        exit("Please provide either design_bed or gene_filter")
    return args

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
                if gchr == chr and start > gstart and end < gend:
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
    logging.debug(f"Done")
    return gene_gtf

def read_mane_gtf(mane_gtf: str, gene_list=None):
    logging.debug(f"Reading gtf: {mane_gtf} ..")
    gene_gtf = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    exon_to_gene = defaultdict(str)
    cds_to_gene = defaultdict(str)
    probe_mapper = defaultdict(dict)
    wgs_gtf = None
    if gene_list:
        wgs_gtf = "filtered.gtf"
    with open(mane_gtf, 'r') as file:
        for line in file:
            original = line
            line = line.rstrip()
            line = line.split('\t')
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
            else:
                if gene_list:
                    if gene_name in gene_list:
                        with open(wgs_gtf, 'a') as out:
                            out.write(original)
                        gene_gtf[gene_name]["covered_by_panel"] = True
                else:
                    gene_gtf[gene_name]["covered_by_panel"] = False
            if annotation_type == "exon":
                exon          = {}
                exon['chr']   = chr
                exon['start'] = start
                exon['end']   = end
                exon['nbr']   = annotation['exon_number']
                gene_gtf[annotation['gene_name']]['exons'][f"{chr}_{start}_{end}"] = exon
                exon_to_gene[f"{chr}_{start}_{end}"] = annotation['gene_name']
            elif annotation_type == "CDS":
                cds          = {}
                cds['chr']   = chr
                cds['start'] = start
                cds['end']   = end
                cds['nbr']   = annotation['exon_number']
                gene_gtf[annotation['gene_name']]['CDS'][f"{chr}_{start}_{end}"] = cds
                cds_to_gene[f"{chr}_{start}_{end}"] = annotation['gene_name']
            elif annotation_type == "transcript":
                ## catch probes in promoters
                start_slop     = int(start) - 1000
                end_slop       = int(end) + 1000
                slopped_coords = f"{chr}_{start_slop}_{end_slop}"
                ## assign values
                probe_mapper[slopped_coords]['gene_name']     = gene_name
                probe_mapper[slopped_coords]['start']         = start
                probe_mapper[slopped_coords]['end']           = end
                probe_mapper[slopped_coords]['chr']           = chr
                probe_mapper[slopped_coords]['transcript_id'] = annotation['transcript_id']
                
                gene_gtf[gene_name]['transcript'] = {"chr":chr,'start':start,'transcript_id':annotation['transcript_id'], 'end': end}
            elif annotation_type == "three_prime_utr":
                pass
            elif annotation_type == "five_prime_utr":
                pass
    logging.debug(f"Done")
    return gene_gtf,exon_to_gene,cds_to_gene,probe_mapper,wgs_gtf

def generate_bed_regions(panel_gtf):
    """
    generates bed-regions from gtf as input for mosdepth stats
    This is used both in panels and WGS
    """
    logging.debug(f"generating exon and CDS input beds")
    cds_file = "cds.bed"
    exon_file = "exons.bed"
    with open(panel_gtf, 'r') as file:
        with open(cds_file, "w") as cdsfile, open(exon_file, "w") as exonfile:
            for line in file:
                line = line.rstrip()
                line = line.split('\t')
                chr = line[0]
                start = line[3]
                end = line[4]
                annotation_type = line[2]
                if annotation_type == "exon":
                    exonfile.write(f"{chr}\t{start}\t{end}\n")
                elif annotation_type == "CDS":
                    cdsfile.write(f"{chr}\t{start}\t{end}\n")
    logging.debug(f"Done, created input beds {cds_file}, {exon_file}")
    return cds_file,exon_file

def intersect_gtf_with_design(mane_gtf,design_bed):
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
    output_file = 'design.gtf'
    command = ['bedtools', 'intersect', '-a', mane_gtf, '-b', design_bed, '-u']
    with open(output_file, "w") as file:
        subprocess.run(command, stdout=file, stderr=subprocess.PIPE, text=True, check=True)
    cds_file,exon_file = generate_bed_regions("design.gtf")
    return cds_file,exon_file

def read_gene_filter(gene_file):
    """
    reads genefilter into list
    """
    gene_list = []
    with open(gene_file, 'r') as file:
        for line in file:
            line = line.rstrip()
            gene_list.append(line)
    return gene_list

if __name__ == "__main__":
    ## logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
    )
    main()
