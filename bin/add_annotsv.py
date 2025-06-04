#!/usr/bin/env python3

from pysam import VariantFile
import cmdvcf
import argparse
import csv
from collections import defaultdict

def main(args: object):
    """
    main function, reads VCF line by line and adds all AnnotSV annotations
    Updates vcf-header with the new annotations being added
    """
    vcf_object = VariantFile(args.input_vcf)
    original_header = vcf_object.header
    # TSV was loaded
    if args.annotsv_keys:
        annotsv_variants, annotsv_headers = read_annotsv_tsv(args.annotsv,args.annotsv_keys)
    else:
        annotsv_variants, annotsv_headers = read_annotsv_vcf(args.annotsv)
    
    vcf_header = join_headers(original_header,annotsv_headers)
    vcf_out = VariantFile(args.out_vcf, "w", header=vcf_header)
    vcf_out.close()
    with open(args.out_vcf, "a") as vcf_out:
        for var in vcf_object.fetch():
            # to fix END flag disapearing from pysam object
            # pysam eats up END from original, don't know why. (PERL vfc2.pm <3)
            var_dict = cmdvcf.parse_variant(var,vcf_object.header)
            svtype = var.info['SVTYPE']
            # AnnotSV treats TDUP as DUP
            if svtype == "TDUP":
                svtype = "DUP"
            simple_id = f"{var.chrom}_{var.pos}_{var.stop}_{svtype}"
            if simple_id in annotsv_variants:
                var_dict['INFO'] = add_annotsv_annotations(var_dict['INFO'],annotsv_variants[simple_id]['INFO'])
            vcf_str = cmdvcf.vcf_string(var_dict,vcf_object.header)
            vcf_out.write(vcf_str)
            
def add_annotsv_annotations(original_info: dict,annotation: dict):
    """
    Add the unique(new) keys from annotsv to original VCF INFO
    """
    for key in annotation:
        if key not in original_info:
            original_info[key] = annotation[key]
    return original_info

def join_headers(original_header: object, new_annotations: dict):
    """
    add new annotation headers to INFO
    """
    for annotation in new_annotations:
        if annotation['key'] not in original_header.info:
            original_header.info.add(
                annotation['key'],
                annotation['number'],
                annotation['type'],
                annotation['description']
            )
    return original_header

def read_annotsv_vcf(file_path: str):
    """
    read VCF from annotsv
    """
    annotsv_variants = {}
    vcf_in = VariantFile(file_path, "r")
    for var in vcf_in:
        var_dict = cmdvcf.parse_variant(var,vcf_in.header)
        # these are not needed, AnnotSV saves input VCF info-field and stores them weirdly
        del var_dict['INFO']['INFO']
        del var_dict['INFO']['FORMAT']
        simple_id_parts = var_dict['INFO']['AnnotSV_ID'].split("_")
        simple_id = f"{simple_id_parts[0]}_{simple_id_parts[1]}_{simple_id_parts[2]}_{simple_id_parts[3]}"
        annotsv_variants[simple_id] = var_dict
    annotsv_headers = remove_cohort_headers(vcf_in)
    vcf_in.close()
    return annotsv_variants,annotsv_headers

def remove_cohort_headers(vcf_in: object):
    """
    get VCF headers of info-field that gets added from AnnotSV

    Ignore [Ss]ample and Cohort specific, not valid VCF headers and not relevant for downsteam analysis
    """
    headers = []
    for key, record in vcf_in.header.info.items():
        if key not in ['ample', 'cohort']:
            info_dict = {}
            info_dict['key'] = key
            info_dict['description'] = record.description
            info_dict['type'] = record.type
            info_dict['number'] = record.number
            headers.append(info_dict)
    return headers

def read_annotsv_tsv(annotsv_tsv: str, keys: str):
    """
    Alternative to reading and merging a full VCF. Pick few fields
    from tsv.
    """
    add_keys = keys.split('-')
    annotsv_variants = defaultdict(dict)
    headers = []
    # Read TSV file into a dictionary
    with open(annotsv_tsv, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            info_dict = {}
            simple_id = f"{row['SV_chrom']}_{row['SV_start']}_{row['SV_end']}_{row['SV_type']}"
            for key in row:
                if key in add_keys:
                    if len(row[key]) > 0:
                        info_dict[key] = row[key]
            annotsv_variants[simple_id]['INFO'] = info_dict
    for key in add_keys:
        info_dict = {}
        info_dict['key'] = key
        info_dict['description'] = "Added annotation via annotsv TSV"
        info_dict['type'] = "String"
        info_dict['number'] = "."
        headers.append(info_dict)
    return annotsv_variants, headers


def parse_arguments():
    parser = argparse.ArgumentParser(description="Adds AnnotSV annotations to a VCF. Can use AnnotSV tsv or tsv converted to VCF via variantconvert")

    parser.add_argument(
        "--input_vcf",
        "-i",
        type=str,
        required=True,  
        help="Path to the VCF that should be annotated"
    )
    parser.add_argument(
        "--annotsv_vcf",
        "-a",
        type=str,
        help="Path to the AnnotSV VCF"
    )
    parser.add_argument(
        "--annotsv_tsv",
        "-t",
        type=str,
        help="Path to the AnnotSV TSV and what keys should be added to input VCF, generic VCF headers will be created. TSV:key1-key2-key3"
    )
    parser.add_argument(
        "--out_vcf",
        "-o",
        type=str,
        required=True,  
        help="name of output.vcf"
    )
    args = parser.parse_args()
    if args.annotsv_tsv:
        file_keys = args.annotsv_tsv.split(':')
        args.annotsv = file_keys[0]
        args.annotsv_keys = file_keys[1]
    elif args.annotsv_vcf:
        args.annotsv = args.annotsv_vcf
    else:
        exit("Need to supply either a annotsv-VCF or a annotsv TSV:key1-key2-key3")
    
    return args

if __name__ == "__main__":
    args = parse_arguments()
    main(args)