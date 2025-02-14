#!/usr/bin/env python3

from pysam import VariantFile
from pprint import pprint
import cmdvcf
import argparse

def main(args):
    """
    main function, reads VCF line by line and adds all AnnotSV annotations
    Updates vcf-header with the new annotations being added
    """
    annotsv_variants, annotsv_headers = read_annotsv_vcf(args.annotsv)
    vcf_object = VariantFile(args.input_vcf)
    
    vcf_header = join_headers(vcf_object.header,annotsv_headers)
    vcf_out = VariantFile(args.out_vcf, "w", header=vcf_header)
    vcf_out.close()
    for var in vcf_object.fetch():
        # to fix END flag disapearing from pysam object
        # pysam eats up END from original, don't know why. (PERL vfc2.pm <3)
        var_dict = cmdvcf.parse_variant(var,vcf_object.header)
        var_dict['INFO']['END'] = var.stop
        svtype = var.info['SVTYPE']
        # AnnotSV treats TDUP as DUP
        if svtype == "TDUP":
            svtype = "DUP"
        simple_id = f"{var.chrom}_{var.pos}_{svtype}"
        if simple_id in annotsv_variants:
            var_dict['INFO'] = add_annotsv_annotations(var_dict['INFO'],annotsv_variants[simple_id]['INFO'])
        vcf_str = cmdvcf.print_vcf_record(var_dict,vcf_object.header)
        with open(args.out_vcf, "a") as vcf_out:
            vcf_out.write(vcf_str)
        
def add_annotsv_annotations(original_info,annotation):
    """
    Add the unique(new) keys from annotsv to original VCF INFO
    """
    for key in annotation:
        if key not in original_info:
            original_info[key] = annotation[key]
    return original_info

def join_headers(original_header, new_annotations):
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

def read_annotsv_vcf(file_path):
    """
    read VCF from annotsv
    """
    vcf_in = VariantFile(file_path, "r")
    annotsv_variants = {}
    for var in vcf_in:
        var_dict = cmdvcf.parse_variant(var,vcf_in.header)
        # these are not needed, AnnotSV saves input VCF info-field and stores them weirdly
        del var_dict['INFO']['INFO']
        del var_dict['INFO']['FORMAT']
        simple_id_parts = var_dict['INFO']['AnnotSV_ID'].split("_")
        simple_id = f"{simple_id_parts[0]}_{simple_id_parts[1]}_{simple_id_parts[3]}"
        annotsv_variants[simple_id] = var_dict

    annotsv_headers = get_relevant_headers(vcf_in)
    return annotsv_variants,annotsv_headers

def get_relevant_headers(vcf_in):
    """
    get VCF headers of info-field that gets added from AnnotSV

    Ignore [Ss]ample and Cohort specific, not valid VCF headers and not relevant for downsteam analysis
    """
    headers = []
    for key, record in vcf_in.header.info.items():
        if 'ample' in key or 'cohort' in key:
            ...
        else:
            info_dict = {}
            info_dict['key'] = key
            info_dict['description'] = record.description
            info_dict['type'] = record.type
            info_dict['number'] = record.number
            headers.append(info_dict)
    return headers

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adds AnnotSV annotations to a VCF. Uses AnnotSV tsv that has been converted to VCF via variantconvert")

    parser.add_argument(
        "--input_vcf",
        "-i",
        type=str,
        required=True,  
        help="Path to the VCF that should be annotated"
    )
    parser.add_argument(
        "--annotsv",
        "-a",
        type=str,
        required=True,  
        help="Path to the AnnotSV VCF"
    )
    parser.add_argument(
        "--out_vcf",
        "-o",
        type=str,
        required=True,  
        help="name of output.vcf"
    )
    args = parser.parse_args()
    main(args)