#!/usr/bin/env python3

from pysam import VariantFile
import cmdvcf
import argparse

def main(args: object):
    """
    Read VCF line by line, add genetic models to INFO-field
    """
    morbid_genes = read_omim_morbid(args.omim_morbid)
    vcf_object = VariantFile(args.input_vcf)
    vcf_object.header.info.add(
        "OMIM_morbid",
        ".",
        "String",
        "Overlaps gene defined in OMIM_morbid or OMIM_morbid_candidate"
    )
    vcf_out = VariantFile(args.out_vcf, "w", header=vcf_object.header)
    vcf_out.close()
    with open(args.out_vcf, "a") as vcf_out:
        for var in vcf_object.fetch():
            var_dict = cmdvcf.parse_variant(var,vcf_object.header)
            morbid = in_omim_morbid(var_dict,morbid_genes)
            if morbid:
                var_dict['INFO']['OMIM_morbid'] = "yes"
            vcf_str = cmdvcf.vcf_string(var_dict,vcf_object.header)
            vcf_out.write(vcf_str)


def read_omim_morbid(file_path:str) -> list:
    morbid_genes = []
    with open(file_path, 'r') as f:
        for line in f:
            morbid_genes.append(line.rstrip())
    return morbid_genes

def in_omim_morbid(var_dict:dict,morbid_genes:list) -> bool:
    """
    return annotated OMIM_morbid in INFO-field
    """
    morbid = False
    for trans in var_dict['INFO']['CSQ']:
        gene = trans['SYMBOL']
        if gene in morbid_genes:
            return True
    return morbid

        

def parse_arguments():
    parser = argparse.ArgumentParser(description="Adds OMIM_morbid=yes to variants overlapping genes in OMIM-AUTO genelist")

    parser.add_argument(
        "--input_vcf",
        "-i",
        type=str,
        required=True,  
        help="Path to the VCF that should be annotated"
    )
    parser.add_argument(
        "--omim_morbid",
        "-m",
        type=str,
        required=True,  
        help="Path omim morbid gene list"
    )
    parser.add_argument(
        "--out_vcf",
        "-o",
        type=str,
        required=True,  
        help="name of output.vcf"
    )
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_arguments()
    main(args)