#!/usr/bin/env python3

from pysam import VariantFile
import cmdvcf
import argparse
from pprint import pprint
from geneticmodels import GM

MANTA_SIZE_FOR_PENALTY = 100000

def main(args: object):
    """
    Read VCF line by line, add genetic models to INFO-field
    """
    vcf_object = VariantFile(args.input_vcf)
    vcf_object.header.info.add(
        "GQC",
        ".",
        "String",
        "GATK genotype quality score when variant is only called by GATK"
    )
    vcf_object.header.info.add(
        "MANTAPENALTY",
        ".",
        "String",
        "big manta events not supported by gatk/cnvnator"
    )
    vcf_out = VariantFile(args.out_vcf, "w", header=vcf_object.header)
    vcf_out.close()
    with open(args.out_vcf, "a") as vcf_out:
        for var in vcf_object.fetch():
            var_dict = cmdvcf.parse_variant(var,vcf_object.header)
            varcall_set = var_dict['INFO'].get('set',None)
            if varcall_set is not None:
                varcallers_list = variant_called_by(varcall_set)
                # GATK only caller, return average genotype qual for variant
                if len(varcallers_list) == 1 and 'gatk' in varcallers_list:
                    var_dict['INFO']['GQC'] = modify_gatk(var_dict)
                # manta only caller, return bool if certain thingie
                elif len(varcallers_list) == 1 and 'manta' in varcallers_list:
                    mantapen = modify_manta(var_dict)
                    if mantapen is not None:
                        var_dict['INFO']['MANTAPENALTY'] = mantapen
            vcf_str = cmdvcf.vcf_string(var_dict,vcf_object.header)
            vcf_out.write(vcf_str)

def variant_called_by(varcallers: str):
    """
    checks what callers called variant
    """
    varcallers_list = []
    if 'manta' in varcallers:
        varcallers_list.append('manta')
    if 'gatk' in varcallers:
        varcallers_list.append('gatk')
    if 'tiddit' in varcallers:
        varcallers_list.append('tiddit')
    return varcallers_list

def modify_gatk(var_dict: dict):
    """
    calculate the average genotype quality of all calls for variant
    """
    genotype_qual_sum = 0
    number_of_qual_values = 0
    for sample in var_dict['GT']:
        genotype_qual = sample.get('QS',None)
        if genotype_qual is not None and genotype_qual != '.':
            genotype_qual_sum = genotype_qual_sum + int(genotype_qual)
            number_of_qual_values = number_of_qual_values + 1
    mean_qual_val = int(genotype_qual_sum/number_of_qual_values)
    return mean_qual_val

def modify_manta(var_dict: dict):
    """
    return stringified boolean if manta called a smaller variant
    than configured cutoff. Rankmodel expects string
    """
    mantapen = None
    svlen = var_dict['INFO'].get('SVLEN',None)
    if svlen is not None:
        if abs(var_dict['INFO']['SVLEN']) > MANTA_SIZE_FOR_PENALTY:
            mantapen = "1"
    return mantapen
        

def parse_arguments():
    parser = argparse.ArgumentParser(description="Adds genetic models to SV-variants, needs a matching pedigree file")

    parser.add_argument(
        "--input_vcf",
        "-i",
        type=str,
        required=True,  
        help="Path to the VCF that should be annotated"
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