#!/usr/bin/env python3

from pysam import VariantFile
import cmdvcf
import argparse
from geneticmodels import GM

def main(args: object):
    """
    Read VCF line by line, add genetic models to INFO-field
    """
    ped, individuals = read_ped(args.pedigree_file)
    vcf_object = VariantFile(args.input_vcf)
    vcf_object.header.info.add(
        "GeneticModel",
        ".",
        "String",
        "Genetic model for variant"
    )
    vcf_object.header.info.add(
        "HOMHEM",
        ".",
        "String",
        "If variant follows geneticmodel and is a homo- or hemi-zygous deletion"
    )
    vcf_out = VariantFile(args.out_vcf, "w", header=vcf_object.header)
    vcf_out.close()
    with open(args.out_vcf, "a") as vcf_out:
        for var in vcf_object.fetch():
            var_dict = cmdvcf.parse_variant(var,vcf_object.header)
            if len(ped.keys()) > 2:
                genetic_model = genetic_model_of_variant(var_dict,ped,individuals)
                var_dict['INFO']['GeneticModel'] = genetic_model
            if var_dict['INFO']['SVTYPE'] == 'DEL':
                hemi_homo_del = is_del_hemi_homo(var_dict,ped,individuals)
                if hemi_homo_del is not None:
                    var_dict['INFO']['HOMHEM'] = hemi_homo_del
            vcf_str = cmdvcf.vcf_string(var_dict,vcf_object.header)
            vcf_out.write(vcf_str)


def genotype_sum_per_individual(gt:dict):
    """
    get sum of genotype, return dict with sample-id and sum
    Hom alt = 2, hom ref = 0, het = 1
    """
    sample_gt_sum = {}
    for sample in gt:
        sample_id = sample.get('_sample_id')
        genotype = sample.get('GT')
        genotype = genotype.replace('.',"0")
        gt_list = genotype.split("/")
        gt_sum = sum(int(x) for x in gt_list)
        sample_gt_sum[sample_id] = gt_sum
    return sample_gt_sum

def genetic_model_of_variant(var_dict:dict,ped:dict,individuals:dict):
    """
    Returns a genetic model that is compliant with disease
    """
    sample_gt_sum = genotype_sum_per_individual(var_dict['GT'])
    chrom = var_dict['CHROM']
    x = 0
    if chrom == 'X':
        x = 1
    mother_gt   = sample_gt_sum[individuals['mother']]
    mother_phe  = ped[individuals['mother']]['PHENO']
    father_gt   = sample_gt_sum[individuals['father']]
    father_phe  = ped[individuals['father']]['PHENO']
    proband_gt  = sample_gt_sum[individuals['proband']]
    proband_sex = ped[individuals['proband']]['SEX']

    gm_search = f"{proband_sex}{proband_gt}{x}{father_gt}{father_phe}{mother_gt}{mother_phe}"
    # ignore homozygous ref proband GTs. Is this a problem?
    genetic_model = GM.get(gm_search,'NA')
    # at the moment rediculous genetic patterns including homozygous dominant traits are marked *, compound finder cannot handle these
    genetic_model = genetic_model.replace('*','')
        
    return genetic_model
    
def is_del_hemi_homo(var_dict:dict,ped,individuals):
    """
    checks if deletion is complete loss of genetic material
    Boys on X special, ladies on X maybe not needed to check
    """
    chrom = var_dict['CHROM']
    x = 0
    if chrom == 'X':
        x = 1
    sample_gt_sum = genotype_sum_per_individual(var_dict['GT'])
    proband_gt  = sample_gt_sum[individuals['proband']]
    proband_sex = ped[individuals['proband']]['SEX']
    is_del_hemi_homo = None
    if x:
        if int(proband_sex) == 1 and int(proband_gt) > 0:
            is_del_hemi_homo = 'LOSS'
        elif int(proband_sex) == 2 and int(proband_gt) == 2:
            is_del_hemi_homo = 'LOSS'
    elif int(proband_gt) == 2:
        is_del_hemi_homo = 'LOSS'
    return is_del_hemi_homo
        

def read_ped(pedfile: str):
    ped = {}
    individuals = {}
    with open(pedfile, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            ind = {
                "FATHER": fields[2],
                "MOTHER": fields[3],
                "SEX": fields[4],
                "PHENO": fields[5]
            }
            ped[fields[1]] = ind
            if fields[2] != "0" or fields[3] != "0":
                individuals["proband"] = fields[1]
                individuals["mother"] = fields[3]
                individuals["father"] = fields[2]
    if len(ped.keys()) == 1:
        for ind in ped:
            individuals['proband'] = ind
    return ped, individuals

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
        "--pedigree_file",
        "-p",
        type=str,
        required=True,  
        help="Path to the trio pedigree file"
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