#!/usr/bin/env python3

from pysam import VariantFile
import cmdvcf
import argparse
import csv
from collections import defaultdict
from pprint import pprint


KEYS_TO_KEEP = [
    "CHR2",
    "END",
    "EVIDENCE",
    "SVLEN",
    "SVTYPE",
    "AF",
    "N_HOMREF",
    "N_HET",
    "N_HOMALT",
    "FREQ_HOMREF",
    "FREQ_HET",
    "FREQ_HOMALT",
]

vcf_object = VariantFile("/fs1/resources/ref/hg38/annotation_dbs/gnomad_4.1/gnomad.v4.1.sv.sites.vcf.gz")
original_header = vcf_object.header

vcf_out = VariantFile("gnomad.v4.1.sv.sites.slim.vcf", "w", header=original_header)
vcf_out.close()
with open("gnomad.v4.1.sv.sites.slim.vcf", "a") as vcf_out:
    for var in vcf_object.fetch():
        new_info = {}
        var_dict = cmdvcf.parse_variant(var,vcf_object.header)
        var_dict['INFO']['END'] = var.stop
        # for key in KEYS_TO_KEEP:
        #     if key in var.info:
        #         new_info[key] = var.info[key]
        #     elif key == 'END':
        #         new_info[key] = var.stop
        for key in KEYS_TO_KEEP:
            if key in var_dict['INFO']:
                new_info[key] = var_dict['INFO'][key]
        var_dict['INFO'] = new_info
        vcf_str = cmdvcf.vcf_string(var_dict,vcf_object.header)
        vcf_out.write(vcf_str)
        
