#!/usr/bin/env python3

import json
import argparse
from typing import Dict


def main():
    parser = argparse.ArgumentParser(description="Convert VCF genotypes into a JSON file.")
    parser.add_argument('genotype_file', type=str, help="Path to the VCF file")
    parser.add_argument('output_file', type=str, help="Path to save the output JSON file")
    
    args = parser.parse_args()

    genotype_data = parse_genotype_file(args.genotype_file)
    idsnp_json_format = {
        "is_paired_sample" : False,
        "id_snp_genotypes" : genotype_data
    }
    
    with open(args.output_file, 'w') as output_file:
        json.dump(idsnp_json_format, output_file, indent=4)
    
    
def parse_genotype_file(genotype_file_path: str) -> Dict[str, str]:
    genotype_dict = {}
    with open(genotype_file_path, 'r') as genotype_file:
        for line in genotype_file:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            var_id = parts[2]
            sample_genotype = parts[9]
            genotype = sample_genotype.split(":")[0]
            genotype_dict[var_id] = genotype
    return genotype_dict


if __name__ == '__main__':
    main()
