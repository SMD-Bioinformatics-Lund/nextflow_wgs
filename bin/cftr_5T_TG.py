#!/usr/bin/env python3

from pysam import VariantFile, AlignmentFile
import cmdvcf
import argparse
import itertools
from pprint import pprint
from collections import defaultdict

SCORE_ADJUSTMENT = 25
CENTER_7T = 117548632
PADDING_7T = 10
REF_T_REPS = 7
BAD_T_COUNT = 5
CENTER_11TG = 117548617
PADDING_11TG = 20
REF_TG_REPS = 11
BAD_TG_COUNT = 13
UPPER_RATIO_HETEROZYGOTE = 1
LOWER_RATIO_HETEROZYGOTE = 0.5

CFTR_REGION = { 
    "chrom": "7",
    "start": 117548601,
    "end" : 117548640,
    "ref": ['T', 'T', 'T', 'T', 'G', 'A', 'T', 'G', 'T', 'G', 'T', 'G', 'T', 'G', 'T', 'G',
            'T', 'G', 'T', 'G', 'T', 'G', 'T', 'G', 'T', 'G', 'T', 'G', 'T', 'T', 'T', 'T',
            'T', 'T', 'T', 'A', 'A', 'C', 'A', 'G']
}
"""
Logic of script:
0: check that there are any variants at all..
1: read the bam-file for the CFTR intron 8 region defined in CFTR_REGION
    * collect reads containing the 7T region
    * collect reads containing the 11TG region
    * check if allele counts could be considered causative BAD_T_COUNT/BAD_TG_COUNTS
    * if not, skip 2
2: read vcf-file and collect variants
    * check if variants could by themselves give rise to a 5T or 14TG repeats
    * if not by itself try combinations of variants
    * assign relevant variants with a higher score
    * prioritize finding 5T and if not, try 14TG (they are often correlated)
3: print new VCF with added scores if relevant, otherwise just print old vcf
"""

def main(args):
    
    vcf_object = VariantFile(args.input_vcf)
    bam = AlignmentFile(args.bam, "rb")
    
    variants = []
    variation_combo_str = None
    for var in vcf_object.fetch(CFTR_REGION['chrom'], CFTR_REGION['start'], CFTR_REGION["end"]):
        variants.append({
            "pos": var.pos,
            "ref": var.ref,
            "alt": var.alts[0],
            "id": f"{CFTR_REGION['chrom']}_{var.pos}_{var.ref}_{var.alts[0]}"
        })
    if len(variants) == 0:
        print_vcf(vcf_object,args.out_vcf,variation_combo_str)
        exit()

    print("Trying to find 5T variants")
    fragments_5t = get_overlapping_reads_from_bam(bam,CENTER_7T,PADDING_7T)
    count = defaultdict(int)
    for fragment in fragments_5t:
        max_count = count_consecutive_ts(fragment)
        count[max_count] = count[max_count] +1
    allele1,allele2 = find_one_or_two_alleles(count)

    if allele1 <= BAD_T_COUNT or allele2 <= BAD_T_COUNT: 
        variation_combo_str = combinatorics_to_find_T_causatives(variants,[allele1,allele2])
    else:
        print("Trying to find 13TG+ variants")
        if allele1 < BAD_TG_COUNT or allele2 < BAD_TG_COUNT:
            print("no bad TG counts either")
        fragments_11tg = get_overlapping_reads_from_bam(bam,CENTER_11TG,PADDING_11TG)
        count = defaultdict(int)
        for fragment in fragments_11tg:
            max_count = count_consecutive_tgs(fragment)
            count[max_count] = count[max_count] +1
        allele1,allele2 = find_one_or_two_alleles(count)
        variation_combo_str = combinatorics_to_find_TG_causatives(variants,[allele1,allele2])
    
    print_vcf(vcf_object,args.out_vcf,variation_combo_str)


def print_vcf(in_vcf,out_vcf,variation_combo_str):
    """
    """
    ids_of_variants_to_modify = []
    if variation_combo_str:
        ids_of_variants_to_modify = variation_combo_str.split('+')
    
    vcf_out = VariantFile(out_vcf, "w", header=in_vcf.header)
    vcf_out.close()
    with open(out_vcf, "a") as vcf_out:
        for var in in_vcf.fetch():
            var_dict = cmdvcf.parse_variant(var,in_vcf.header)
            var_id = f"{var_dict['CHROM']}_{var_dict['POS']}_{var_dict['REF']}_{var_dict['ALT']}"
            if var_id in ids_of_variants_to_modify:
                if var_dict['INFO'].get('RankScore'):
                    family_id = str(var_dict['INFO']['RankScore'].split(":")[0])
                    current_score = int(var_dict['INFO']['RankScore'].split(":")[1])
                    adjusted_score = current_score + SCORE_ADJUSTMENT
                    var_dict['INFO']['RankScore'] = f"{family_id}:{str(adjusted_score)}"
                else:
                    exit("Needs a scored VCF")
            vcf_str = cmdvcf.vcf_string(var_dict,in_vcf.header)
            vcf_out.write(vcf_str)
        

def parse_arguments():
    parser = argparse.ArgumentParser(description="Adds variantcaller specific penalties to variants with only one variant caller")

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
    parser.add_argument(
        "--bam",
        "-b",
        type=str,
        required=True,  
        help="name of sample bam-file"
    )
    args = parser.parse_args()
    return args

def count_consecutive_ts(string):
    count = 0
    max_count = 0
    for char in string:
        if char == 'T':
            count += 1
            max_count = max(max_count, count)
        else:
            count = 0  # reset when non-'T' encountered
    return max_count

def count_consecutive_tgs(s):
    count = 0
    max_count = 0
    i = 0

    while i < len(s) - 1:
        if s[i:i+2] == "TG":
            count += 1
            max_count = max(max_count, count)
            i += 2  # skip ahead since we've consumed one TG pair
        else:
            count = 0
            i += 1
    return max_count

def get_overlapping_reads_from_bam(bam,ref_pos,padding):
    reads = []
    for read in bam.fetch("7", ref_pos, ref_pos+1 ):
        if read.is_unmapped or read.is_secondary:
            continue
        # Find the read index that aligns to the reference position
        read_index = None
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if rpos == ref_pos:
                read_index = qpos
                break

        # If the read covers that position
        if read_index is not None:
            # Extract bases from that position onward
            seq_fragment = read.query_sequence[read_index-padding:read_index+padding]
            if len(seq_fragment) == 2*padding:
                reads.append(seq_fragment)
    return reads

def find_one_or_two_alleles(allele_counts):
    
    sorted_items = sorted(allele_counts.items(), key=lambda item: item[1], reverse=True)
    if len(sorted_items) == 0:
        return None,None
    elif len(sorted_items) == 1:
        return sorted_items[0][0], None
    else:
        top1 = sorted_items[0][1]
        top2 = sorted_items[1][1]

        ratio = top2 / top1 if top1 != 0 else 0

        if LOWER_RATIO_HETEROZYGOTE <= ratio <= UPPER_RATIO_HETEROZYGOTE:
            return sorted_items[0][0], sorted_items[1][0]
        else:
            return sorted_items[0][0], None

def apply_variants(base_seq, variants):
    """Apply a list of variants to a copy of the reference motif."""
    seq = list(base_seq)
    for var in variants:
        ref = var["ref"]
        alt = var["alt"]
        pos = var["pos"]
        relative_position = pos - CFTR_REGION['start']

        # Deletion
        if len(ref) > 1:
            del_length = len(ref) - len(alt)
            relative_position += 1
            for p in range(del_length):
                new = relative_position + p
                if new < len(seq):
                    seq[new] = "-"
        # Insertion
        elif len(alt) > 1:
            ins_length = len(alt) - len(ref)
            mutation = alt[1:]
            for p in range(ins_length):
                new = relative_position + p
                if new < len(seq):
                    seq[new] = f"{seq[new]}{mutation}"
        # SNP
        else:
            if 0 <= relative_position < len(seq):
                seq[relative_position] = alt
    return seq

def combinatorics_to_find_T_causatives(variants, alleles):
    # Iterate over all variant combinations
    for r in range(1, len(variants) + 1):
        for combo in itertools.combinations(variants, r):
            ref_motif_variant = apply_variants(CFTR_REGION['ref'], combo)
            variant_ids = [v['id'] for v in combo]
            combo_label = "+".join(variant_ids)

            variant_change_T_count = count_consecutive_ts(''.join(ref_motif_variant[-20:]))
            if variant_change_T_count <= BAD_T_COUNT and variant_change_T_count in alleles:
                print(f"Combination ({combo_label}) is causative ")
                return combo_label
    return None

def combinatorics_to_find_TG_causatives(variants,alleles):
    # Iterate over all variant combinations
    for r in range(1, len(variants) + 1):
        for combo in itertools.combinations(variants, r):
            ref_motif_variant = apply_variants(CFTR_REGION['ref'], combo)
            variant_ids = [v['id'] for v in combo]
            combo_label = "+".join(variant_ids)
            variant_change_TG_count = count_consecutive_tgs(''.join(ref_motif_variant))
            if variant_change_TG_count >= BAD_TG_COUNT and variant_change_TG_count in alleles:
                print(f"Combination ({combo_label}) is causative ")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
