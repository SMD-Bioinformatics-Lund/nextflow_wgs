#!/usr/bin/env python3
import csv
import json
import argparse
from datetime import datetime, timezone

"""
This script takes output files from peddy and creates json input for register_sample.py --peddy
Needs sex_check.csv and ped_check.csv and sample ids + sequencing run ids

peddy2cdm.py --sex sex_check.csv --ped ped_check.csv --sample s1:run_idX

In addition the script also creates singaling files for middleman to load json into CDM
"""

def load_sex_check(file_path, sample):
    """
    read sex_check csv from peddy
    return given sex and if correct
    """
    sex_ok = True
    found = False
    with open(file_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row["sample_id"]
            if sample_id == sample:
                found = True
                sex_ok = row["error"].lower() == "false"
                ped_sex = row["ped_sex"]
    if not found:
        exit(f"{sample} is not in pedigree")
    return sex_ok, ped_sex


def load_ped_check(file_path):
    ped_rows = []
    with open(file_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ped_rows.append(row)
    return ped_rows


def evaluate_kinship(ped_rows, sample_id, samples_dict):
    """
    returns all matchings to other samples in trio
    """
    rel_status = {}

    for row in ped_rows:
        if row["sample_a"] == sample_id:
            other = row["sample_b"]
            rel_status[other] = {
                "is_correct": relation_is_correct(row),
                "sequencing_run": samples_dict.get(other),
            }

        if row["sample_b"] == sample_id:
            other = row["sample_a"]
            rel_status[other] = {
                "is_correct": relation_is_correct(row),
                "sequencing_run": samples_dict.get(other),
            }

    return rel_status


def relation_is_correct(row):
    parent_error = row["parent_error"].lower() == "true"
    is_correct = not parent_error
    return is_correct


def build_per_sample_outputs(sex_ok, ped_sex, ped_rows, sample, samples_dict):
    """
    for a given sample, check sex correctness and relationship status for others in family
    """
    result = {}

    is_trio = len(ped_rows) > 0
    if is_trio:
        rel_status = evaluate_kinship(ped_rows, sample, samples_dict)
        result["kinship"] = rel_status

    sex_check = {
        "is_correct_sex": sex_ok,
        "pedigree_sex": ped_sex,
    }

    result["sex_check"] = sex_check
    result["analysis_date"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    return result


def write_individual_jsons(results, output):

    filename = output
    with open(filename, "w") as f:
        json.dump(results, f, indent=2)


def write_cdm_load(sample, cdmassay, output_file, samples_dict, results_dir):

    seqrun = samples_dict.get(sample)
    filename = f"{sample}.peddy2cdm"
    with open(filename, "w") as f:
        f.write(
            f"--sequencing-run {seqrun} --assay {cdmassay} --sample-id {sample} --peddy {results_dir}/{output_file}"
        )


def parse_samples(sample_args):
    """
    expects  s1:run_idX for single samples
    and s1:run_idX&s2:run_idX&s3:run_idX for trios
    Returns:
    {
        "sample1": "seq1",
        "sample2": "seq2"
    }
    """
    samples = {}

    for part in sample_args:
        if ":" in part:
            sample, seqrun = part.split(":", 1)
        else:
            exit("Need format s1:run_idX")

        samples[sample] = seqrun

    return samples

def main(ped_file, sex_file, sample_arg, cdmassay, results_dir):
    ped_rows = load_ped_check(ped_file)

    samples_dict = parse_samples(sample_arg)

    for sample in samples_dict:
        sex_ok, ped_sex = load_sex_check(sex_file, sample)

        result = build_per_sample_outputs(
            sex_ok, ped_sex, ped_rows, sample, samples_dict
        )

        output_file = f"{sample}_peddy.json"

        write_individual_jsons(result, output_file)
        write_cdm_load(sample, cdmassay, output_file, samples_dict, results_dir)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process PED and SEX check files")
    parser.add_argument("--ped", required=True, help="id.ped_check.csv")
    parser.add_argument("--sex", required=True, help="id.sex_check.csv")
    parser.add_argument("--sample", required=True, action="append", help="sample:run_id (can be used multiple times)")
    parser.add_argument("--cdmassay", required=True, help="cdm assay of sample")
    parser.add_argument("--results_dir", required=True, help="cdm assay of sample")

    args = parser.parse_args()

    main(args.ped, args.sex, args.sample, args.cdmassay, args.results_dir)
