#!/usr/bin/env python3
import csv
import json
from datetime import datetime, timezone

def load_sex_check(file_path, sample):
    sex_ok = True
    found = False
    with open(file_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row["sample_id"]
            if sample_id == sample:
                found = True
                sex_ok = row["error"].lower() == "false"
    if not found:
        exit(f"{sample} is not in pedigree")
    return sex_ok


def load_ped_check(file_path):
    ped_rows = []
    with open(file_path, newline='') as f:
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
        if row['sample_a'] == sample_id:
            other = row['sample_b']
            rel_status[other] = {
                "is_correct": is_correct(row),
                "sequencing_run": samples_dict.get(other)
            }

        if row['sample_b'] == sample_id:
            other = row['sample_a']
            rel_status[other] = {
                "is_correct": is_correct(row),
                "sequencing_run": samples_dict.get(other)
            }

    return rel_status

def is_correct(row):
    parent_error = row["parent_error"].lower() == "true"
    is_correct = not parent_error
    return is_correct
    

def build_per_sample_outputs(sex_ok, ped_rows, sample, samples_dict):
    is_trio = len(ped_rows) > 0

    if is_trio:
        rel_status = evaluate_kinship(ped_rows, sample, samples_dict)

    result = {
        "trio": is_trio,
        "sex": sex_ok,
        'analysis_date': datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }

    if is_trio:
        result['kinship'] = rel_status

    return result


def write_individual_jsons(results, output):


    filename = output
    with open(filename, "w") as f:
        json.dump(results, f, indent=2)

def write_cdm_load(sample,cdmassay,output_file,samples_dict,results_dir):

    seqrun = samples_dict.get(sample)
    filename = f"{sample}.peddy2cdm"
    with open(filename, "w") as f:
        f.write(f"--sequencing_run {seqrun} --assay {cdmassay} --sample-id {sample} --peddy {results_dir}/{output_file}")

def parse_samples(sample_arg):
    """
    Returns:
    {
        "sample1": "seq1",
        "sample2": "seq2"
    }
    """
    samples = {}
    parts = sample_arg.split("-")

    for part in parts:
        if ":" in part:
            sample, seqrun = part.split(":", 1)
        else:
            sample, seqrun = part, None

        samples[sample] = seqrun

    return samples

def main(ped_file, sex_file, sample_arg, cdmassay,results_dir):
    ped_rows = load_ped_check(ped_file)

    samples_dict = parse_samples(sample_arg)

    for sample in samples_dict:
        sex_ok = load_sex_check(sex_file, sample)

        result = build_per_sample_outputs(
            sex_ok,
            ped_rows,
            sample,
            samples_dict
        )


        output_file = f"{sample}.json"
        
        write_individual_jsons(result, output_file)
        write_cdm_load(sample,cdmassay,output_file,samples_dict,results_dir)
        


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process PED and SEX check files")
    parser.add_argument("--ped", required=True, help="id.ped_check.csv")
    parser.add_argument("--sex", required=True, help="id.sex_check.csv")
    parser.add_argument("--sample", required=True, help="id of sample")
    parser.add_argument("--cdmassay", required=True, help="cdm assay of sample")
    parser.add_argument("--results_dir", required=True, help="cdm assay of sample")

    args = parser.parse_args()

    main(args.ped, args.sex, args.sample, args.cdmassay, args.results_dir)