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


def evaluate_kinship(ped_rows, sample_id):
    """
    returns all matchings to other samples in trio
    """
    rel_status = {}

    for row in ped_rows:
        if row['sample_a'] == sample_id:
            rel_status[row['sample_b']] = is_correct(row)
        if row['sample_b'] == sample_id:
            rel_status[row['sample_a']] = is_correct(row)

    return rel_status

def is_correct(row):
    parent_error = row["parent_error"].lower() == "true"
    is_correct = not parent_error
    return is_correct
    

def build_per_sample_outputs(sex_ok, ped_rows, sample):
    is_trio = len(ped_rows) > 0

    if is_trio:
        rel_status = evaluate_kinship(ped_rows, sample)

    result = {
        "trio": is_trio,
        "sex": sex_ok,
        'analysis_date' : datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }

    if is_trio:
        result['kinship'] = rel_status

    return result


def write_individual_jsons(results, output):


    filename = output
    with open(filename, "w") as f:
        json.dump(results, f, indent=2)


def main(ped_file, sex_file, sample, output):
    sex_ok = load_sex_check(sex_file,sample)
    ped_rows = load_ped_check(ped_file)

    result = build_per_sample_outputs(sex_ok, ped_rows, sample)

    write_individual_jsons(result, output)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process PED and SEX check files")
    parser.add_argument("--ped", required=True, help="id.ped_check.csv")
    parser.add_argument("--sex", required=True, help="id.sex_check.csv")
    parser.add_argument("--sample", required=True, help="id of sample")
    parser.add_argument("--output", default=None, help="id of sample")


    args = parser.parse_args()

    if args.output is None:
        args.output = f"{args.sample}.json"

    main(args.ped, args.sex, args.sample, args.output)