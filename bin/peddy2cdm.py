#!/usr/bin/env python3
import csv
import json

def load_sex_check(file_path):
    sex_data = {}

    with open(file_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row["sample_id"]
            sex_ok = row["error"].lower() == "false"
            sex_data[sample_id] = sex_ok

    return sex_data


def load_ped_check(file_path):
    ped_rows = []
    with open(file_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            ped_rows.append(row)
    return ped_rows


def evaluate_kinship(ped_rows, proband, father, mother):
    """
    return mother father errors
    """
    rel_status = {}

    for row in ped_rows:
        if row['sample_a'] == proband and row['sample_b'] == mother:
            rel_status['mother'] = is_correct(row)
            rel_status['mother_id'] = mother
        if row['sample_a'] == proband and row['sample_b'] == father:
            rel_status['father'] = is_correct(row)
            rel_status['father_id'] = father

    return rel_status

def is_correct(row):
    parent_error = row["parent_error"].lower() == "true"
    is_correct = not parent_error
    return is_correct
    

def build_per_sample_outputs(sex_data, ped_rows, proband, father, mother):
    is_trio = len(ped_rows) > 0

    results = {}

    if is_trio:
        rel_status = evaluate_kinship(ped_rows, proband, father, mother)

    for sample_id, sex_ok in sex_data.items():
        result = {
            "trio": is_trio,
            "sex": sex_ok
        }

        if is_trio:
            if sample_id == proband:
                result["kinship"] = rel_status
            elif sample_id in (father, mother):
                role = "father" if sample_id == father else "mother"
                result["kinship"] = {"proband": rel_status.get(role, False), "proband_id" : proband}

        results[sample_id] = result

    return results


def write_individual_jsons(results):

    for sample_id, data in results.items():
        filename = f"{sample_id}.json"
        with open(filename, "w") as f:
            json.dump(data, f, indent=2)


def main(ped_file, sex_file, proband, father, mother):
    sex_data = load_sex_check(sex_file)
    ped_rows = load_ped_check(ped_file)

    results = build_per_sample_outputs(sex_data, ped_rows, proband, father, mother)

    write_individual_jsons(results)

    print(f"Wrote {len(results)} JSON files")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process PED and SEX check files")
    parser.add_argument("--ped", required=True, help="id.ped_check.csv")
    parser.add_argument("--sex", required=True, help="id.sex_check.csv")
    parser.add_argument("--proband", required=True, help="proband-id")
    parser.add_argument("--mother", default=None, help="mother-id")
    parser.add_argument("--father", default=None, help="father-id")

    args = parser.parse_args()

    main(args.ped, args.sex, args.proband, args.father, args.mother)