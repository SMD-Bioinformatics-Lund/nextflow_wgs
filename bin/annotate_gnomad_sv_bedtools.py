#!/usr/bin/env python3
"""Populate VEP's gnomAD SV CSQ fields using same-class reciprocal BED overlap."""

import argparse
import gzip
import json
from pathlib import Path
import re
import subprocess
import tempfile


FIELDS = ("gnomad", "gnomad_AF", "gnomad_FREQ_HOMALT", "gnomad_SVTYPE")
MAX_RECORDS = 50  # Match VEP's default custom-annotation record limit.
SV_CLASSES = {
    "DEL": "deletion", "DUP": "duplication", "TDUP": "tandem_duplication",
    "INS": "insertion", "INV": "inversion", "CNV": "copy_number_variation",
    "TREP": "tandem_repeat", "BND": "chromosome_breakpoint",
}
for _kind, _name in (("DEL", "deletion"), ("INS", "insertion")):
    for _suffix, _prefix in (("ME", "mobile_element"), ("ALU", "Alu"),
                             ("HERV", "HERV"), ("LINE1", "LINE1"), ("SVA", "SVA")):
        SV_CLASSES[_kind + "_" + _suffix] = _prefix + "_" + _name


def open_vcf(path):
    return gzip.open(str(path), "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def records(path):
    with open_vcf(path) as handle:
        row = 0
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            row += 1
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 8:
                raise ValueError("{}: record {} has fewer than eight columns".format(path, row))
            info = dict(field.split("=", 1) if "=" in field else (field, None)
                        for field in columns[7].split(";"))
            yield row, columns, info


def csq_schema(line):
    match = re.search(r'^##INFO=<ID=CSQ,.*?Format: ([^"]+)', line)
    return match.group(1).split("|") if match else None


def chromosome_names(synonyms):
    def basic(name):
        name = name[3:] if name.startswith("chr") else name
        return "MT" if name == "M" else name

    aliases = {}
    with open(synonyms) as handle:
        for line in handle:
            names = [basic(name) for name in line.split()] if not line.startswith("#") else []
            if names:
                aliases.update((name, names[0]) for name in names)
    return lambda name: aliases.get(basic(name), basic(name))


def symbolic_class(alt):
    # Follow VEP 114's ALT-to-SO classes, including tandem/mobile-element types.
    # https://github.com/Ensembl/ensembl-vep/blob/release/114.2/modules/Bio/EnsEMBL/VEP/Parser.pm
    mobile = re.search(r'(INS|DEL):ME:?([A-Z0-9]+)?', alt, re.I)
    if mobile:
        element = (mobile.group(2) or "ME").upper()
        element = "LINE1" if element == "L1" else element
        element = element if element in ("ALU", "HERV", "LINE1", "SVA") else "ME"
        return SV_CLASSES[mobile.group(1).upper() + "_" + element]
    if "DUP:TANDEM" in alt:
        return SV_CLASSES["TDUP"]
    if "CNV:TR" in alt:
        return SV_CLASSES["TREP"]
    if re.search(r'CN=?[0-9]', alt):
        return SV_CLASSES["DEL" if re.fullmatch(r'<?CN=?0>?', alt) else
                          "DUP" if re.fullmatch(r'<?CN=?2>?', alt) else "CNV"]
    if "CNV" in alt:
        return SV_CLASSES["CNV"]
    if any(char in alt for char in "[]") or alt.startswith(".") or alt.endswith("."):
        return SV_CLASSES["BND"]
    return SV_CLASSES.get(alt.strip("<>").split(":", 1)[0])


def interval(columns, info):
    """Return an anchor-free half-open interval and a reference SO class."""
    pos, ref, alt = int(columns[1]), columns[3], columns[4]
    # VEP uses SVTYPE/SVLEN for alleles containing unknown bases, including
    # Manta's combined insertion/deletion records with N-padded REF alleles.
    structural = info.get("SVTYPE") and not re.fullmatch(r'[ACGT]+', ref + alt)
    if not structural and re.fullmatch(r'[ACGTN]+', ref, re.I) and re.fullmatch(r'[ACGTN]+', alt, re.I):
        ref, alt = ref.upper(), alt.upper()
        prefix = 0
        while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
            prefix += 1
        ref, alt, start = ref[prefix:], alt[prefix:], pos - 1 + prefix
        suffix = 0
        while suffix < min(len(ref), len(alt)) and ref[-suffix - 1] == alt[-suffix - 1]:
            suffix += 1
        if suffix:
            ref, alt = ref[:-suffix], alt[:-suffix]
        if ref and not alt:
            return start, start + len(ref), "deletion"
        if alt and not ref:
            # VEP's symbolic INS interval uses inserted length when SVLEN exists.
            return start, start + len(alt), "insertion"
        return None

    kind = symbolic_class(alt)
    if not kind and structural and "," not in alt:
        kind = symbolic_class(info["SVTYPE"])
    if not kind or kind == "chromosome_breakpoint":
        return None
    length = info.get("SVLEN", ".").split(",")[0]
    end = pos + abs(int(length)) if length != "." else int(info.get("END", pos + 1))
    if end <= pos:
        if "insertion" not in kind:
            return None
        return pos - 1, pos + 1, kind  # VEP's un-sized insertion boundary interval.
    return pos, end, kind


def prepare_beds(vcf, database, synonyms, query_bed, reference_bed):
    chrom = chromosome_names(synonyms)
    with open_vcf(vcf) as handle:
        schema = next((csq_schema(line) for line in handle if line.startswith("##INFO=<ID=CSQ,")), None)
    if schema is None or "VARIANT_CLASS" not in schema:
        raise ValueError("Expected VEP CSQ header with VARIANT_CLASS")
    class_index = schema.index("VARIANT_CLASS")
    with open(query_bed, "w") as output:
        for row, columns, info in records(vcf):
            if "CSQ" not in info:
                continue  # Do not rescue records VEP deliberately left unannotated.
            classes = {entry.split("|")[class_index] for entry in info["CSQ"].split(",")}
            if len(classes) != 1:
                raise ValueError("Record {} has multiple VEP variant classes".format(row))
            span = interval(columns, info)
            if span:
                output.write("{}\t{}\t{}\t{}\t0\t.\t{}\n".format(chrom(columns[0]), span[0], span[1], row, classes.pop()))
            elif classes.intersection(set(SV_CLASSES.values()) - {"chromosome_breakpoint"}):
                raise ValueError("Cannot represent SV record {} ({}) as a BED interval".format(row, columns[2]))
    with open(reference_bed, "w") as output:
        for row, columns, info in records(database):
            span = interval(columns, info)
            if span:
                # VEP reports the first ID, or source coordinates when it is absent.
                name = columns[2].split(";")[0]
                if name in ("", ".", "0"):
                    name = "{}:{}-{}".format(columns[0], span[0] + 1, span[1])
                data = []
                for field in ("AF", "FREQ_HOMALT", "SVTYPE"):
                    value = info.get(field, "")
                    value = "1" if value is None else value  # VEP reports INFO flags as 1.
                    parts = value.split(",")
                    if len(parts) == len(columns[4].split(",")) + 1:
                        value = ",".join(parts[1:])  # Discard a value for REF, as VEP does.
                    data.append(value)
                values = (chrom(columns[0]), span[0], span[1], row, "0", ".", name, span[2], json.dumps(data))
                output.write("\t".join(map(str, values)) + "\n")
    return schema


def find_matches(query_bed, reference_bed, overlaps):
    with open(overlaps, "w") as output:
        # bedtools uses single-precision fractions: at large SV sizes it can
        # round a true 70% match down, or a sub-70% match up. Get candidates
        # with a small tolerance, then enforce exactly 7/10 using integers.
        subprocess.run(["bedtools", "intersect", "-a", str(query_bed), "-b", str(reference_bed),
                        "-f", "0.699999", "-r", "-wa", "-wb"], stdout=output, check=True)
    matches = {}
    with open(overlaps) as handle:
        for line in handle:
            values = line.rstrip("\n").split("\t")
            q_start, q_end, r_start, r_end = map(int, (values[1], values[2], values[8], values[9]))
            overlap = min(q_end, r_end) - max(q_start, r_start)
            if (values[6] == values[14] and 10 * overlap >= 7 * (q_end - q_start)
                    and 10 * overlap >= 7 * (r_end - r_start)):
                fields = [values[13]] + json.loads(values[15])
                matches.setdefault(int(values[3]), []).append((int(values[10]), fields))
    return matches


def annotation_values(matches):
    ordered = sorted(matches)
    selected = [entry[1] for entry in ordered[:MAX_RECORDS]]
    # Preserve every reported match, so existing CSQ/max-AF scoring remains valid.
    # VEP omits absent INFO fields from each list, but keeps literal '.' values.
    values = ["&".join(entry[index] for entry in selected if entry[index] != "") for index in range(4)]
    if len(ordered) > MAX_RECORDS:
        values[0] += "&..."
    return [re.sub(r"\s+", "_", value.replace(",", "&").replace(";", "%3B")).replace("|", "&")
            if value != "-" else "" for value in values]


def write_annotations(vcf, output, schema, matches, database):
    fields = schema + [field for field in FIELDS if field not in schema]
    indexes = [fields.index(field) for field in FIELDS]
    header_ids = set()
    row = 0
    with open_vcf(vcf) as source, open(output, "w") as target:
        for line in source:
            if line.startswith("#"):
                info_header = re.match(r"##INFO=<ID=([^,]+),", line)
                if info_header:
                    header_ids.add(info_header.group(1))
                if csq_schema(line):
                    line = line.replace("Format: " + "|".join(schema), "Format: " + "|".join(fields), 1)
                if line.startswith(("##VEP-command-line=", "#CHROM")):
                    # Restore VEP's custom INFO headers. The configured database
                    # uses an absolute path, which VEP masks to [PATH]/basename.
                    for field in FIELDS:
                        if field not in header_ids:
                            description = (field[len("gnomad_"):] + " field from " if field != "gnomad" else "")
                            description += "[PATH]/" + Path(database).name
                            target.write('##INFO=<ID={},Number=.,Type=String,Description="{}">\n'.format(field, description))
                    header_ids.update(FIELDS)
                target.write(line)
                continue
            if not line.strip():
                continue
            row += 1
            columns = line.rstrip("\n").split("\t")
            values = annotation_values(matches.get(row, []))
            info = columns[7].split(";")
            for index, entry in enumerate(info):
                if entry.startswith("CSQ="):
                    consequences = []
                    for consequence in entry[4:].split(","):
                        parts = consequence.split("|")
                        parts += [""] * (len(fields) - len(parts))
                        for field, value in zip(indexes, values):
                            parts[field] = value
                        consequences.append("|".join(parts))
                    info[index] = "CSQ=" + ",".join(consequences)
            columns[7] = ";".join(info)
            target.write("\t".join(columns) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--synonyms", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="gnomad_sv_") as work:
        query, reference, overlaps = (Path(work) / name for name in ("query.bed", "reference.bed", "overlaps.tsv"))
        schema = prepare_beds(args.vcf, args.database, args.synonyms, query, reference)
        matches = find_matches(query, reference, overlaps)
        write_annotations(args.vcf, args.output, schema, matches, args.database)


if __name__ == "__main__":
    main()
