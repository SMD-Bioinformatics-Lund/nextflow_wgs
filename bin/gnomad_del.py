#!/usr/bin/env python3
"""Prepare deletion intervals for bcftools and restore annotations to the original VCF."""

import argparse
import gzip
import re
import sys
from urllib.parse import quote, unquote


FIELDS = ("gnomad", "gnomad_AF", "gnomad_FREQ_HOMALT", "gnomad_SVTYPE")
MATCHES = "GNOMAD_DEL_MATCHES"
MAX_RECORDS = 50  # VEP's default custom-annotation limit, in reference source order.


def open_vcf(path):
    return gzip.open(str(path), "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def records(path):
    with open_vcf(path) as handle:
        row = 0
        for line in handle:
            if line.startswith("#"):
                continue
            row += 1
            columns = line.rstrip("\n").split("\t")
            info = dict(entry.split("=", 1) if "=" in entry else (entry, None)
                        for entry in columns[7].split(";"))
            yield row, columns, info


def schema(path):
    with open_vcf(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ,"):
                return line.split("Format: ", 1)[1].split('"', 1)[0].split("|")
            if line.startswith("#CHROM"):
                return None


def chromosome(name):
    name = name[3:] if name.startswith("chr") else name
    return "MT" if name == "M" else name


def deletion_span(columns, info):
    """Return an anchor-free, zero-based half-open deletion interval."""
    pos, ref, alt = int(columns[1]), columns[3], columns[4]
    if alt in ("<DEL>", "<CN0>", "<CN=0>") or (
            info.get("SVTYPE") == "DEL" and not re.fullmatch(r"[ACGT]+", ref + alt) and "," not in alt):
        length = info.get("SVLEN", ".").split(",")[0]
        end = pos + abs(int(length)) if length != "." else int(info.get("END", pos))
        return (pos, end) if end > pos else None
    if not re.fullmatch(r"[ACGTN]+", ref + alt, re.I):
        return None
    ref, alt = ref.upper(), alt.upper()
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    ref, alt = ref[prefix:], alt[prefix:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    start = pos - 1 + prefix
    return (start, start + len(ref)) if ref and not alt else None


def reference(vcf, output):
    destination = gzip.open(output, "wt") if str(output).endswith(".gz") else open(output, "w")
    with destination as handle:
        for row, columns, info in records(vcf):
            # These ALT classes correspond to VEP's ordinary 'deletion' class.
            # Mobile-element deletions have a different class and are excluded.
            if columns[4] not in ("<DEL>", "<CN0>", "<CN=0>"):
                continue
            span = deletion_span(columns, info)
            if span is None:
                raise ValueError("Incomplete reference deletion: " + columns[2])
            name = columns[2].split(";", 1)[0]
            if name in ("", ".", "0"):
                name = "{}:{}-{}".format(columns[0], span[0] + 1, span[1])
            values = []
            for key in ("AF", "FREQ_HOMALT", "SVTYPE"):
                value = info.get(key, "")
                value = "1" if value is None else value
                parts = value.split(",")
                values.append(parts[1] if len(parts) == 2 else value)  # VEP's REF/ALT-value handling.
            # One payload per match keeps missing fields aligned while bcftools
            # appends hits. Source row and coordinates preserve VEP order and
            # allow an exact integer check at the 70% boundary.
            payload = "|".join(quote(str(value), safe="") for value in (row, *span, name, *values))
            handle.write("{}\t{}\t{}\t{}\n".format(chromosome(columns[0]), *span, payload))


def query(vcf, output, synonyms=None):
    fields = schema(vcf)
    class_index = fields.index("VARIANT_CLASS") if fields else None
    aliases = {}
    if synonyms:
        with open(synonyms) as handle:
            for line in handle:
                names = [chromosome(name) for name in line.split()] if not line.startswith("#") else []
                if names:
                    aliases.update((name, names[0]) for name in names)
    intervals = []
    unclassified = 0
    for row, columns, info in records(vcf):
        if fields:
            if "CSQ" not in info:
                continue  # Keep VEP's skipped records unannotated.
            classes = {entry.split("|")[class_index] for entry in info["CSQ"].split(",")}
            if "deletion" not in classes:
                continue
            if classes != {"deletion"}:
                raise ValueError("Multiple VEP classes for " + columns[2])
        span = deletion_span(columns, info)
        if span is None:
            if fields:
                raise ValueError("Cannot represent deletion as an interval: {}:{} ({})".format(*columns[:3]))
            unclassified += info.get("SVTYPE") == "DEL"
            continue
        chrom = chromosome(columns[0])
        intervals.append((aliases.get(chrom, chrom), span[0], span[1], row))
    with open(output, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chrom in sorted({item[0] for item in intervals}):
            handle.write("##contig=<ID={}>\n".format(chrom))
        handle.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Matching interval end">\n')
        handle.write('##INFO=<ID={},Number=.,Type=String,Description="Internal gnomAD deletion matches">\n'.format(MATCHES))
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, start, end, row in sorted(intervals):
            handle.write("{}\t{}\t{}\tN\t<DEL>\t.\tPASS\tEND={}\n".format(chrom, start + 1, row, end))
    if unclassified:
        print("Left {} complex DEL records unchanged; VEP classification is needed.".format(unclassified), file=sys.stderr)


def restore(vcf, matches_vcf, output):
    annotations = {}
    for _, columns, info in records(matches_vcf):
        start, end = int(columns[1]) - 1, int(info["END"])
        hits = []
        for entry in (info.get(MATCHES) or "").split(","):
            if not entry:
                continue
            row, rstart, rend, *values = [unquote(value) for value in entry.split("|")]
            rstart, rend = int(rstart), int(rend)
            overlap = min(end, rend) - max(start, rstart)
            if 10 * overlap >= 7 * (end - start) and 10 * overlap >= 7 * (rend - rstart):
                hits.append((int(row), values))
        selected = [values for _, values in sorted(hits)[:MAX_RECORDS]]
        values = ["&".join(hit[i] for hit in selected if hit[i] != "") for i in range(4)]
        if len(hits) > MAX_RECORDS:
            values[0] += "&..."
        annotations[int(columns[2])] = [
            re.sub(r"\s+", "_", value.replace(",", "&").replace(";", "%3B")).replace("|", "&")
            if value != "-" else "" for value in values]
    fields = schema(vcf)
    indexes = [fields.index(field) for field in FIELDS] if fields else []
    row = 0
    with open_vcf(vcf) as source, open(output, "w") as target:
        for line in source:
            if line.startswith("#"):
                if line.startswith("#CHROM") and fields is None:
                    for field in FIELDS:
                        target.write('##INFO=<ID={},Number=.,Type=String,Description="gnomAD deletion annotation">\n'.format(field))
                target.write(line)
                continue
            row += 1
            if row not in annotations:
                target.write(line)
                continue
            columns = line.rstrip("\n").split("\t")
            info = columns[7].split(";")
            if fields:
                for i, entry in enumerate(info):
                    if entry.startswith("CSQ="):
                        consequences = []
                        for consequence in entry[4:].split(","):
                            parts = consequence.split("|")
                            for index, value in zip(indexes, annotations[row]):
                                parts[index] = value
                            consequences.append("|".join(parts))
                        info[i] = "CSQ=" + ",".join(consequences)
            else:
                info = [entry for entry in info if entry.split("=", 1)[0] not in FIELDS and entry != "."]
                info.extend("{}={}".format(field, value) for field, value in zip(FIELDS, annotations[row]) if value)
            columns[7] = ";".join(info) or "."
            target.write("\t".join(columns) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("reference", "query", "restore"))
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--synonyms")
    parser.add_argument("--matches")
    args = parser.parse_args()
    if args.mode == "reference":
        reference(args.vcf, args.output)
    elif args.mode == "query":
        query(args.vcf, args.output, args.synonyms)
    else:
        if not args.matches:
            parser.error("restore requires --matches")
        restore(args.vcf, args.matches, args.output)


if __name__ == "__main__":
    main()
