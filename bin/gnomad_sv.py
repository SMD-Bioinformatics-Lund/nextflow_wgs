#!/usr/bin/env python3
"""Prepare DEL/INS matching records and restore annotations to the original VCF."""

import argparse
from collections import Counter
from contextlib import ExitStack
import gzip
from pathlib import Path
import re
import sys
from urllib.parse import quote, unquote


FIELDS = ("gnomad", "gnomad_AF", "gnomad_FREQ_HOMALT", "gnomad_SVTYPE")
MATCHES = "GNOMAD_SV_MATCHES"
MODES = ("del", "ins_site", "ins_length")
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


def metadata_end(pos, info, kind):
    length = info.get("SVLEN", ".").split(",")[0]
    return pos + abs(int(length)) if length != "." else int(info.get("END", pos if kind == "DEL" else pos + 1))


def event(columns, info):
    """Return the type, matching mode, BED span and VEP's original scan bounds."""
    pos, ref, alt = int(columns[1]), columns[3], columns[4]
    if "," in alt:
        return None
    kind = {"<DEL>": "DEL", "<CN0>": "DEL", "<CN=0>": "DEL", "<INS>": "INS"}.get(alt)
    if not kind and not alt.startswith("<") and not re.fullmatch(r"[ACGT]+", ref + alt):
        kind = info.get("SVTYPE")
    if kind in ("DEL", "INS"):
        end = metadata_end(pos, info, kind)
        if kind == "DEL":
            return (kind, "del", pos, end, pos, end) if end > pos else None
        # VEP swaps reversed insertion boundaries before checking overlap.
        return kind, "ins_length", min(pos + 1, end) - 1, max(pos + 1, end), pos, end
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
    if ref and not alt:
        return "DEL", "del", start, start + len(ref), start, start + len(ref)
    if alt and not ref:
        # Preserve VEP's sequence-INS path: compare the insertion site, without
        # comparing inserted lengths. SVLEN does not affect this path.
        return "INS", "ins_site", start - 1, start, start, start
    return None


def reference(vcf, output):
    directory = Path(output)
    directory.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        handles = {mode: stack.enter_context(gzip.open(str(directory / (mode + ".bed.gz")), "wt")) for mode in MODES}
        for row, columns, info in records(vcf):
            kind = {"<DEL>": "DEL", "<CN0>": "DEL", "<CN=0>": "DEL", "<INS>": "INS"}.get(columns[4])
            if kind is None:
                continue
            pos = int(columns[1])
            end = metadata_end(pos, info, kind)
            # Original tabix bounds matter for legacy INS matching: VEP cannot
            # match a reference record that its original VCF index did not fetch.
            index_end = int(info["END"]) if info.get("END") else pos - 1 + len(columns[3])
            if kind == "DEL":
                if end <= pos:
                    raise ValueError("Incomplete reference deletion: " + columns[2])
                spans = {"del": (pos, end)}
            else:
                spans = {"ins_site": (pos - 1, pos - 1 + len(columns[3]))}
                if end > pos:
                    spans["ins_length"] = (pos, end)
            name = columns[2].split(";", 1)[0]
            if name in ("", ".", "0"):
                name = "{}:{}-{}".format(columns[0], pos + 1, end)
            values = []
            for key in ("AF", "FREQ_HOMALT", "SVTYPE"):
                value = info.get(key, "")
                value = "1" if value is None else value
                parts = value.split(",")
                values.append(parts[1] if len(parts) == 2 else value)  # VEP's REF/ALT-value handling.
            # One payload per match keeps missing fields aligned while bcftools
            # appends hits. Source row and coordinates preserve VEP order and
            # allow an exact integer check at the 70% boundary.
            for mode, span in spans.items():
                payload = "|".join(quote(str(value), safe="") for value in (row, *span, index_end, name, *values))
                handles[mode].write("{}\t{}\t{}\t{}\n".format(chromosome(columns[0]), *span, payload))


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
    intervals = {mode: [] for mode in MODES}
    unclassified = Counter()
    for row, columns, info in records(vcf):
        if fields:
            if "CSQ" not in info:
                continue  # Keep VEP's skipped records unannotated.
            classes = {entry.split("|")[class_index] for entry in info["CSQ"].split(",")}
            selected = classes.intersection({"deletion", "insertion"})
            if not selected:
                continue
            if len(classes) != 1:
                raise ValueError("Multiple VEP classes for " + columns[2])
        parsed = event(columns, info)
        if parsed is None:
            if fields:
                raise ValueError("Cannot represent DEL/INS for matching: {}:{} ({})".format(*columns[:3]))
            if info.get("SVTYPE") in ("DEL", "INS"):
                unclassified[info["SVTYPE"]] += 1
            continue
        kind, mode, start, end, seek_start, native_end = parsed
        expected_class = "deletion" if kind == "DEL" else "insertion"
        if fields and classes != {expected_class}:
            raise ValueError("VEP class disagrees with alleles for " + columns[2])
        chrom = chromosome(columns[0])
        intervals[mode].append((aliases.get(chrom, chrom), start, end, row, seek_start, native_end))
    for mode, rows in intervals.items():
        with open("{}.{}.intervals.vcf".format(output, mode), "w") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            for chrom in sorted({item[0] for item in rows}):
                handle.write("##contig=<ID={}>\n".format(chrom))
            for key, datatype in (("END", "Integer"), ("SEEK_START", "Integer"), ("QUERY_END", "Integer"), ("MATCH_MODE", "String")):
                handle.write('##INFO=<ID={},Number=1,Type={},Description="Internal matching coordinate or mode">\n'.format(key, datatype))
            handle.write('##INFO=<ID={},Number=.,Type=String,Description="Internal gnomAD matches">\n'.format(MATCHES))
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for chrom, start, end, row, seek_start, native_end in sorted(rows):
                kind = "DEL" if mode == "del" else "INS"
                handle.write("{}\t{}\t{}\tN\t<{}>\t.\tPASS\tEND={};MATCH_MODE={};SEEK_START={};QUERY_END={}\n".format(
                    chrom, start + 1, row, kind, end, mode, seek_start, native_end))
    if unclassified:
        print("Left complex records unchanged ({}); VEP classification is needed.".format(dict(unclassified)), file=sys.stderr)


def matched_values(matches_vcf):
    annotations = {}
    for _, columns, info in records(matches_vcf):
        start, end = int(columns[1]) - 1, int(info["END"])
        hits = []
        for entry in (info.get(MATCHES) or "").split(","):
            if not entry:
                continue
            row, rstart, rend, index_end, *values = [unquote(value) for value in entry.split("|")]
            rstart, rend = int(rstart), int(rend)
            overlap = min(end, rend) - max(start, rstart)
            mode = info["MATCH_MODE"]
            if mode != "ins_site" and (10 * overlap < 7 * (end - start) or 10 * overlap < 7 * (rend - rstart)):
                continue
            if mode.startswith("ins_"):
                reference_pos = rstart + 1 if mode == "ins_site" else rstart
                if int(index_end) <= int(info["SEEK_START"]) or reference_pos > int(info["QUERY_END"]):
                    continue
            hits.append((int(row), values))
        selected = [values for _, values in sorted(hits)[:MAX_RECORDS]]
        values = ["&".join(hit[i] for hit in selected if hit[i] != "") for i in range(4)]
        if len(hits) > MAX_RECORDS:
            values[0] += "&..."
        annotations[int(columns[2])] = [
            re.sub(r"\s+", "_", value.replace(",", "&").replace(";", "%3B")).replace("|", "&")
            if value != "-" else "" for value in values]
    return annotations


def restore(vcf, matches_vcfs, output):
    annotations = {}
    for path in matches_vcfs:
        annotations.update(matched_values(path))
    fields = schema(vcf)
    indexes = [fields.index(field) for field in FIELDS] if fields else []
    row = 0
    with open_vcf(vcf) as source, open(output, "w") as target:
        for line in source:
            if line.startswith("#"):
                if line.startswith("#CHROM") and fields is None:
                    for field in FIELDS:
                        target.write('##INFO=<ID={},Number=.,Type=String,Description="gnomAD DEL/INS annotation">\n'.format(field))
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
    parser.add_argument("--matches", nargs="+")
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
