#!/usr/bin/env python3
"""Prepare SV matching records and add gnomAD annotations to VCF INFO fields."""

import argparse
from collections import Counter
from contextlib import ExitStack
import gzip
from pathlib import Path
import re
import sys
from typing import Dict, Iterable, Iterator, List, Literal, Optional, Sequence, TextIO, Tuple, Union, cast
from urllib.parse import quote, unquote


FilePath = Union[str, Path]
VcfInfo = Dict[str, Optional[str]]  # Bare INFO flags have no value.
MatchMode = Literal["overlap", "ins_site", "ins_length"]
Event = Tuple[str, MatchMode, int, int, int, int]
QueryInterval = Tuple[str, int, int, int, int, int, str]

FIELDS = ("gnomad", "gnomad_AF", "gnomad_FREQ_HOMALT", "gnomad_SVTYPE")
FIELD_TYPES = ("String", "Float", "Float", "String")
MATCHES = "GNOMAD_SV_MATCHES"
MODES: Tuple[MatchMode, ...] = ("overlap", "ins_site", "ins_length")
MAX_RECORDS = 50  # VEP's default custom-annotation limit, in reference source order.
SV_CLASSES = {
    "DEL": "deletion", "DUP": "duplication", "TDUP": "tandem_duplication",
    "INS": "insertion", "INV": "inversion", "CNV": "copy_number_variation",
    "TREP": "tandem_repeat", "BND": "chromosome_breakpoint",
}
for _kind, _name in (("DEL", "deletion"), ("INS", "insertion")):
    for _suffix, _prefix in (("ME", "mobile_element"), ("ALU", "Alu"),
                             ("HERV", "HERV"), ("LINE1", "LINE1"), ("SVA", "SVA")):
        SV_CLASSES[_kind + "_" + _suffix] = _prefix + "_" + _name


def open_vcf(path: FilePath) -> TextIO:
    """
    Detect if VCF is compressed or not and read it
    """
    return gzip.open(str(path), "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def records(path: FilePath) -> Iterator[Tuple[int, List[str], VcfInfo]]:
    """
    Open VCF file, iterate through, return:
    row: The raw row
    columns: The split row
    info: Info col parsed into key-value dict format
    """
    with open_vcf(path) as handle:
        row = 0
        for line in handle:
            if line.startswith("#"):
                continue
            row += 1
            columns = line.rstrip("\n").split("\t")
            info_fields = columns[7].split(";")
            info: VcfInfo = {}
            for entry in info_fields:
                key, separator, value = entry.partition("=")
                info[key] = value if separator else None
            yield row, columns, info


def schema(path: FilePath) -> Optional[List[str]]:
    """
    Detect the CSQ fields
    Return None if not present
    """
    with open_vcf(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ,"):
                return line.split("Format: ", 1)[1].split('"', 1)[0].split("|")
            if line.startswith("#CHROM"):
                break
    return None


def normalize_chr_name(name: str) -> str:
    """
    Normalize chromosome names
    chr3 -> 3
    M -> MT
    """
    name = name[3:] if name.startswith("chr") else name
    return "MT" if name == "M" else name


def symbolic_class(alt: str) -> Optional[str]:
    """Keep VEP's same-class matching, including tandem and mobile-element types."""
    mobile = re.search(r"(INS|DEL):ME:?([A-Z0-9]+)?", alt, re.I)
    if mobile:
        element = (mobile.group(2) or "ME").upper()
        element = "LINE1" if element == "L1" else element
        element = element if element in ("ALU", "HERV", "LINE1", "SVA") else "ME"
        return SV_CLASSES[mobile.group(1).upper() + "_" + element]
    if "DUP:TANDEM" in alt:
        return SV_CLASSES["TDUP"]
    if "CNV:TR" in alt:
        return SV_CLASSES["TREP"]
    if re.search(r"CN=?[0-9]", alt):
        return SV_CLASSES["DEL" if re.fullmatch(r"<?CN=?0>?", alt) else
                          "DUP" if re.fullmatch(r"<?CN=?2>?", alt) else "CNV"]
    if "CNV" in alt:
        return SV_CLASSES["CNV"]
    if any(char in alt for char in "[]") or alt.startswith(".") or alt.endswith("."):
        return SV_CLASSES["BND"]
    return SV_CLASSES.get(alt.strip("<>").split(":", 1)[0])


def metadata_end(pos: int, info: VcfInfo, kind: str) -> int:
    length = cast(str, info.get("SVLEN", ".")).split(",")[0]
    if length != ".":
        return pos + abs(int(length))
    default_end = pos if kind.endswith("deletion") else pos + 1
    return int(cast(Union[str, int], info.get("END", default_end)))


def event(columns: Sequence[str], info: VcfInfo) -> Optional[Event]:
    """Return the type, matching mode, BED span and VEP's original scan bounds."""
    pos, ref, alt = int(columns[1]), columns[3], columns[4]
    if "," in alt:
        return None
    kind = None
    if not re.fullmatch(r"[ACGT]+", ref + alt):
        # VEP falls back to SVTYPE for non-standard ALT spellings, including TDUP.
        variant_type = alt
        if not re.match(r"^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)", alt):
            variant_type = cast(str, info.get("SVTYPE", alt))
        kind = symbolic_class(variant_type)
    if kind:
        end = metadata_end(pos, info, kind)
        if kind == "chromosome_breakpoint" and info.get("SVLEN", ".") == ".":
            end = pos + 1  # END can identify the mate; VEP ignores it for query breakends.
        if kind.endswith("insertion"):
            # VEP swaps reversed insertion boundaries before checking overlap.
            return kind, "ins_length", min(pos + 1, end) - 1, max(pos + 1, end), pos, end
        return (kind, "overlap", pos, end, pos, end) if end > pos else None
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
        return "deletion", "overlap", start, start + len(ref), start, start + len(ref)
    if alt and not ref:
        # Preserve VEP's sequence-INS path: compare the insertion site, without
        # comparing inserted lengths. SVLEN does not affect this path.
        return "insertion", "ins_site", start - 1, start, start, start
    return None


def reference(vcf: FilePath, output: FilePath) -> None:
    directory = Path(output)
    directory.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        handles = {mode: stack.enter_context(gzip.open(str(directory / (mode + ".bed.gz")), "wt")) for mode in MODES}
        for row, columns, info in records(vcf):
            kind = symbolic_class(columns[4])
            if kind is None:
                continue
            pos = int(columns[1])
            end = metadata_end(pos, info, kind)
            # Original tabix bounds matter for legacy SV matching: VEP cannot
            # match a reference record that its original VCF index did not fetch.
            reference_end = info.get("END")
            index_end = int(reference_end) if reference_end else pos - 1 + len(columns[3])
            spans: Dict[MatchMode, Tuple[int, int]]
            if not kind.endswith("insertion"):
                if end <= pos:
                    raise ValueError("Incomplete reference SV: " + columns[2])
                spans = {"overlap": (pos, end)}
            else:
                spans = {"ins_site": (pos - 1, pos - 1 + len(columns[3]))}
                if end > pos:
                    spans["ins_length"] = (pos, end)
            name = columns[2].split(";", 1)[0]
            if name in ("", ".", "0"):
                name = "{}:{}-{}".format(columns[0], pos + 1, end)
            values: List[str] = []
            for key in ("AF", "FREQ_HOMALT", "SVTYPE"):
                value = info.get(key, "")
                value = "1" if value is None else value
                parts = value.split(",")
                values.append(parts[1] if len(parts) == 2 else value)  # VEP's REF/ALT-value handling.
            # One payload per match keeps missing fields aligned while bcftools
            # appends hits. Source row and coordinates preserve VEP order and
            # allow an exact integer check at the 70% boundary.
            for mode, span in spans.items():
                payload = "|".join(quote(str(value), safe="") for value in (row, *span, index_end, kind, name, *values))
                handles[mode].write("{}\t{}\t{}\t{}\n".format(normalize_chr_name(columns[0]), *span, payload))


def query(vcf: FilePath, output: FilePath, synonyms: Optional[FilePath] = None) -> None:
    fields = schema(vcf)
    class_index = fields.index("VARIANT_CLASS") if fields else None
    aliases: Dict[str, str] = {}
    if synonyms:
        with open(synonyms) as handle:
            for line in handle:
                names = [normalize_chr_name(name) for name in line.split()] if not line.startswith("#") else []
                if names:
                    aliases.update((name, names[0]) for name in names)
    intervals: Dict[MatchMode, List[QueryInterval]] = {mode: [] for mode in MODES}
    unclassified: Dict[str, int] = Counter()
    for row, columns, info in records(vcf):
        classes = None
        if class_index is not None:
            if "CSQ" not in info:
                continue  # Keep VEP's skipped records unannotated.
            csq = cast(str, info["CSQ"])
            classes = {entry.split("|")[class_index] for entry in csq.split(",")}
            selected = classes.intersection(SV_CLASSES.values())
            if not selected:
                continue
            if len(classes) != 1:
                raise ValueError("Multiple VEP classes for " + columns[2])
        parsed = event(columns, info)
        if parsed is None:
            if fields:
                raise ValueError("Cannot represent SV for matching: {}:{} ({})".format(*columns[:3]))
            svtype = info.get("SVTYPE")
            if svtype == "DEL" or svtype == "INS":
                unclassified[svtype] += 1
            continue
        kind, mode, start, end, seek_start, native_end = parsed
        if classes is not None and classes != {kind}:
            raise ValueError("VEP class disagrees with alleles for " + columns[2])
        chrom = normalize_chr_name(columns[0])
        intervals[mode].append((aliases.get(chrom, chrom), start, end, row, seek_start, native_end, kind))
    for mode, rows in intervals.items():
        with open("{}.{}.intervals.vcf".format(output, mode), "w") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            for chrom in sorted({item[0] for item in rows}):
                handle.write("##contig=<ID={}>\n".format(chrom))
            for key, datatype in (("END", "Integer"), ("SEEK_START", "Integer"), ("QUERY_END", "Integer"),
                                  ("MATCH_MODE", "String"), ("SVCLASS", "String")):
                handle.write('##INFO=<ID={},Number=1,Type={},Description="Internal matching coordinate or mode">\n'.format(key, datatype))
            handle.write('##INFO=<ID={},Number=.,Type=String,Description="Internal gnomAD matches">\n'.format(MATCHES))
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for chrom, start, end, row, seek_start, native_end, kind in sorted(rows):
                handle.write("{}\t{}\t{}\tN\t<SV>\t.\tPASS\tEND={};MATCH_MODE={};SEEK_START={};QUERY_END={};SVCLASS={}\n".format(
                    chrom, start + 1, row, end, mode, seek_start, native_end, kind))
    if unclassified:
        print("Left complex records unchanged ({}); VEP classification is needed.".format(dict(unclassified)), file=sys.stderr)


def matched_values(matches_vcf: FilePath) -> Dict[int, List[str]]:
    annotations: Dict[int, List[str]] = {}
    for _, columns, info in records(matches_vcf):
        start, end = int(columns[1]) - 1, int(cast(str, info["END"]))
        hits: List[Tuple[int, List[str]]] = []
        for entry in (info.get(MATCHES) or "").split(","):
            if not entry:
                continue
            row, rstart, rend, index_end, kind, *values = [unquote(value) for value in entry.split("|")]
            if kind != info["SVCLASS"]:
                continue
            rstart, rend = int(rstart), int(rend)
            overlap = min(end, rend) - max(start, rstart)
            mode = info["MATCH_MODE"]
            if mode != "ins_site" and (10 * overlap < 7 * (end - start) or 10 * overlap < 7 * (rend - rstart)):
                continue
            if kind != "deletion":
                reference_pos = rstart + 1 if mode == "ins_site" else rstart
                if (
                    int(index_end) <= int(cast(str, info["SEEK_START"]))
                    or reference_pos > int(cast(str, info["QUERY_END"]))
                ):
                    continue
            hits.append((int(row), values))
        selected = [values for _, values in sorted(hits)[:MAX_RECORDS]]
        values = [",".join(quote(hit[i], safe=":._-|") if i in (0, 3) and hit[i] else hit[i] or "."
                           for hit in selected) for i in range(4)]
        if len(hits) > MAX_RECORDS:
            values[0] += ",..."
        annotations[int(columns[2])] = values
    return annotations


def restore(vcf: FilePath, matches_vcfs: Iterable[FilePath], output: FilePath) -> None:
    annotations: Dict[int, List[str]] = {}
    for path in matches_vcfs:
        annotations.update(matched_values(path))
    row = 0
    with open_vcf(vcf) as source, open(output, "w") as target:
        for line in source:
            if line.startswith("#"):
                if line.startswith(tuple("##INFO=<ID=" + field + "," for field in FIELDS)):
                    continue
                if line.startswith("#CHROM"):
                    for field, datatype in zip(FIELDS, FIELD_TYPES):
                        target.write('##INFO=<ID={},Number=.,Type={},Description="gnomAD SV annotation">\n'.format(field, datatype))
                target.write(line)
                continue
            row += 1
            if row not in annotations:
                target.write(line)
                continue
            columns = line.rstrip("\n").split("\t")
            info = columns[7].split(";")
            info = [entry for entry in info if entry.split("=", 1)[0] not in FIELDS and entry != "."]
            info.extend("{}={}".format(field, value) for field, value in zip(FIELDS, annotations[row]) if value)
            columns[7] = ";".join(info) or "."
            target.write("\t".join(columns) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode", choices=("reference", "query", "restore"),
        help="Build reference BEDs, prepare query intervals, or restore annotations to the original VCF.",
    )
    parser.add_argument("--vcf", required=True, help="Input VCF, plain or gzip-compressed.")
    parser.add_argument("--output", required=True, help="Reference directory, query prefix, or annotated VCF path.")
    parser.add_argument("--synonyms", help="Chromosome alias file for query mode.")
    parser.add_argument("--matches", nargs="+", help="VCFs containing bcftools matches. Required in restore mode.")
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
