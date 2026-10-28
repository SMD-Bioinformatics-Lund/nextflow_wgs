"""Join original AnnotSV TSV rows without confusing net length with REF span.

Full AnnotSV rows retain the original ID/REF/ALT and INFO. Those identify the
input allele even when AnnotSV derives a different SV_end for a literal indel.
Symbolic alleles still require the original END. Legacy coordinate-only tables
are supported by exact, unambiguous coordinate matching only.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re
import sys


MISSING = (None, "", ".", "NA")


def _chrom(value):
    return value[3:] if value.startswith("chr") else value


def _kind(value):
    return "DUP" if value in ("TDUP", "DUP:TANDEM") else value


@dataclass
class _Row:
    key: tuple
    annotations: dict


class AnnotSVTable:
    def __init__(self, fields):
        self.fields = tuple(fields)
        self.by_id = {}
        self.by_allele = {}
        self.legacy_coordinates = {}

    @classmethod
    def read(cls, path, fields=("ACMG_class",)):
        table = cls(fields)
        previous = csv.field_size_limit()
        csv.field_size_limit(sys.maxsize)
        try:
            with Path(path).open(newline="", encoding="utf-8-sig") as stream:
                reader = csv.DictReader(stream, delimiter="\t")
                # Explicit empty outputs from the existing workflow and stubs.
                if reader.fieldnames in (["1"], ["AnnotSV_ID", "ACMG_class"]):
                    if next(reader, None) is not None:
                        raise ValueError("Unexpected rows in an empty AnnotSV table")
                    return table
                required = {"SV_chrom", "SV_start", "SV_end", "SV_type", *fields}
                if not required.issubset(reader.fieldnames or ()):
                    raise ValueError("AnnotSV TSV lacks required columns: " +
                                     ", ".join(sorted(required - set(reader.fieldnames or ()))))
                for row in reader:
                    if None in row:
                        raise ValueError("AnnotSV TSV row has more columns than its header")
                    if row.get("Annotation_mode", "full") != "full":
                        continue
                    table._add(row)
        finally:
            csv.field_size_limit(previous)
        return table

    def _add(self, row):
        chrom, pos, end, kind = (_chrom(row["SV_chrom"]), int(row["SV_start"]),
                                  int(row["SV_end"]), _kind(row["SV_type"]))
        annotations = {key: row[key] for key in self.fields if row.get(key) not in MISSING}
        if "ACMG_class" in annotations:
            if annotations["ACMG_class"] not in ("1", "2", "3", "4", "5"):
                raise ValueError("Invalid AnnotSV ACMG_class")
        ref, alt = row.get("REF"), row.get("ALT")
        if ref in MISSING or alt in MISSING:
            key = (chrom, pos, end, kind)
            self.legacy_coordinates.setdefault(key, []).append(_Row(key, annotations))
            return
        explicit_end = re.search(r"(?:^|;)END=([0-9]+)(?:;|$)", row.get("INFO") or "")
        if explicit_end:
            original_end = int(explicit_end[1])
        elif re.fullmatch("[ACGTNacgtn]+", ref) and re.fullmatch("[ACGTNacgtn]+", alt):
            original_end = pos + len(ref) - 1
        else:
            original_end = end
        key = (chrom, pos, original_end, kind, ref, alt)
        entry = _Row(key, annotations)
        identifier = row.get("ID")
        if identifier not in MISSING:
            previous = self.by_id.get(identifier)
            if previous is not None and (previous.key != key or previous.annotations != annotations):
                raise ValueError("Conflicting full AnnotSV rows for ID " + identifier)
            self.by_id[identifier] = entry
        self.by_allele.setdefault(key, []).append(entry)

    @staticmethod
    def _concordant(rows):
        annotations = rows[0].annotations
        if any(row.annotations != annotations for row in rows[1:]):
            raise ValueError("Ambiguous AnnotSV annotations for the same input allele/coordinates")
        return dict(annotations)

    def lookup(self, record):
        if len(record.alts or ()) != 1:
            raise ValueError("AnnotSV joins require biallelic input records")
        chrom, pos, end, kind = (_chrom(record.contig), record.pos, record.stop,
                                 _kind(record.info.get("SVTYPE")))
        key = (chrom, pos, end, kind, record.ref, record.alts[0])
        entry = self.by_id.get(record.id)
        if entry is not None:
            if entry.key != key:
                raise ValueError("AnnotSV ID describes a different allele or coordinates: " + record.id)
            return dict(entry.annotations)
        candidates = self.by_allele.get(key)
        if candidates:
            return self._concordant(candidates)
        candidates = self.legacy_coordinates.get((chrom, pos, end, kind))
        return self._concordant(candidates) if candidates else {}
