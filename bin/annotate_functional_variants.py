#!/usr/bin/env python3
"""Create a structured functional annotation table from full variant calls."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SEGMENT_NAMES = {
    "1": "PB2",
    "2": "PB1",
    "3": "PA",
    "4": "HA",
    "5": "NP",
    "6": "NA",
    "7": "MP",
    "8": "NS",
}
SEGMENT_NUMBERS = {name: number for number, name in SEGMENT_NAMES.items()}
SEGMENT_NUMBERS.update({"M": "7", "M1": "7", "M2": "7", "NEP": "8", "NS1": "8", "NS2": "8"})

FIELDS = [
    "sample",
    "type",
    "subtype",
    "segment",
    "segment_name",
    "gene",
    "aa_position",
    "ref_aa",
    "alt_aa",
    "mutation",
    "effect",
    "impact",
    "frequency",
    "depth",
    "annotation_source",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def typing_lookup(path: Path) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        sample = row.get("sample") or row.get("Sample") or ""
        if sample:
            lookup[sample] = row
    return lookup


def parse_segment(segment: str) -> tuple[str, str, str]:
    segment = segment or ""
    if segment in SEGMENT_NAMES:
        return "", "", segment
    if segment.upper() in SEGMENT_NUMBERS:
        return "", "", SEGMENT_NUMBERS[segment.upper()]

    parts = segment.rsplit("_", 2)
    if len(parts) == 3:
        sample, flu_type, segment_num = parts
        if segment_num.upper() in SEGMENT_NUMBERS:
            segment_num = SEGMENT_NUMBERS[segment_num.upper()]
        return sample, flu_type, segment_num
    return "", "", segment or "-"


def subtype_for(sample: str, flu_type: str, typing: dict[str, dict[str, str]]) -> str:
    row = typing.get(sample, {})
    if flu_type == "B":
        return "B"
    ha = row.get("subtype_HA") or row.get("Subtype Ha") or ""
    na = row.get("subtype_NA") or row.get("Subtype Na") or ""
    subtype = f"{ha}{na}".replace("-", "").replace("nd", "")
    return subtype or "-"


def classify_effect(gene: str, ref_aa: str, alt_aa: str, mutation: str) -> tuple[str, str]:
    gene = (gene or "").strip()
    ref_aa = (ref_aa or "").strip()
    alt_aa = (alt_aa or "").strip()
    mutation = (mutation or "").strip()

    if gene.upper() == "UTR":
        return "non_coding", "MODIFIER"
    if not mutation or mutation == ".":
        return "no_amino_acid_change", "LOW"
    if alt_aa == "*":
        return "stop_gained", "HIGH"
    if ref_aa == "*":
        return "stop_lost", "HIGH"
    if ref_aa == alt_aa:
        return "synonymous", "LOW"
    return "amino_acid_change", "MODERATE"


def annotate_rows(fullvar_rows: list[dict[str, str]], typing: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    annotated: list[dict[str, str]] = []
    for row in fullvar_rows:
        sample = row.get("sample", "")
        flu_type = row.get("type", "")
        _, parsed_type, segment_num = parse_segment(row.get("segment", ""))
        if not flu_type:
            flu_type = parsed_type or "-"
        effect, impact = classify_effect(row.get("gene", ""), row.get("ref_aa", ""), row.get("alt_aa", ""), row.get("mutation", ""))
        annotated.append(
            {
                "sample": sample,
                "type": flu_type or "-",
                "subtype": subtype_for(sample, flu_type, typing),
                "segment": segment_num,
                "segment_name": SEGMENT_NAMES.get(segment_num, "-"),
                "gene": row.get("gene", "-") or "-",
                "aa_position": row.get("aa_position", "-") or "-",
                "ref_aa": row.get("ref_aa", "-") or "-",
                "alt_aa": row.get("alt_aa", "-") or "-",
                "mutation": row.get("mutation", "-") or "-",
                "effect": effect,
                "impact": impact,
                "frequency": row.get("frequency", "-") or "-",
                "depth": row.get("depth", "-") or "-",
                "annotation_source": "RefSeq GFF3 + iVar protein table",
            }
        )
    return annotated


def main() -> int:
    parser = argparse.ArgumentParser(description="Build MK-FluPipe functional annotation table")
    parser.add_argument("--fullvar-tsv", required=True)
    parser.add_argument("--typing-tsv", required=True)
    parser.add_argument("--output-tsv", required=True)
    args = parser.parse_args()

    fullvar_rows = read_tsv(Path(args.fullvar_tsv))
    typing = typing_lookup(Path(args.typing_tsv))
    write_tsv(Path(args.output_tsv), annotate_rows(fullvar_rows, typing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
