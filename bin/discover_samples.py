#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
PAIRED_PATTERNS = [
    ("_R1", "_R2"),
    (".R1", ".R2"),
    ("_1", "_2"),
    (".1", ".2"),
]


def parse_args():
    parser = argparse.ArgumentParser(description="Discover MK Flu-Pipe samples")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--seq-type", default="auto")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--summary", required=True)
    return parser.parse_args()


def is_fastq(path: Path) -> bool:
    return any(str(path).endswith(suffix) for suffix in FASTQ_SUFFIXES)


def strip_fastq_suffix(name: str) -> str:
    for suffix in FASTQ_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def paired_match(stem: str):
    for read1, read2 in PAIRED_PATTERNS:
        if read1 in stem:
            return stem.replace(read1, ""), read1, read2
    return None


def discover(input_dir: Path, seq_type: str):
    files = sorted([p for p in input_dir.iterdir() if p.is_file() and is_fastq(p)])
    consumed = set()
    rows = []

    for path in files:
        if path in consumed:
            continue

        stem = strip_fastq_suffix(path.name)

        if seq_type in {"auto", "short_paired"}:
            match = paired_match(stem)
            if match:
                sample_id, read1_token, read2_token = match
                mate_name = path.name.replace(read1_token, read2_token, 1)
                mate_path = input_dir / mate_name
                if mate_path.exists():
                    consumed.add(path)
                    consumed.add(mate_path)
                    rows.append(
                        {
                            "sample_id": sample_id,
                            "r1": str(path.resolve()),
                            "r2": str(mate_path.resolve()),
                            "layout": "paired",
                            "seq_type": "short_paired",
                        }
                    )
                    continue

        layout = "single"
        inferred = seq_type if seq_type != "auto" else "short_single"
        if seq_type == "long":
            inferred = "long"

        consumed.add(path)
        rows.append(
            {
                "sample_id": stem,
                "r1": str(path.resolve()),
                "r2": "",
                "layout": layout,
                "seq_type": inferred,
            }
        )

    rows.sort(key=lambda item: item["sample_id"])
    return rows


def write_samplesheet(path: Path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "r1", "r2", "layout", "seq_type"],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, input_dir: Path, seq_type: str, rows):
    summary = {
        "input_dir": str(input_dir.resolve()),
        "requested_seq_type": seq_type,
        "sample_count": len(rows),
        "paired_count": sum(1 for row in rows if row["layout"] == "paired"),
        "single_count": sum(1 for row in rows if row["layout"] == "single"),
        "samples": rows,
    }
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        raise SystemExit(f"Input directory not found: {input_dir}")
    if not input_dir.is_dir():
        raise SystemExit(f"Input path is not a directory: {input_dir}")

    rows = discover(input_dir, args.seq_type)
    if not rows:
        raise SystemExit(f"No FASTQ files found in: {input_dir}")

    write_samplesheet(Path(args.samplesheet), rows)
    write_summary(Path(args.summary), input_dir, args.seq_type, rows)


if __name__ == "__main__":
    main()

