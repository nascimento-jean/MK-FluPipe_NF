#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def parse_args():
    parser = argparse.ArgumentParser(description="Discover MK Flu-Pipe samples")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--seq-type", default="auto")
    parser.add_argument("--seq-mode", default="")
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


def classify_illumina_paired(path: Path):
    name = path.name

    for suffix in FASTQ_SUFFIXES:
        token = f"_L001_R1_001{suffix}"
        if name.endswith(token):
            sample_id = name[: -len(token)]
            mate_name = f"{sample_id}_L001_R2_001{suffix}"
            return sample_id, mate_name
    return None


def classify_sra_paired(path: Path):
    name = path.name

    for suffix in FASTQ_SUFFIXES:
        token = f"_1{suffix}"
        if name.endswith(token):
            sample_id = name[: -len(token)]
            mate_name = f"{sample_id}_2{suffix}"
            return sample_id, mate_name
    return None


def classify_generic_paired(path: Path):
    name = path.name

    for suffix in FASTQ_SUFFIXES:
        token = f"_R1{suffix}"
        if name.endswith(token):
            sample_id = name[: -len(token)]
            mate_name = f"{sample_id}_R2{suffix}"
            return sample_id, mate_name

        mid = "_R1_"
        if mid in name and name.endswith(suffix):
            prefix, rest = name.rsplit(mid, 1)
            sample_id = prefix
            mate_name = f"{prefix}_R2_{rest}"
            return sample_id, mate_name

    return None


def classify_short_fastq(path: Path, seq_mode: str):
    if seq_mode == "illumina_paired":
        return classify_illumina_paired(path)
    if seq_mode == "sra_paired":
        return classify_sra_paired(path)
    if seq_mode == "generic_paired":
        return classify_generic_paired(path)

    for classifier in (classify_illumina_paired, classify_sra_paired, classify_generic_paired):
        match = classifier(path)
        if match:
            return match
    return None


def discover(input_dir: Path, seq_type: str, seq_mode: str):
    files = sorted([p for p in input_dir.iterdir() if p.is_file() and is_fastq(p)])
    consumed = set()
    rows = []

    for path in files:
        if path in consumed:
            continue

        stem = strip_fastq_suffix(path.name)

        forced_single = seq_mode == "single" or seq_type in {"short_single", "long"}
        if not forced_single and seq_type in {"auto", "short_paired"}:
            match = classify_short_fastq(path, seq_mode)
            if match:
                sample_id, mate_name = match
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


def write_summary(path: Path, input_dir: Path, seq_type: str, seq_mode: str, rows):
    summary = {
        "input_dir": str(input_dir.resolve()),
        "requested_seq_type": seq_type,
        "requested_seq_mode": seq_mode,
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

    rows = discover(input_dir, args.seq_type, args.seq_mode)
    if not rows:
        raise SystemExit(f"No FASTQ files found in: {input_dir}")

    write_samplesheet(Path(args.samplesheet), rows)
    write_summary(Path(args.summary), input_dir, args.seq_type, args.seq_mode, rows)


if __name__ == "__main__":
    main()
