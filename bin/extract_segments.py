#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def parse_fasta(path: Path):
    header = None
    seq_chunks = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


def normalize_sequence(sequence: str) -> str:
    return re.sub(r"[^ACGTacgt]", "N", sequence)


def detect_segment(header: str):
    token = header[1:].split()[0]
    match = re.search(r"_(\d)$", token)
    if not match:
        return None
    value = int(match.group(1))
    if 1 <= value <= 8:
        return value
    return None


def main():
    parser = argparse.ArgumentParser(description="Extract influenza segments 1-8 from a consensus FASTA")
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    input_fasta = Path(args.input_fasta)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    segment_entries = {segment: [] for segment in range(1, 9)}

    for header, sequence in parse_fasta(input_fasta):
        segment = detect_segment(header)
        if segment is None:
            continue
        segment_entries[segment].append((header, normalize_sequence(sequence)))

    for segment in range(1, 9):
        segment_dir = output_dir / f"single_segment{segment}"
        segment_dir.mkdir(parents=True, exist_ok=True)
        sample_file = segment_dir / f"{args.sample_id}_segment_{segment}.fasta"
        entries = segment_entries[segment]
        with sample_file.open("w", encoding="utf-8") as handle:
            for header, sequence in entries:
                handle.write(f"{header}\n{sequence}\n")


if __name__ == "__main__":
    main()
