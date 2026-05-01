#!/usr/bin/env python3

import argparse
import gzip
import statistics
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Generate aggregate FASTQ set statistics")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--layout", required=True)
    parser.add_argument("--seq-type", required=True)
    parser.add_argument("--input-fastq", nargs="+", required=True)
    parser.add_argument("--output-fastq", nargs="+", required=True)
    parser.add_argument("--output-tsv", required=True)
    return parser.parse_args()


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def n50(lengths):
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length
    return 0


def collect_lengths(paths):
    lengths = []
    for path in paths:
        with open_maybe_gzip(path) as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline().rstrip("\n")
                plus = handle.readline()
                qual = handle.readline()
                if not plus or not qual:
                    raise SystemExit(f"Malformed FASTQ: {path}")
                lengths.append(len(seq))
    return lengths


def stats_for(paths):
    lengths = collect_lengths(paths)
    read_count = len(lengths)
    total_bases = sum(lengths)
    mean_len = round(total_bases / read_count, 2) if read_count else 0
    median_len = round(statistics.median(lengths), 2) if lengths else 0
    min_len = min(lengths) if lengths else 0
    max_len = max(lengths) if lengths else 0
    return {
        "reads": read_count,
        "bases": total_bases,
        "mean_len": mean_len,
        "median_len": median_len,
        "min_len": min_len,
        "max_len": max_len,
        "n50": n50(lengths),
    }


def retention(before, after):
    if before == 0:
        return 0
    return round((after / before) * 100, 2)


def main():
    args = parse_args()
    input_paths = [Path(p) for p in args.input_fastq]
    output_paths = [Path(p) for p in args.output_fastq]
    output_tsv = Path(args.output_tsv)

    before = stats_for(input_paths)
    after = stats_for(output_paths)

    fields = [
        "sample_id",
        "layout",
        "seq_type",
        "input_files",
        "output_files",
        "input_reads",
        "output_reads",
        "removed_reads",
        "read_retention_pct",
        "input_bases",
        "output_bases",
        "removed_bases",
        "base_retention_pct",
        "input_mean_len",
        "output_mean_len",
        "input_median_len",
        "output_median_len",
        "input_n50",
        "output_n50",
        "input_min_len",
        "output_min_len",
        "input_max_len",
        "output_max_len",
    ]

    values = {
        "sample_id": args.sample_id,
        "layout": args.layout,
        "seq_type": args.seq_type,
        "input_files": ",".join(str(p) for p in input_paths),
        "output_files": ",".join(str(p) for p in output_paths),
        "input_reads": before["reads"],
        "output_reads": after["reads"],
        "removed_reads": before["reads"] - after["reads"],
        "read_retention_pct": retention(before["reads"], after["reads"]),
        "input_bases": before["bases"],
        "output_bases": after["bases"],
        "removed_bases": before["bases"] - after["bases"],
        "base_retention_pct": retention(before["bases"], after["bases"]),
        "input_mean_len": before["mean_len"],
        "output_mean_len": after["mean_len"],
        "input_median_len": before["median_len"],
        "output_median_len": after["median_len"],
        "input_n50": before["n50"],
        "output_n50": after["n50"],
        "input_min_len": before["min_len"],
        "output_min_len": after["min_len"],
        "input_max_len": before["max_len"],
        "output_max_len": after["max_len"],
    }

    output_tsv.write_text(
        "\t".join(fields) + "\n" + "\t".join(str(values[field]) for field in fields) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
