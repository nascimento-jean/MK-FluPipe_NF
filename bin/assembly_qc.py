#!/usr/bin/env python3

import argparse
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
                    yield header[1:].split()[0], "".join(seq_chunks)
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header[1:].split()[0], "".join(seq_chunks)


def read_segment_coverages(tables_dir: Path):
    coverage = {}
    if not tables_dir.exists():
        return coverage

    for coverage_file in tables_dir.glob("*-coverage.txt"):
        values = []
        with coverage_file.open("r", encoding="utf-8", errors="ignore") as handle:
            next(handle, None)
            for raw_line in handle:
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    values.append(float(fields[2]))
                except ValueError:
                    continue
        average = round(sum(values) / len(values)) if values else 0
        coverage[coverage_file.name.replace("-coverage.txt", "")] = int(average)
    return coverage


def main():
    parser = argparse.ArgumentParser(description="Compute post-assembly QC for an IRMA run")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--consensus-fasta", required=True)
    parser.add_argument("--irma-dir", required=True)
    parser.add_argument("--min-coverage", type=int, required=True)
    parser.add_argument("--max-n-pct", type=float, required=True)
    parser.add_argument("--min-segments", type=int, required=True)
    parser.add_argument("--output-tsv", required=True)
    args = parser.parse_args()

    consensus_fasta = Path(args.consensus_fasta)
    irma_dir = Path(args.irma_dir)
    output_tsv = Path(args.output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    segment_count = 0
    high_n_segments = []
    for segment_name, sequence in parse_fasta(consensus_fasta):
        segment_count += 1
        sequence_length = len(sequence)
        if sequence_length == 0:
            continue
        n_count = sum(1 for base in sequence if base in {"N", "n"})
        n_pct = (n_count / sequence_length) * 100
        if n_pct > args.max_n_pct:
            high_n_segments.append(f"{segment_name}({n_pct:.1f}%N)")

    low_cov_segments = []
    for segment_name, avg_cov in read_segment_coverages(irma_dir / "tables").items():
        if avg_cov < args.min_coverage:
            low_cov_segments.append(f"{segment_name}({avg_cov}x)")

    qc_status = "PASS"
    if segment_count < args.min_segments:
        qc_status = "FAIL"
    elif high_n_segments or low_cov_segments:
        qc_status = "WARN"

    detail = f"segs:{segment_count}/8"
    if low_cov_segments:
        detail += f"|low_cov:{' '.join(low_cov_segments)}"
    if high_n_segments:
        detail += f"|high_N:{' '.join(high_n_segments)}"

    with output_tsv.open("w", encoding="utf-8") as handle:
        handle.write("sample\tqc_assembly\tqc_detail\n")
        handle.write(f"{args.sample_id}\t{qc_status}\t{detail}\n")


if __name__ == "__main__":
    main()
