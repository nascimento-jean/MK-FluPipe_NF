#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


HEADER = ["sample", "coinfection_status", "details"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze IRMA allAlleles outputs for potential coinfection/mixed infection"
    )
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--minority-freq", type=float, required=True)
    parser.add_argument("--coinfection-pct", type=float, required=True)
    parser.add_argument("irma_dirs", nargs="*")
    return parser.parse_args()


def analyze_allele_file(path: Path, minority_freq: float):
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return 0, 0, 0.0

        col_index = {name: idx for idx, name in enumerate(header)}
        required = ["Position", "Frequency", "Total", "Allele_Type"]
        if any(name not in col_index for name in required):
            return 0, 0, 0.0

        position_totals = {}
        minority_sum = {}

        for row in reader:
            try:
                pos = row[col_index["Position"]].strip()
                total = float(row[col_index["Total"]])
                freq = float(row[col_index["Frequency"]])
                allele_type = row[col_index["Allele_Type"]].strip()
            except (IndexError, ValueError):
                continue

            if not pos or total < 10:
                continue

            position_totals[pos] = total
            if allele_type == "Minority":
                minority_sum[pos] = minority_sum.get(pos, 0.0) + freq

        total_positions = len(position_totals)
        flagged_positions = sum(
            1 for pos in position_totals if minority_sum.get(pos, 0.0) >= minority_freq
        )
        flagged_pct = (flagged_positions / total_positions) * 100.0 if total_positions else 0.0
        return flagged_positions, total_positions, flagged_pct


def analyze_sample(sample: str, irma_dir: Path, minority_freq: float, coinfection_pct: float):
    tables_dir = irma_dir / "tables"
    allele_files = sorted(tables_dir.glob("*-allAlleles.txt"))

    details = []
    total_segments = 0
    flagged_segments = 0

    for allele_file in allele_files:
        total_segments += 1
        segment_name = allele_file.name[: -len("-allAlleles.txt")]
        flagged_pos, total_pos, flagged_pct = analyze_allele_file(allele_file, minority_freq)
        if flagged_pct >= coinfection_pct:
            flagged_segments += 1
            details.append(f"{segment_name}({flagged_pos}/{total_pos}pos,{flagged_pct:.2f}%)")

    status = "WARN" if flagged_segments else "OK"
    suffix = " ".join(details) if details else "no_alert"
    detail_text = f"segs_analyzed:{total_segments}|segs_flagged:{flagged_segments}|{suffix}"
    return status, detail_text, total_segments, flagged_segments


def main():
    args = parse_args()
    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    log_lines = []
    rows_written = 0

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(HEADER)

        for irma_dir_arg in args.irma_dirs:
            irma_dir = Path(irma_dir_arg)
            sample = irma_dir.name
            if not irma_dir.exists():
                log_lines.append(f"Skipping {sample}: IRMA directory not found")
                continue

            status, details, total_segments, flagged_segments = analyze_sample(
                sample, irma_dir, args.minority_freq, args.coinfection_pct
            )
            writer.writerow([sample, status, details])
            rows_written += 1
            log_lines.append(
                f"{sample}: {status} | analyzed={total_segments} | flagged={flagged_segments} | {details}"
            )

    if rows_written == 0:
        log_lines.append("No IRMA directories were available for coinfection analysis.")

    log_lines.append(f"Total samples written: {rows_written}")
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
