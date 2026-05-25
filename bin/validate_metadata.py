#!/usr/bin/env python3

import argparse
import csv
import json
import shutil
from datetime import date
from pathlib import Path


REQUIRED_COLUMNS = ("sample_name", "collection_date")


def parse_args():
    parser = argparse.ArgumentParser(description="Validate sample metadata for MK Flu-Pipe")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--validated-samplesheet", required=True)
    return parser.parse_args()


def read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise SystemExit("Metadata CSV is empty or has no header")
        fieldnames = [name.strip() for name in reader.fieldnames]
        if len(fieldnames) != len(set(fieldnames)):
            raise SystemExit("Metadata CSV has duplicated column names")
        reader.fieldnames = fieldnames
        rows = [{key: str(value or "").strip() for key, value in row.items()} for row in reader]
    return fieldnames, rows


def read_sample_ids(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [row["sample_id"].strip() for row in csv.DictReader(handle)]


def validate_rows(fieldnames, rows, sample_ids):
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
    if missing_columns:
        raise SystemExit(
            "Metadata CSV is missing required column(s): " + ", ".join(missing_columns)
        )
    if not rows:
        raise SystemExit("Metadata CSV contains no sample rows")

    errors = []
    rows_by_sample = {}
    for line_number, row in enumerate(rows, start=2):
        sample_name = row["sample_name"]
        collection_date = row["collection_date"]
        if not sample_name:
            errors.append(f"line {line_number}: sample_name is empty")
        elif sample_name in rows_by_sample:
            errors.append(f"line {line_number}: duplicate sample_name '{sample_name}'")
        else:
            rows_by_sample[sample_name] = row
        try:
            date.fromisoformat(collection_date)
        except ValueError:
            errors.append(
                f"line {line_number}: collection_date '{collection_date}' must use YYYY-MM-DD"
            )

    expected = set(sample_ids)
    provided = set(rows_by_sample)
    missing_samples = sorted(expected - provided)
    extra_samples = sorted(provided - expected)
    if missing_samples:
        errors.append("samples missing from metadata: " + ", ".join(missing_samples))
    if extra_samples:
        errors.append("metadata samples absent from FASTQ input: " + ", ".join(extra_samples))
    if errors:
        raise SystemExit("Invalid metadata CSV: " + "; ".join(errors))

    return [rows_by_sample[sample_id] for sample_id in sample_ids]


def main():
    args = parse_args()
    metadata_path = Path(args.metadata)
    samplesheet_path = Path(args.samplesheet)
    output_path = Path(args.output)
    summary_path = Path(args.summary)
    validated_samplesheet_path = Path(args.validated_samplesheet)

    fieldnames, rows = read_csv(metadata_path)
    sample_ids = read_sample_ids(samplesheet_path)
    ordered_rows = validate_rows(fieldnames, rows, sample_ids)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ordered_rows)

    summary_path.write_text(
        json.dumps(
            {
                "metadata_file": str(metadata_path),
                "sample_count": len(ordered_rows),
                "columns": fieldnames,
                "required_columns": list(REQUIRED_COLUMNS),
                "optional_location_columns_present": [
                    name for name in ("country", "state", "city", "location") if name in fieldnames
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    shutil.copyfile(samplesheet_path, validated_samplesheet_path)


if __name__ == "__main__":
    main()
