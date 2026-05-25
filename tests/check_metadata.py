#!/usr/bin/env python3
"""Smoke tests for metadata validation without running analysis tools."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATE = REPO_ROOT / "bin" / "validate_metadata.py"


def run_validation(metadata: Path, expect_success: bool) -> subprocess.CompletedProcess[str]:
    with tempfile.TemporaryDirectory(prefix="mkflupipe_metadata_") as tmp:
        tmpdir = Path(tmp)
        samplesheet = tmpdir / "samplesheet.csv"
        samplesheet.write_text(
            "sample_id,r1,r2,layout,seq_type\n"
            "SAMPLE_A,/tmp/a_R1.fastq.gz,/tmp/a_R2.fastq.gz,paired,short_paired\n",
            encoding="utf-8",
        )
        command = [
            sys.executable,
            str(VALIDATE),
            "--metadata",
            str(metadata),
            "--samplesheet",
            str(samplesheet),
            "--output",
            str(tmpdir / "validated_metadata.csv"),
            "--summary",
            str(tmpdir / "metadata_validation.json"),
            "--validated-samplesheet",
            str(tmpdir / "validated_samplesheet.csv"),
        ]
        result = subprocess.run(command, cwd=REPO_ROOT, capture_output=True, text=True)
        if expect_success:
            assert result.returncode == 0, result.stderr
            with (tmpdir / "validated_metadata.csv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            summary = json.loads((tmpdir / "metadata_validation.json").read_text(encoding="utf-8"))
            assert rows[0]["sample_name"] == "SAMPLE_A", rows
            assert rows[0]["collection_date"] == "2026-05-24", rows
            assert summary["sample_count"] == 1, summary
            assert summary["optional_location_columns_present"] == ["country", "state", "city"], summary
        else:
            assert result.returncode != 0, result.stdout
        return result


def main() -> int:
    metadata_dir = REPO_ROOT / "tests" / "data" / "metadata"
    run_validation(metadata_dir / "valid_short.csv", expect_success=True)
    missing = run_validation(metadata_dir / "invalid_missing_column.csv", expect_success=False)
    mismatch = run_validation(metadata_dir / "invalid_unmatched_sample.csv", expect_success=False)
    assert "missing required column" in missing.stderr, missing.stderr
    assert "metadata samples absent from FASTQ input" in mismatch.stderr, mismatch.stderr
    print("metadata validation smoke tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
