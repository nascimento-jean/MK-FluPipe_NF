#!/usr/bin/env python3
"""Smoke tests for sample discovery without running containerized tools."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DISCOVER = REPO_ROOT / "bin" / "discover_samples.py"


def run_discover(input_dir: Path, seq_type: str) -> tuple[list[dict[str, str]], dict]:
    with tempfile.TemporaryDirectory(prefix="mkflupipe_discover_") as tmp:
        tmpdir = Path(tmp)
        samplesheet = tmpdir / "samplesheet.csv"
        summary = tmpdir / "summary.json"
        subprocess.run(
            [
                sys.executable,
                str(DISCOVER),
                "--input-dir",
                str(input_dir),
                "--seq-type",
                seq_type,
                "--samplesheet",
                str(samplesheet),
                "--summary",
                str(summary),
            ],
            cwd=REPO_ROOT,
            check=True,
        )
        with samplesheet.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        with summary.open(encoding="utf-8") as handle:
            metadata = json.load(handle)
    return rows, metadata


def assert_short_alias() -> None:
    rows, metadata = run_discover(REPO_ROOT / "tests" / "data" / "short_paired", "short")
    assert len(rows) == 1, rows
    assert rows[0]["sample_id"] == "SAMPLE_A", rows
    assert rows[0]["layout"] == "paired", rows
    assert rows[0]["seq_type"] == "short_paired", rows
    assert metadata["sample_count"] == 1, metadata
    assert metadata["paired_count"] == 1, metadata


def assert_long_discovery() -> None:
    rows, metadata = run_discover(REPO_ROOT / "tests" / "data" / "long_single", "long")
    assert len(rows) == 1, rows
    assert rows[0]["sample_id"] == "LONG_A", rows
    assert rows[0]["layout"] == "single", rows
    assert rows[0]["seq_type"] == "long", rows
    assert metadata["sample_count"] == 1, metadata
    assert metadata["single_count"] == 1, metadata


def main() -> int:
    assert_short_alias()
    assert_long_discovery()
    print("discover_samples smoke tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
