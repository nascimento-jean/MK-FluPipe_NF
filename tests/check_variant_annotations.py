#!/usr/bin/env python3
"""Smoke tests for functional and SnpEff variant annotation helpers."""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def run(command: list[str], *, env: dict[str, str] | None = None) -> None:
    result = subprocess.run(command, cwd=REPO_ROOT, env=env, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise AssertionError(f"Command failed: {' '.join(command)}")


def check_functional_annotation(tmp: Path) -> None:
    fullvar = tmp / "fullvar.tsv"
    typing = tmp / "typing.tsv"
    output = tmp / "functional_annotation.tsv"

    write_text(
        fullvar,
        "\t".join(["sample", "type", "segment", "gene", "aa_position", "ref_aa", "alt_aa", "mutation", "frequency", "depth"])
        + "\n"
        + "\t".join(["SAMPLE_A", "A", "SAMPLE_A_A_4", "HA", "190", "K", "R", "K190R", "0.42", "120"])
        + "\n"
        + "\t".join(["SAMPLE_A", "A", "SAMPLE_A_A_4", "UTR", "-", "-", "-", "-", "0.10", "50"])
        + "\n",
    )
    write_text(
        typing,
        "sample\ttype_blast\tsubtype_HA\tsubtype_NA\n"
        "SAMPLE_A\tA\tH3\tN2\n",
    )

    run(
        [
            sys.executable,
            "bin/annotate_functional_variants.py",
            "--fullvar-tsv",
            str(fullvar),
            "--typing-tsv",
            str(typing),
            "--output-tsv",
            str(output),
        ]
    )

    rows = read_tsv(output)
    assert len(rows) == 2, rows
    aa_change = rows[0]
    assert aa_change["sample"] == "SAMPLE_A"
    assert aa_change["subtype"] == "H3N2"
    assert aa_change["segment"] == "4"
    assert aa_change["segment_name"] == "HA"
    assert aa_change["effect"] == "amino_acid_change"
    assert aa_change["impact"] == "MODERATE"
    assert aa_change["annotation_source"] == "RefSeq GFF3 + iVar protein table"

    utr = rows[1]
    assert utr["effect"] == "non_coding"
    assert utr["impact"] == "MODIFIER"


def check_snpeff_tool_missing_fallback(tmp: Path) -> None:
    sample_dir = tmp / "sample_variants"
    refseq_dir = tmp / "refseq"
    output_dir = tmp / "snpeff"
    log_file = tmp / "snpeff.log"
    blast = tmp / "blast.tsv"

    write_text(
        sample_dir / "SAMPLE_A_4_ivar.tsv",
        "POS\tREF\tALT\tALT_FREQ\tTOTAL_DP\tPASS\n"
        "52\tA\tG\t0.25\t80\tTRUE\n"
        "53\tT\tC\t0.10\t40\tFALSE\n",
    )
    write_text(refseq_dir / "NC_007366.fa", ">NC_007366\n" + "A" * 120 + "\n")
    write_text(
        refseq_dir / "NC_007366.gff3",
        "##gff-version 3\n"
        "NC_007366\tRefSeq\tCDS\t1\t90\t.\t+\t0\tID=cds-HA;gene=HA;product=hemagglutinin\n",
    )
    write_text(
        blast,
        "sample\ttype_blast\tsubtype_HA\tsubtype_NA\n"
        "SAMPLE\tA\tH3\tN2\n",
    )

    env = os.environ.copy()
    fake_path = tmp / "empty_path"
    fake_path.mkdir()
    env["PATH"] = str(fake_path)

    run(
        [
            sys.executable,
            "bin/run_snpeff_annotation.py",
            "--sample-dirs",
            str(sample_dir),
            "--refseq-dir",
            str(refseq_dir),
            "--blast-summary",
            str(blast),
            "--output-dir",
            str(output_dir),
            "--log-file",
            str(log_file),
        ],
        env=env,
    )

    rows = read_tsv(output_dir / "snpeff_annotation.tsv")
    assert len(rows) == 1, rows
    row = rows[0]
    assert row["sample"] == "SAMPLE"
    assert row["type"] == "A"
    assert row["subtype"] == "H3N2"
    assert row["segment"] == "4"
    assert row["segment_name"] == "HA"
    assert row["accession"] == "NC_007366"
    assert row["pos"] == "52"
    assert row["frequency"] == "0.25"
    assert row["depth"] == "80"
    assert row["status"] == "tool_missing"
    assert "snpEff executable was not found" in row["message"]
    assert "snpEff executable was not found" in log_file.read_text(encoding="utf-8")
    assert (output_dir / "vcf" / "SAMPLE_A_4_ivar.vcf").exists()


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mkflupipe_variant_annotation_") as tmpdir:
        tmp = Path(tmpdir)
        check_functional_annotation(tmp / "functional")
        check_snpeff_tool_missing_fallback(tmp / "snpeff")
    print("variant annotation smoke tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
