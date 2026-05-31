#!/usr/bin/env python3
"""Run a tiny stubbed Nextflow workflow and verify final deliverables."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = REPO_ROOT / "tests" / "output"


def nextflow_executable() -> str:
    path_nextflow = shutil.which("nextflow")
    if path_nextflow:
        return path_nextflow

    conda_nextflow = Path.home() / "miniconda3" / "envs" / "nextflow" / "bin" / "nextflow"
    if conda_nextflow.exists():
        return str(conda_nextflow)

    raise FileNotFoundError("nextflow executable was not found in PATH or ~/miniconda3/envs/nextflow/bin")


def configure_nextflow_environment(env: dict[str, str]) -> None:
    conda_java = Path.home() / "miniconda3" / "envs" / "nextflow" / "lib" / "jvm" / "bin" / "java"
    if conda_java.exists():
        env["JAVA_CMD"] = str(conda_java)


def run_nextflow(name: str, args: list[str]) -> Path:
    output_dir = OUTPUT_ROOT / name
    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"mkflupipe_{name}_nxf_") as nxf_home:
        env = os.environ.copy()
        env["NXF_HOME"] = nxf_home
        env.setdefault("NXF_ANSI_LOG", "false")
        configure_nextflow_environment(env)
        command = [
            nextflow_executable(),
            "run",
            "main.nf",
            "-stub-run",
            "-profile",
            "linux",
            "--output_dir",
            str(output_dir),
            "--max_cpus",
            "2",
            "--queue_size",
            "1",
            *args,
        ]
        result = subprocess.run(command, cwd=REPO_ROOT, env=env, capture_output=True, text=True)
        if result.returncode != 0:
            sys.stderr.write(result.stdout)
            sys.stderr.write(result.stderr)
            raise AssertionError(f"{name} mini workflow failed with exit code {result.returncode}")

    return output_dir


def read_text(path: Path) -> str:
    assert path.exists(), f"Missing expected file: {path}"
    return path.read_text(encoding="utf-8")


def check_common_outputs(output_dir: Path, sample_id: str) -> None:
    surveillance = output_dir / "Surveillance_Outputs"
    dashboard = read_text(surveillance / "surveillance_report.html")
    run_summary = json.loads(read_text(surveillance / "run_summary.json"))
    typing = read_text(surveillance / "typing_results.tsv")
    consensus = read_text(surveillance / "multisample_consensus.fasta")
    gisaid = read_text(surveillance / "GISAID_ready" / "gisaid_sequences.fasta")
    phylogeny_summary = read_text(surveillance / "phylogeny" / "phylogeny_summary.tsv")

    assert run_summary["samples"], run_summary
    assert run_summary["samples"][0]["sample"] == sample_id, run_summary
    assert sample_id in dashboard, f"{sample_id} missing from dashboard"
    assert sample_id in typing, f"{sample_id} missing from typing results"
    assert consensus.strip(), "multisample_consensus.fasta is empty"
    assert consensus == gisaid, "GISAID FASTA should mirror multisample consensus FASTA"
    assert "segment" in phylogeny_summary, "phylogeny summary header is missing"


def main() -> int:
    short_out = run_nextflow(
        "mini-short-pipeline",
        [
            "--input_dir",
            str(REPO_ROOT / "tests" / "data" / "mini_short"),
            "--irma_module",
            "FLU-utr",
            "--seq_type",
            "short",
            "--host_depletion",
            "false",
            "--run_ivar",
            "true",
            "--run_antiviral",
            "true",
            "--run_h5_virulence",
            "true",
            "--run_fullvarcall",
            "true",
            "--metadata_csv",
            str(REPO_ROOT / "tests" / "data" / "mini_short" / "metadata.csv"),
            "--gisaid_location",
            "Brazil-AL",
            "--gisaid_year",
            "2026",
            "--run_phylogeny",
            "true",
            "--phylogeny_min_sequences",
            "99",
        ],
    )
    check_common_outputs(short_out, "MINI-SHORT-001")
    short_dashboard = read_text(short_out / "Surveillance_Outputs" / "surveillance_report.html")
    assert "MINI-SHORT-001_S1" not in short_dashboard, "Illumina suffix leaked into sample name"
    assert (short_out / "variant_calls" / "MINI-SHORT-001").exists(), "short iVar output missing"
    assert (short_out / "full_variant_calls" / "MINI-SHORT-001").exists(), "short full variant output missing"
    assert (short_out / "Surveillance_Outputs" / "antiviral_resistance" / "antiviral_resistance.tsv").exists()
    assert (short_out / "Surveillance_Outputs" / "h5_virulence" / "h5_virulence_markers.tsv").exists()
    assert (short_out / "Surveillance_Outputs" / "functional_annotation" / "functional_annotation.tsv").exists()

    long_out = run_nextflow(
        "mini-long-pipeline",
        [
            "--input_dir",
            str(REPO_ROOT / "tests" / "data" / "mini_long"),
            "--irma_module",
            "FLU-minion",
            "--seq_type",
            "long",
            "--host_depletion",
            "false",
            "--run_medaka",
            "true",
            "--run_antiviral",
            "true",
            "--run_h5_virulence",
            "false",
            "--run_fullvarcall",
            "false",
            "--metadata_csv",
            str(REPO_ROOT / "tests" / "data" / "mini_long" / "metadata.csv"),
            "--gisaid_location",
            "Brazil-AL",
            "--gisaid_year",
            "2026",
            "--run_phylogeny",
            "true",
            "--phylogeny_min_sequences",
            "99",
        ],
    )
    check_common_outputs(long_out, "MINI_LONG")
    assert (long_out / "variant_calls" / "MINI_LONG").exists(), "long Medaka output missing"
    assert (long_out / "variant_calls_canonical_long" / "MINI_LONG").exists(), "long canonical Medaka output missing"
    assert (long_out / "Surveillance_Outputs" / "antiviral_resistance" / "antiviral_resistance.tsv").exists()

    print("mini pipeline integration tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
