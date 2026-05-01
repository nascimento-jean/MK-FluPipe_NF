#!/usr/bin/env python3

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


CANONICAL_ACCESSIONS = {
    "N1_NA": "CY121683",
    "N2_NA": "JN976844",
    "N1_PA": "CY121685",
    "N2_PA": "JN976847",
    "N1_MP": "CY121682",
    "N2_MP": "JN976843",
}


def write_log(log_path: Path, text: str):
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(text)
        if text and not text.endswith("\n"):
            handle.write("\n")


def normalize_header(fasta_text: str, key: str):
    lines = fasta_text.splitlines()
    if not lines:
        return ""
    lines[0] = f">{key}"
    return "\n".join(lines) + "\n"


def fetch_with_efetch(accession: str):
    if not shutil.which("efetch"):
        return None
    result = subprocess.run(
        ["efetch", "-db", "nuccore", "-id", accession, "-format", "fasta"],
        capture_output=True,
        text=True,
    )
    if result.returncode == 0 and result.stdout.strip():
        return result.stdout
    return None


def fetch_with_http(accession: str):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        f"efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    )
    if shutil.which("wget"):
        result = subprocess.run(["wget", "-q", "-O", "-", url], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout
    if shutil.which("curl"):
        result = subprocess.run(["curl", "-L", url], capture_output=True, text=True)
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout
    return None


def main():
    parser = argparse.ArgumentParser(description="Prepare canonical influenza references for iVar resistance analysis")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.write_text("", encoding="utf-8")

    downloaded = 0
    for key, accession in CANONICAL_ACCESSIONS.items():
        fasta_path = output_dir / f"{key}.fa"
        if fasta_path.exists() and fasta_path.stat().st_size > 0:
            continue

        write_log(log_path, f"Preparing canonical reference {key} ({accession})")
        fasta_text = fetch_with_efetch(accession)
        if fasta_text is None:
            fasta_text = fetch_with_http(accession)

        if fasta_text:
            fasta_path.write_text(normalize_header(fasta_text, key), encoding="utf-8")
            downloaded += 1
            write_log(log_path, f"  -> {fasta_path}")
        else:
            if fasta_path.exists():
                fasta_path.unlink()
            write_log(log_path, f"  -> failed to download {key} ({accession})")

    write_log(log_path, f"Canonical references downloaded in this run: {downloaded}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
