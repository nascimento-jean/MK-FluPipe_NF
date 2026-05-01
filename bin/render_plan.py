#!/usr/bin/env python3

import argparse
import csv
import json
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Render a migration execution plan")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seq-type", required=True)
    parser.add_argument("--irma-module", required=True)
    parser.add_argument("--run-fastqc", required=True)
    parser.add_argument("--host-depletion", required=True)
    parser.add_argument("--run-ivar", required=True)
    parser.add_argument("--run-medaka", required=True)
    parser.add_argument("--run-antiviral", required=True)
    parser.add_argument("--run-h5-virulence", required=True)
    parser.add_argument("--run-fullvarcall", required=True)
    parser.add_argument("--run-legacy-bridge", required=True)
    return parser.parse_args()


def as_bool(text: str) -> bool:
    return str(text).strip().lower() == "true"


def stage_rows(args):
    stages = [
        ("discover_samples", True, "Create bootstrap samplesheet from FASTQ input"),
        ("plan_run", True, "Write execution plan and migration summary"),
        ("fastqc", as_bool(args.run_fastqc), "Raw-read QC"),
        ("host_depletion", as_bool(args.host_depletion), "Remove human host reads before assembly"),
        ("irma_assembly", True, "Core influenza assembly with IRMA"),
        ("blast_typing", True, "HA/NA subtype support with BLAST"),
        ("nextclade", True, "Clade assignment and QC"),
        ("ivar", as_bool(args.run_ivar), "Canonical short-read variant calling"),
        ("medaka", as_bool(args.run_medaka), "ONT variant calling"),
        ("coinfection", True, "Minority allele and mixed infection analysis"),
        ("antiviral", as_bool(args.run_antiviral), "Antiviral resistance screening"),
        ("h5_virulence", as_bool(args.run_h5_virulence), "Conditional H5 marker screening"),
        ("fullvarcall", as_bool(args.run_fullvarcall), "All-segment RefSeq + GFF3 protein variant calling"),
        ("surveillance_outputs", True, "Consolidated TSV/JSON/HTML/GISAID outputs"),
        ("legacy_bridge", as_bool(args.run_legacy_bridge), "Run current Bash pipeline from Nextflow while migrating"),
    ]
    return stages


def main():
    args = parse_args()
    samplesheet = Path(args.samplesheet)
    summary_json = Path(args.summary)
    output = Path(args.output)

    rows = list(csv.DictReader(samplesheet.open(encoding="utf-8")))
    summary = json.loads(summary_json.read_text(encoding="utf-8"))
    enabled_stages = [row for row in stage_rows(args) if row[1]]

    lines = []
    lines.append("MK Flu-Pipe Nextflow MVP")
    lines.append("")
    lines.append(f"IRMA module: {args.irma_module}")
    lines.append(f"Requested seq_type: {args.seq_type}")
    lines.append(f"Detected samples: {summary['sample_count']}")
    lines.append(f"Paired samples: {summary['paired_count']}")
    lines.append(f"Single samples: {summary['single_count']}")
    lines.append("")
    lines.append("Samples")
    for row in rows:
        mate = row["r2"] if row["r2"] else "-"
        lines.append(
            f"- {row['sample_id']} | layout={row['layout']} | seq_type={row['seq_type']} | r1={row['r1']} | r2={mate}"
        )
    lines.append("")
    lines.append("Enabled stages")
    for name, _, description in enabled_stages:
        lines.append(f"- {name}: {description}")
    lines.append("")
    lines.append("Migration note")
    lines.append("- This MVP currently bootstraps sample discovery and planning.")
    lines.append("- Full modular migration will replace the legacy Bash logic step by step.")

    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

