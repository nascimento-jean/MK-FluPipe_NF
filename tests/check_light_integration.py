#!/usr/bin/env python3
"""Lightweight integration test across discovery, metadata, and final outputs."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DISCOVER = REPO_ROOT / "bin" / "discover_samples.py"
VALIDATE_METADATA = REPO_ROOT / "bin" / "validate_metadata.py"
SURVEILLANCE_OUTPUTS = REPO_ROOT / "bin" / "run_surveillance_outputs.py"


def write(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def read_csv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mkflupipe_light_integration_") as tmp:
        work = Path(tmp)
        input_dir = work / "fastq"
        deps = work / "deps"
        consensus_dir = work / "consensus"
        irma_dir = work / "irma_runs" / "261118000051"
        deps.mkdir()
        irma_dir.mkdir(parents=True)

        write(input_dir / "261118000051_S40_L001_R1_001.fastq.gz", "synthetic-r1\n")
        write(input_dir / "261118000051_S40_L001_R2_001.fastq.gz", "synthetic-r2\n")

        samplesheet = work / "samplesheet.csv"
        discovery_summary = work / "discovery_summary.json"
        subprocess.run(
            [
                sys.executable,
                str(DISCOVER),
                "--input-dir",
                str(input_dir),
                "--seq-type",
                "short",
                "--samplesheet",
                str(samplesheet),
                "--summary",
                str(discovery_summary),
            ],
            cwd=REPO_ROOT,
            check=True,
        )
        samples = read_csv(samplesheet)
        assert samples[0]["sample_id"] == "261118000051", samples
        assert samples[0]["r1"].endswith("261118000051_S40_L001_R1_001.fastq.gz"), samples

        metadata_csv = work / "metadata.csv"
        write(
            metadata_csv,
            "sample_name,collection_date,country,state,city\n"
            "261118000051,2026-02-20,Brazil,Alagoas,Maceio\n",
        )
        validated_metadata = deps / "validated_metadata.csv"
        subprocess.run(
            [
                sys.executable,
                str(VALIDATE_METADATA),
                "--metadata",
                str(metadata_csv),
                "--samplesheet",
                str(samplesheet),
                "--output",
                str(validated_metadata),
                "--summary",
                str(work / "metadata_validation.json"),
                "--validated-samplesheet",
                str(work / "validated_samplesheet.csv"),
            ],
            cwd=REPO_ROOT,
            check=True,
        )
        assert read_csv(validated_metadata)[0]["sample_name"] == "261118000051"

        sample = "261118000051"
        write(
            deps / "blast_typing_summary.tsv",
            "sample\ttype_blast\tsubtype_HA\tsubtype_NA\thit_HA\thit_NA\n"
            f"{sample}\tA\tH3\tN2\tgi|1|gb|LK054740|Influenza\tgi|2|gb|KY925263|Influenza\n",
        )
        write(deps / "nextclade_summary.tsv", f"sample\tclade_display\tqc_status\n{sample}\t3C.2a1b.2a.2\tgood\n")
        write(deps / "assembly_qc_report.tsv", f"sample\tqc_assembly\tqc_detail\n{sample}\tPASS\tsegs:8/8\n")
        write(deps / "depth_summary.tsv", f"sample\tsegment\tcov_mean\tcov_min\tcov_max\tpositions_covered\tref_length\n{sample}\tA_HA\t95\t20\t150\t1700\t1701\n")
        write(deps / "coinfection_report.tsv", f"sample\tcoinfection_status\tdetails\n{sample}\tOK\tsegs_analyzed:8|segs_flagged:0|no_alert\n")
        write(deps / "antiviral_resistance.tsv", "sample\tgene\taa_position\twt_who\tmut_who\talt_observed\tfrequency\tdepth_total\tdrug\tsignificance\tnomenclature\n")
        write(deps / "h5_virulence_markers.tsv", "sample\tmarker\tstatus\n")
        write(deps / "all_samples_protein_mutations.tsv", "sample\tsegment\tgene\tmutation\n")
        write(
            deps / f"{sample}.fastp.json",
            json.dumps(
                {
                    "summary": {
                        "before_filtering": {"total_reads": 50},
                        "after_filtering": {"total_reads": 45, "q30_rate": 0.97, "read1_mean_length": 140},
                    },
                    "filtering_result": {"passed_filter_reads": 45},
                }
            ),
        )
        write(
            deps / f"{sample}.host_depletion.stats.tsv",
            "sample_id\tlayout\tseq_type\tinput_reads\toutput_reads\tread_retention_pct\tinput_mean_len\toutput_mean_len\n"
            f"{sample}\tpaired\tshort_paired\t45\t43\t95.56\t140\t139\n",
        )
        write(deps / "multiqc_report.html", "<html><body>multiqc</body></html>\n")

        consensus = consensus_dir / f"{sample}.fasta"
        write(consensus, f">{sample}_4\nACTGACTG\n>{sample}_6\nACTGACAA\n")

        out_dir = work / "Surveillance_Outputs"
        command = [
            sys.executable,
            str(SURVEILLANCE_OUTPUTS),
            "--output-dir",
            str(out_dir),
            "--log-file",
            str(work / "reports" / "surveillance_outputs.log"),
            "--irma-module",
            "FLU-utr",
            "--pipeline-version",
            "0.1.0",
            "--gisaid-location",
            "Brazil-AL",
            "--gisaid-year",
            "2026",
            "--dependencies",
            *[str(path) for path in sorted(deps.iterdir())],
            "--consensus-fastas",
            str(consensus),
            "--irma-dirs",
            str(irma_dir),
        ]
        result = subprocess.run(command, cwd=REPO_ROOT, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr

        run_summary = json.loads((out_dir / "run_summary.json").read_text(encoding="utf-8"))
        assert run_summary["samples"][0]["sample"] == sample, run_summary
        dashboard = (out_dir / "surveillance_report.html").read_text(encoding="utf-8")
        assert "261118000051" in dashboard, dashboard
        assert "261118000051_S40" not in dashboard, dashboard
        assert "Sample metadata" in dashboard, dashboard

    print("light integration smoke test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
