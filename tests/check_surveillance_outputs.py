#!/usr/bin/env python3
"""Smoke tests for final surveillance output rendering."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "bin" / "run_surveillance_outputs.py"


def write(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mkflupipe_outputs_") as tmp:
        work = Path(tmp)
        deps = work / "deps"
        out_dir = work / "Surveillance_Outputs"
        cons_dir = work / "consensus"
        irma_dir = work / "irma_runs" / "SAMPLE_A"
        irma_dir.mkdir(parents=True)

        write(
            deps / "blast_typing_summary.tsv",
            "sample\ttype_blast\tsubtype_HA\tsubtype_NA\thit_HA\thit_NA\n"
            "SAMPLE_A\tA\tH3\tN2\tgi|1|gb|LK054740|Influenza\tgi|2|gb|KY925263|Influenza\n",
        )
        write(
            deps / "nextclade_summary.tsv",
            "sample\tclade_display\tqc_status\n"
            "SAMPLE_A\t3C.2a1b.2a.2\tgood\n",
        )
        write(
            deps / "assembly_qc_report.tsv",
            "sample\tqc_assembly\tqc_detail\n"
            "SAMPLE_A\tPASS\tsegs:8/8|low_cov:A_NP(5x) A_PB1(7x)\n",
        )
        write(
            deps / "depth_summary.tsv",
            "sample\tsegment\tcov_mean\tcov_min\tcov_max\tpositions_covered\tref_length\n"
            "SAMPLE_A\tA_HA\t100\t10\t250\t1700\t1701\n",
        )
        write(
            deps / "coinfection_report.tsv",
            "sample\tcoinfection_status\tdetails\n"
            "SAMPLE_A\tWARN\tsegs_analyzed:8|segs_flagged:2|A_NP(23/280pos,8.21%) A_PB1(34/515pos,6.60%)\n",
        )
        write(
            deps / "antiviral_resistance.tsv",
            "sample\tgene\taa_position\twt_who\tmut_who\talt_observed\tfrequency\tdepth_total\tdrug\tsignificance\tnomenclature\n",
        )
        write(
            deps / "h5_virulence_markers.tsv",
            "sample\tmarker\tstatus\n",
        )
        write(
            deps / "all_samples_protein_mutations.tsv",
            "sample\tsegment\tgene\tmutation\n"
            "SAMPLE_A\tA_HA\tHA\tK123N\n",
        )
        write(
            deps / "validated_metadata.csv",
            "sample_name,collection_date,country,state,city\n"
            "SAMPLE_A,2026-05-24,Brazil,Alagoas,Maceio\n",
        )
        write(
            deps / "SAMPLE_A.fastp.json",
            json.dumps(
                {
                    "summary": {
                        "before_filtering": {"total_reads": 100, "total_bases": 10000},
                        "after_filtering": {"total_reads": 80, "q30_rate": 0.95, "read1_mean_length": 125},
                    },
                    "filtering_result": {"passed_filter_reads": 80},
                }
            ),
        )
        write(
            deps / "SAMPLE_A.host_depletion.stats.tsv",
            "sample_id\tlayout\tseq_type\tinput_reads\toutput_reads\tread_retention_pct\tinput_mean_len\toutput_mean_len\n"
            "SAMPLE_A\tpaired\tshort_paired\t80\t75\t93.75\t125\t124\n",
        )
        write(deps / "multiqc_report.html", "<html><body>multiqc</body></html>\n")

        phylogeny_dir = deps / "phylogeny"
        write(
            phylogeny_dir / "phylogeny_summary.tsv",
            "group\ttype\tsegment\tsubtype\tpipeline_sequences\tcontext_sequences\ttotal_sequences\tstatus\tmessage\tauspice_json\ttree_html\n"
            "A_H3_HA\tA\tHA\tH3\t1\t2\t3\tPASS\tTime-scaled Augur tree generated\tA_H3_HA/A_H3_HA.json\tA_H3_HA/A_H3_HA.html\n",
        )
        write(phylogeny_dir / "A_H3_HA" / "A_H3_HA.html", "<html><body>tree html</body></html>\n")
        write(phylogeny_dir / "A_H3_HA" / "A_H3_HA.json", "{}\n")
        write(phylogeny_dir / "A_H3_HA" / "tree.nwk", "(SAMPLE_A:0.1,CTX:0.2);\n")

        consensus = cons_dir / "SAMPLE_A.fasta"
        write(consensus, ">SAMPLE_A_4\nACTGNN\n>SAMPLE_A_6\nACTGAA\n")

        command = [
            sys.executable,
            str(SCRIPT),
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

        typing = read_tsv(out_dir / "typing_results.tsv")
        assert typing[0]["hit_blast_HA"] == "Accession: LK054740 | Influenza", typing
        assert typing[0]["hit_blast_NA"] == "Accession: KY925263 | Influenza", typing
        assert typing[0]["qc_detail"] == "Segments: 8/8\nLow Coverage: A_NP(5x), A_PB1(7x)", typing

        coinfection = read_tsv(out_dir / "coinfection" / "coinfection_report.tsv")
        assert coinfection[0]["details"] == (
            "8 segments analyzed and alerts found in: "
            "A_NP(23/280 positions, 8.21%) and A_PB1(34/515 positions, 6.60%)"
        ), coinfection

        assert (out_dir / "metadata.csv").exists()
        assert (out_dir / "phylogeny" / "A_H3_HA" / "A_H3_HA.html").exists()
        assert read_tsv(out_dir / "phylogeny" / "phylogeny_summary.tsv")[0]["tree_html"] == "A_H3_HA/A_H3_HA.html"

        dashboard = (out_dir / "surveillance_report.html").read_text(encoding="utf-8")
        assert "Open tree" in dashboard, dashboard
        assert "phylogeny/A_H3_HA/A_H3_HA.html" in dashboard, dashboard
        assert "Sample metadata" in dashboard, dashboard

        assert (out_dir / "GISAID_ready" / "gisaid_sequences.fasta").read_text(encoding="utf-8") == (
            out_dir / "multisample_consensus.fasta"
        ).read_text(encoding="utf-8")

        summary = json.loads((out_dir / "run_summary.json").read_text(encoding="utf-8"))
        assert summary["total_samples"] == 1, summary

    print("surveillance output smoke tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
