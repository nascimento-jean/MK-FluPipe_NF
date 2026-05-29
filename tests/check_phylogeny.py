#!/usr/bin/env python3
"""Smoke-test HA/NA grouping without requiring Augur binaries."""

from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_PHYLOGENY = REPO_ROOT / "bin" / "run_phylogeny.py"


def write(path: Path, text: str):
    path.write_text(text, encoding="utf-8")


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mkflupipe_phylogeny_") as tmp:
        work = Path(tmp)
        write(
            work / "blast_typing_summary.tsv",
            "sample\ttype_blast\tsubtype_HA\tsubtype_NA\thit_HA\thit_NA\n"
            "A_SAMPLE\tA\tH3\tN2\t-\t-\n"
            "B_SAMPLE\tB\t-\t-\t-\t-\n",
        )
        write(
            work / "metadata.csv",
            "sample_name,collection_date,country,state,city\n"
            "A_SAMPLE,2026-01-01,Brazil,Alagoas,Maceio\n"
            "B_SAMPLE,2026-01-02,Brazil,Alagoas,Maceio\n",
        )
        write(work / "A_SAMPLE_segment_4.fasta", ">A_SAMPLE_HA\nACTGACTGACTG\n")
        write(work / "A_SAMPLE_segment_6.fasta", ">A_SAMPLE_NA\nACTGACTGACAA\n")
        write(work / "B_SAMPLE_segment_4.fasta", ">B_SAMPLE_HA\nACTGACTGTTTG\n")
        write(work / "B_SAMPLE_segment_6.fasta", ">B_SAMPLE_NA\nACTGACTGTTAA\n")
        write(
            work / "context.fasta",
            ">CTX_A_HA\nACTGACTGACTA\n>CTX_B_NA\nACTGACTGTTAC\n",
        )
        write(
            work / "context.csv",
            "strain,collection_date,type,segment,subtype_HA,subtype_NA,country,state,source\n"
            "CTX_A_HA,2025-12-01,A,HA,H3,N2,Brazil,Estado_de_Sao_Paulo,GISAID\n"
            "CTX_B_NA,2025-12-02,B,NA,-,-,Brazil,Bahia,GISAID\n",
        )
        command = [
            sys.executable,
            str(RUN_PHYLOGENY),
            "--blast-summary",
            str(work / "blast_typing_summary.tsv"),
            "--metadata",
            str(work / "metadata.csv"),
            "--context-fasta",
            str(work / "context.fasta"),
            "--context-metadata",
            str(work / "context.csv"),
            "--output-dir",
            str(work / "phylogeny"),
            "--log-file",
            str(work / "phylogeny.log"),
            "--min-sequences",
            "99",
            str(work / "A_SAMPLE_segment_4.fasta"),
            str(work / "A_SAMPLE_segment_6.fasta"),
            str(work / "B_SAMPLE_segment_4.fasta"),
            str(work / "B_SAMPLE_segment_6.fasta"),
        ]
        result = subprocess.run(command, cwd=REPO_ROOT, capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        with (work / "phylogeny" / "phylogeny_summary.tsv").open(encoding="utf-8", newline="") as handle:
            rows = {row["group"]: row for row in csv.DictReader(handle, delimiter="\t")}
        assert set(rows) == {"A_H3_HA", "A_N2_NA", "B_HA", "B_NA"}, rows
        assert rows["A_H3_HA"]["pipeline_sequences"] == "1", rows
        assert rows["A_H3_HA"]["context_sequences"] == "1", rows
        assert rows["B_NA"]["context_sequences"] == "1", rows
        assert all(row["status"] == "SKIPPED" for row in rows.values()), rows
        with (work / "phylogeny" / "A_H3_HA" / "metadata.tsv").open(encoding="utf-8", newline="") as handle:
            metadata = {row["strain"]: row for row in csv.DictReader(handle, delimiter="\t")}
        assert metadata["A_SAMPLE"]["display_group"] == "User Sequences", metadata
        assert metadata["CTX_A_HA"]["display_group"] == "Estado de Sao Paulo", metadata
        spec = importlib.util.spec_from_file_location("run_phylogeny", RUN_PHYLOGENY)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        colors_path = work / "colors.tsv"
        module.write_display_group_colors(
            colors_path,
            [
                {"display_group": "User Sequences"},
                {"display_group": "Estado de Sao Paulo"},
                {"display_group": "Bahia"},
            ],
        )
        colors = colors_path.read_text(encoding="utf-8")
        assert "display_group\tUser Sequences\t#8B0000\n" in colors, colors
        assert "display_group\tEstado de Sao Paulo\t" in colors, colors
        assert "display_group\tBahia\t" in colors, colors
        assert "display_group\tEstado de Sao Paulo\t#8B0000\n" not in colors, colors
        color_rows = [line.split("\t") for line in colors.strip().splitlines()]
        state_colors = [row[2] for row in color_rows if row[1] != "User Sequences"]
        assert len(state_colors) == len(set(state_colors)), colors
        assert all(color != "#8B0000" for color in state_colors), colors
        html_group = work / "html_group"
        html_group.mkdir()
        write(html_group / "tree.nwk", "(A_SAMPLE:0.1,CTX_A_HA:0.2):0.0;\n")
        colors_path.replace(html_group / "colors.tsv")
        tree_html = module.write_tree_html(
            html_group,
            "A_H3_HA",
            [
                {"strain": "A_SAMPLE", "display_group": "User Sequences", "date": "2026-01-01"},
                {"strain": "CTX_A_HA", "display_group": "Estado de Sao Paulo", "date": "2025-12-01"},
            ],
        )
        assert tree_html == "A_H3_HA/A_H3_HA.html", tree_html
        html_text = (html_group / "A_H3_HA.html").read_text(encoding="utf-8")
        assert "MK Flu-Pipe A_H3_HA phylogeny" in html_text, html_text
        assert "#8B0000" in html_text, html_text
        assert "<svg" in html_text, html_text
    print("phylogeny grouping smoke tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
