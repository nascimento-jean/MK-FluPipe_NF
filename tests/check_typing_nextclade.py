#!/usr/bin/env python3
"""Smoke tests for BLAST typing and Nextclade parsing helpers."""

from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


blast = load_module("run_blast_typing", REPO_ROOT / "bin" / "run_blast_typing.py")
nextclade = load_module("run_nextclade_batch", REPO_ROOT / "bin" / "run_nextclade_batch.py")


def check_blast_typing_helpers(tmp: Path) -> None:
    tmp.mkdir(parents=True, exist_ok=True)
    assert blast.parse_blast_hit("Influenza A virus (A/Brazil/1/2026(H3N2)) hemagglutinin") == ("A", "H3N2")
    assert blast.parse_blast_hit("Influenza B virus B/Victoria lineage segment 4") == ("B", "Victoria")
    assert blast.parse_blast_hit("Influenza B virus B/Yamagata lineage segment 4") == ("B", "Yamagata")
    assert blast.parse_blast_hit("Unclassified orthomyxovirus") == ("Not determined", "nd")

    blast_out = tmp / "blast.tsv"
    blast_out.write_text(
        "q1\tlow\t98.0\t900\t1e-50\t200\tInfluenza A virus H1N1\n"
        "q1\tbest\t99.0\t800\t1e-20\t300\tInfluenza A virus H3N2\n"
        "q1\tbad\tmissing\tcolumns\n",
        encoding="utf-8",
    )
    results = blast.parse_blast_results(blast_out)
    assert len(results) == 2, results
    assert blast.select_best_hit(results)["sseqid"] == "best"

    sample, segment = blast.sample_from_filename(Path("SAMPLE-001_segment_6.fasta"))
    assert (sample, segment) == ("SAMPLE-001", "6")
    assert blast.sample_from_filename(Path("SAMPLE-001.fasta")) == (None, None)


def check_nextclade_helpers(tmp: Path) -> None:
    tmp.mkdir(parents=True, exist_ok=True)
    assert nextclade.get_nextclade_key("A", "H1") == "H1N1"
    assert nextclade.get_nextclade_key("A", "H3") == "H3N2"
    assert nextclade.get_nextclade_key("A", "H5") == "H5"
    assert nextclade.get_nextclade_key("B", "-") == "B"
    assert nextclade.get_nextclade_key("Not determined", "-") == ""

    a_summary = tmp / "a_nextclade.tsv"
    a_summary.write_text(
        "seqName\tclade\tshort-clade\tlegacy-clade\tqc.overallStatus\n"
        "SAMPLE_A\t3C.2a1b.2a.2\t3C.2a1b\t3C.2a\tgood\n",
        encoding="utf-8",
    )
    parsed_a = nextclade.parse_nextclade_summary(a_summary, "A")
    assert parsed_a["clade_display"] == "3C.2a1b.2a.2/3C.2a"
    assert parsed_a["qc_status"] == "good"

    b_summary = tmp / "b_nextclade.tsv"
    b_summary.write_text(
        "seqName\tclade\tshort-clade\tlegacy-clade-vic\tlegacy-clade-yam\tlineage\tqc.overallStatus\n"
        "SAMPLE_B\tV1A.3a.2\t-\tVic\t-\tVictoria\tmediocre\n",
        encoding="utf-8",
    )
    parsed_b = nextclade.parse_nextclade_summary(b_summary, "B")
    assert parsed_b["clade_display"] == "V1A.3a.2/Vic"
    assert parsed_b["qc_status"] == "mediocre"

    missing_summary = tmp / "missing_nextclade.tsv"
    missing_summary.write_text(
        "seqName\tclade\tshort-clade\tlegacy-clade\tqc.overallStatus\n"
        "SAMPLE_ND\tN/A\tN/A\tN/A\tN/A\n",
        encoding="utf-8",
    )
    parsed_missing = nextclade.parse_nextclade_summary(missing_summary, "A")
    assert parsed_missing["clade_display"] == "-"
    assert parsed_missing["qc_status"] == "N/A"


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mkflupipe_typing_nextclade_") as tmpdir:
        tmp = Path(tmpdir)
        check_blast_typing_helpers(tmp / "blast")
        check_nextclade_helpers(tmp / "nextclade")
    print("typing and Nextclade helper tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
