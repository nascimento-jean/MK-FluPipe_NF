#!/usr/bin/env python3
"""Smoke tests for optional H5 avian GenoFLU/FluMut wrapper."""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "bin" / "run_h5_avian_analysis.py"


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_executable(path: Path, text: str) -> None:
    write(path, text)
    path.chmod(path.stat().st_mode | 0o111)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="mkflupipe_h5avian_") as tmp:
        work = Path(tmp)
        tools = work / "tools"
        deps = work / "deps"
        out_dir = work / "h5_avian"

        make_executable(
            tools / "genoflu.py",
            "#!/usr/bin/env python3\n"
            "print('GenoFLU synthetic run')\n"
            "print('Genotype --> GsGD-like reassortant')\n",
        )
        make_executable(
            tools / "flumut",
            "#!/usr/bin/env python3\n"
            "import argparse\n"
            "parser = argparse.ArgumentParser()\n"
            "parser.add_argument('-m', dest='markers')\n"
            "parser.add_argument('-M', dest='mutations')\n"
            "parser.add_argument('-l', dest='literature')\n"
            "parser.add_argument('fasta')\n"
            "args = parser.parse_args()\n"
            "open(args.markers, 'w').write('marker\\tgene\\nPB2_E627K\\tPB2\\n')\n"
            "open(args.mutations, 'w').write('mutation\\tgene\\nE627K\\tPB2\\n')\n"
            "open(args.literature, 'w').write('reference\\tpmid\\nExample reference\\t123\\n')\n"
            "print('FluMut synthetic run')\n",
        )

        write(
            deps / "blast_typing_summary.tsv",
            "sample\ttype_blast\tsubtype_HA\tsubtype_NA\n"
            "SAMPLE_H5\tA\tH5\tN1\n"
            "SAMPLE_H3\tA\tH3\tN2\n",
        )
        write(deps / "SAMPLE_H5_segment_4.fasta", ">SAMPLE_H5_4\nACTGACTG\n")
        write(deps / "SAMPLE_H5_segment_6.fasta", ">SAMPLE_H5_6\nACTGACAA\n")
        write(deps / "SAMPLE_H3_segment_4.fasta", ">SAMPLE_H3_4\nACTGACTG\n")

        env = os.environ.copy()
        env["PATH"] = str(tools) + os.pathsep + env.get("PATH", "")
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--blast-summary",
                str(deps / "blast_typing_summary.tsv"),
                "--output-dir",
                str(out_dir),
                "--log-file",
                str(work / "reports" / "h5_avian_analysis.log"),
                str(deps / "SAMPLE_H5_segment_4.fasta"),
                str(deps / "SAMPLE_H5_segment_6.fasta"),
                str(deps / "SAMPLE_H3_segment_4.fasta"),
            ],
            cwd=REPO_ROOT,
            env=env,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr

        summary = read_tsv(out_dir / "h5_avian_summary.tsv")
        by_sample = {row["sample"]: row for row in summary}
        assert by_sample["SAMPLE_H5"]["selected_for_h5_avian"] == "yes", summary
        assert by_sample["SAMPLE_H5"]["genoflu_status"] == "ok", summary
        assert by_sample["SAMPLE_H5"]["genoflu_genotype"] == "GsGD-like reassortant", summary
        assert by_sample["SAMPLE_H5"]["flumut_status"] == "ok", summary
        assert by_sample["SAMPLE_H5"]["flumut_markers_detected"] == "1", summary
        assert by_sample["SAMPLE_H3"]["selected_for_h5_avian"] == "no", summary
        assert by_sample["SAMPLE_H3"]["genoflu_status"] == "skipped", summary

        markers = read_tsv(out_dir / "flumut_markers.tsv")
        assert markers[0]["sample"] == "SAMPLE_H5", markers
        assert markers[0]["marker"] == "PB2_E627K", markers
        assert (out_dir / "SAMPLE_H5.h5_avian.fasta").exists()

        result_missing = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--blast-summary",
                str(deps / "blast_typing_summary.tsv"),
                "--output-dir",
                str(work / "h5_avian_missing_tools"),
                "--log-file",
                str(work / "reports" / "h5_avian_missing_tools.log"),
                str(deps / "SAMPLE_H5_segment_4.fasta"),
            ],
            cwd=REPO_ROOT,
            env={**env, "PATH": ""},
            capture_output=True,
            text=True,
        )
        assert result_missing.returncode == 0, result_missing.stderr
        missing_summary = read_tsv(work / "h5_avian_missing_tools" / "h5_avian_summary.tsv")
        assert missing_summary[0]["genoflu_status"] == "tool_missing", missing_summary
        assert missing_summary[0]["flumut_status"] == "tool_missing", missing_summary

    print("H5 avian GenoFLU/FluMut smoke tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
