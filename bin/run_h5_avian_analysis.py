#!/usr/bin/env python3

"""Optional H5 avian-focused analysis with GenoFLU and FluMut.

The wrapper is intentionally defensive: missing tools or tool failures are
reported in TSV/log outputs instead of stopping the whole surveillance run.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import textwrap
from collections import defaultdict
from pathlib import Path


SEGMENT_NAMES = {
    1: "PB2",
    2: "PB1",
    3: "PA",
    4: "HA",
    5: "NP",
    6: "NA",
    7: "MP",
    8: "NS",
}

SUMMARY_FIELDS = [
    "sample",
    "type",
    "subtype_HA",
    "subtype_NA",
    "segments_available",
    "selected_for_h5_avian",
    "genoflu_status",
    "genoflu_genotype",
    "genoflu_message",
    "flumut_status",
    "flumut_markers_detected",
    "flumut_message",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def parse_fasta(path: Path):
    header = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if header is not None:
        yield header, "".join(chunks)


def infer_sample_segment(path: Path) -> tuple[str, int] | None:
    stem = path.name
    marker = "_segment_"
    if marker not in stem:
        return None
    sample, tail = stem.rsplit(marker, 1)
    segment_txt = tail.split(".", 1)[0]
    if not segment_txt.isdigit():
        return None
    segment = int(segment_txt)
    if segment not in SEGMENT_NAMES:
        return None
    return sample, segment


def collect_segment_files(paths: list[Path]) -> dict[str, dict[int, Path]]:
    by_sample: dict[str, dict[int, Path]] = defaultdict(dict)
    for path in paths:
        inferred = infer_sample_segment(path)
        if inferred is None:
            continue
        sample, segment = inferred
        if path.exists() and path.stat().st_size > 0:
            by_sample[sample][segment] = path
    return by_sample


def is_h5_candidate(row: dict[str, str]) -> bool:
    flu_type = (row.get("type_blast") or row.get("type") or "").strip().upper()
    subtype_ha = (row.get("subtype_HA") or row.get("subtype_HA".lower()) or "").strip().upper()
    subtype_na = (row.get("subtype_NA") or row.get("subtype_NA".lower()) or "").strip().upper()
    return flu_type == "A" and subtype_ha.startswith("H5") and subtype_na == "N1"


def write_sample_fasta(sample: str, segment_map: dict[int, Path], output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    out_fasta = output_dir / f"{sample}.h5_avian.fasta"
    with out_fasta.open("w", encoding="utf-8") as handle:
        for segment in range(1, 9):
            path = segment_map.get(segment)
            if not path:
                continue
            for _header, seq in parse_fasta(path):
                if not seq:
                    continue
                handle.write(f">{sample}_{SEGMENT_NAMES[segment]}_segment_{segment}\n{seq}\n")
    return out_fasta


def run_command(cmd: list[str], cwd: Path) -> tuple[int, str, str]:
    proc = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, check=False)
    return proc.returncode, proc.stdout, proc.stderr


def parse_genoflu_stdout(stdout: str) -> tuple[str, str]:
    text = " ".join(line.strip() for line in stdout.splitlines() if line.strip())
    genotype = ""
    if "Genotype -->" in text:
        genotype = text.split("Genotype -->", 1)[1].strip()
    elif "Genotype:" in text:
        genotype = text.split("Genotype:", 1)[1].strip()
    if genotype:
        genotype = genotype.split("\t", 1)[0].strip()
    return genotype or "-", text or "-"


def run_genoflu(sample: str, fasta: Path, work_dir: Path) -> tuple[str, str, str]:
    executable = shutil.which("genoflu.py") or shutil.which("genoflu")
    if not executable:
        return "tool_missing", "-", "GenoFLU executable not found in PATH"

    code, stdout, stderr = run_command([executable, "-f", str(fasta)], work_dir)
    (work_dir / f"{sample}.genoflu.stdout.txt").write_text(stdout, encoding="utf-8")
    (work_dir / f"{sample}.genoflu.stderr.txt").write_text(stderr, encoding="utf-8")
    genotype, message = parse_genoflu_stdout(stdout)
    if code != 0:
        return "failed", genotype, (stderr.strip() or message or f"GenoFLU exited with status {code}")
    return "ok", genotype, message


def normalize_flumut_table(path: Path, sample: str, fallback_fields: list[str]) -> tuple[list[dict[str, str]], list[str]]:
    if not path.exists() or path.stat().st_size == 0:
        return [], ["sample"] + fallback_fields
    rows = read_tsv(path)
    if not rows:
        return [], ["sample"] + fallback_fields
    fields = ["sample"] + [field for field in rows[0].keys() if field != "sample"]
    out_rows = []
    for row in rows:
        new_row = {"sample": row.get("sample", sample) or sample}
        for field in fields:
            if field != "sample":
                new_row[field] = row.get(field, "")
        out_rows.append(new_row)
    return out_rows, fields


def run_flumut(sample: str, fasta: Path, work_dir: Path) -> tuple[str, str, str, list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[str], list[str], list[str]]:
    executable = shutil.which("flumut")
    if not executable:
        return (
            "tool_missing",
            "0",
            "FluMut executable not found in PATH",
            [],
            [],
            [],
            ["sample", "marker"],
            ["sample", "mutation"],
            ["sample", "reference"],
        )

    markers = work_dir / f"{sample}.flumut.markers.tsv"
    mutations = work_dir / f"{sample}.flumut.mutations.tsv"
    literature = work_dir / f"{sample}.flumut.literature.tsv"
    cmd = [
        executable,
        "-m",
        str(markers),
        "-M",
        str(mutations),
        "-l",
        str(literature),
        str(fasta),
    ]
    code, stdout, stderr = run_command(cmd, work_dir)
    (work_dir / f"{sample}.flumut.stdout.txt").write_text(stdout, encoding="utf-8")
    (work_dir / f"{sample}.flumut.stderr.txt").write_text(stderr, encoding="utf-8")

    marker_rows, marker_fields = normalize_flumut_table(markers, sample, ["marker"])
    mutation_rows, mutation_fields = normalize_flumut_table(mutations, sample, ["mutation"])
    literature_rows, literature_fields = normalize_flumut_table(literature, sample, ["reference"])
    if code != 0:
        return "failed", str(len(marker_rows)), stderr.strip() or f"FluMut exited with status {code}", marker_rows, mutation_rows, literature_rows, marker_fields, mutation_fields, literature_fields
    return "ok", str(len(marker_rows)), stdout.strip() or "FluMut completed", marker_rows, mutation_rows, literature_rows, marker_fields, mutation_fields, literature_fields


def main() -> None:
    parser = argparse.ArgumentParser(description="Run optional GenoFLU/FluMut H5N1 avian analysis")
    parser.add_argument("--blast-summary", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("segment_files", nargs="*")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    work_dir = output_dir / "_work"
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    log_lines = ["MK Flu-Pipe H5 avian analysis", "Tools: GenoFLU + FluMut"]

    blast_rows = read_tsv(Path(args.blast_summary))
    segment_files = [Path(item) for item in args.segment_files]
    segments_by_sample = collect_segment_files(segment_files)

    summary_rows: list[dict[str, str]] = []
    all_marker_rows: list[dict[str, str]] = []
    all_mutation_rows: list[dict[str, str]] = []
    all_literature_rows: list[dict[str, str]] = []
    marker_fields = ["sample", "marker"]
    mutation_fields = ["sample", "mutation"]
    literature_fields = ["sample", "reference"]

    for row in blast_rows:
        sample = (row.get("sample") or "").strip()
        if not sample:
            continue
        segment_map = segments_by_sample.get(sample, {})
        subtype_ha = (row.get("subtype_HA") or "-").strip() or "-"
        subtype_na = (row.get("subtype_NA") or "-").strip() or "-"
        selected = is_h5_candidate(row)
        base = {
            "sample": sample,
            "type": (row.get("type_blast") or "-").strip() or "-",
            "subtype_HA": subtype_ha,
            "subtype_NA": subtype_na,
            "segments_available": ",".join(str(seg) for seg in sorted(segment_map)) or "-",
            "selected_for_h5_avian": "yes" if selected else "no",
            "genoflu_status": "skipped",
            "genoflu_genotype": "-",
            "genoflu_message": "Sample is not influenza A/H5N1",
            "flumut_status": "skipped",
            "flumut_markers_detected": "0",
            "flumut_message": "Sample is not influenza A/H5N1",
        }
        if not selected:
            summary_rows.append(base)
            continue

        if not segment_map:
            base.update(
                {
                    "genoflu_status": "skipped",
                    "genoflu_message": "No segment FASTA files available",
                    "flumut_status": "skipped",
                    "flumut_message": "No segment FASTA files available",
                }
            )
            summary_rows.append(base)
            continue

        sample_dir = work_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        sample_fasta = write_sample_fasta(sample, segment_map, sample_dir)
        shutil.copy2(sample_fasta, output_dir / f"{sample}.h5_avian.fasta")
        log_lines.append(f"Analyzing {sample}: HA={subtype_ha} NA={subtype_na} segments={base['segments_available']}")

        genoflu_status, genotype, genoflu_message = run_genoflu(sample, sample_fasta, sample_dir)
        (
            flumut_status,
            markers_count,
            flumut_message,
            marker_rows,
            mutation_rows,
            literature_rows,
            marker_fields,
            mutation_fields,
            literature_fields,
        ) = run_flumut(sample, sample_fasta, sample_dir)

        all_marker_rows.extend(marker_rows)
        all_mutation_rows.extend(mutation_rows)
        all_literature_rows.extend(literature_rows)
        base.update(
            {
                "genoflu_status": genoflu_status,
                "genoflu_genotype": genotype,
                "genoflu_message": genoflu_message,
                "flumut_status": flumut_status,
                "flumut_markers_detected": markers_count,
                "flumut_message": flumut_message,
            }
        )
        summary_rows.append(base)

    if not summary_rows:
        log_lines.append("No BLAST typing rows were available for H5 avian analysis.")

    write_tsv(output_dir / "h5_avian_summary.tsv", summary_rows, SUMMARY_FIELDS)
    write_tsv(output_dir / "genoflu_summary.tsv", summary_rows, SUMMARY_FIELDS[:7] + ["genoflu_genotype", "genoflu_message"])
    write_tsv(output_dir / "flumut_markers.tsv", all_marker_rows, marker_fields)
    write_tsv(output_dir / "flumut_mutations.tsv", all_mutation_rows, mutation_fields)
    write_tsv(output_dir / "flumut_literature.tsv", all_literature_rows, literature_fields)

    Path(args.log_file).parent.mkdir(parents=True, exist_ok=True)
    Path(args.log_file).write_text(
        "\n".join(log_lines)
        + "\n\n"
        + textwrap.dedent(
            """\
            Notes:
            - GenoFLU is optimized for HPAI H5 clade 2.3.4.4b genotyping.
            - FluMut is designed for mutation surveillance in A(H5N1) nucleotide sequences.
            - Missing tools are reported as tool_missing instead of failing the pipeline.
            """
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
