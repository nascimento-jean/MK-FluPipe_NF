#!/usr/bin/env python3

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


NEXTCLADE_DATASETS = {
    "H1N1": "nextstrain/flu/h1n1pdm/ha/MW626062",
    "H3N2": "nextstrain/flu/h3n2/ha/EPI1857216",
    "B": "nextstrain/flu/b/ha/KX058884",
    "H5": "community/moncla-lab/iav-h5/ha/all-clades",
}


def get_nextclade_key(flu_type: str, subtype_ha: str):
    subtype = (subtype_ha or "").upper()
    if flu_type == "A":
        if subtype.startswith("H1N1") or subtype == "H1":
            return "H1N1"
        if subtype.startswith("H3N2") or subtype == "H3":
            return "H3N2"
        if subtype.startswith("H5"):
            return "H5"
        return ""
    if flu_type == "B":
        return "B"
    return ""


def days_since_modified(path: Path):
    if not path.exists():
        return 999
    return int((time.time() - path.stat().st_mtime) / 86400)


def ensure_dataset(dataset_key: str, datasets_root: Path, max_days: int, log_path: Path):
    dataset_name = NEXTCLADE_DATASETS[dataset_key]
    dataset_dir = datasets_root / dataset_key
    timestamp = dataset_dir / ".last_update"

    needs_refresh = (
        (not dataset_dir.exists())
        or (not (dataset_dir / "pathogen.json").exists())
        or (days_since_modified(timestamp) > max_days)
    )
    if not needs_refresh:
        return dataset_dir

    dataset_dir.mkdir(parents=True, exist_ok=True)
    command = [
        "nextclade",
        "dataset",
        "get",
        "--name",
        dataset_name,
        "--output-dir",
        str(dataset_dir),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    with log_path.open("a", encoding="utf-8") as handle:
        if result.stdout:
            handle.write(result.stdout)
            if not result.stdout.endswith("\n"):
                handle.write("\n")
        if result.stderr:
            handle.write(result.stderr)
            if not result.stderr.endswith("\n"):
                handle.write("\n")
    if result.returncode != 0:
        raise RuntimeError(f"nextclade dataset get failed for {dataset_key}")

    timestamp.write_text(time.strftime("%Y-%m-%d"), encoding="utf-8")
    return dataset_dir


def parse_nextclade_summary(tsv_path: Path, flu_type: str):
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return {
            "clade_display": "â€”",
            "qc_status": "N/A",
        }

    with tsv_path.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, None)
        if not row:
            return {
                "clade_display": "â€”",
                "qc_status": "N/A",
            }

    clade = row.get("clade", "") or ""
    legacy = row.get("legacy-clade", "") or ""
    legacy_vic = row.get("legacy-clade-vic", "") or ""
    legacy_yam = row.get("legacy-clade-yam", "") or ""
    lineage = row.get("lineage", "") or ""
    qc_status = row.get("qc.overallStatus", "") or "N/A"

    if clade in {"", "N/A", "ERRO"}:
        return {
            "clade_display": "â€”",
            "qc_status": qc_status,
        }

    if flu_type == "B":
        b_legacy = ""
        if legacy_vic not in {"", "unassigned", "N/A"}:
            b_legacy = legacy_vic
        elif legacy_yam not in {"", "unassigned", "N/A"}:
            b_legacy = legacy_yam

        if b_legacy:
            clade_display = f"{b_legacy}/{lineage}" if lineage not in {"", "N/A", b_legacy} else b_legacy
        else:
            clade_display = f"{clade}/{lineage}" if lineage not in {"", "N/A"} else clade
    else:
        if clade not in {"", "N/A", "unassigned"} and legacy not in {"", "N/A", "unassigned"}:
            clade_display = f"{clade}/{legacy}"
        elif clade not in {"", "N/A", "unassigned"}:
            clade_display = clade
        else:
            clade_display = legacy if legacy not in {"", "N/A", "unassigned"} else "â€”"

    return {
        "clade_display": clade_display,
        "qc_status": qc_status,
    }


def main():
    parser = argparse.ArgumentParser(description="Run Nextclade per sample based on BLAST typing")
    parser.add_argument("--blast-summary", required=True)
    parser.add_argument("--datasets-root", required=True)
    parser.add_argument("--max-days", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("segment4_files", nargs="+")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    per_sample_dir = output_dir / "per_sample"
    per_sample_dir.mkdir(parents=True, exist_ok=True)
    datasets_root = Path(args.datasets_root)
    datasets_root.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.write_text("", encoding="utf-8")

    segment4_map = {}
    for item in args.segment4_files:
        file_path = Path(item)
        if file_path.name.endswith("_segment_4.fasta"):
            sample_id = file_path.name[: -len("_segment_4.fasta")]
            segment4_map[sample_id] = file_path

    rows = []
    with Path(args.blast_summary).open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)

    summary_path = output_dir / "nextclade_summary.tsv"
    with summary_path.open("w", encoding="utf-8") as summary:
        summary.write("sample\ttype_blast\tsubtype_HA\tsubtype_NA\tnextclade_dataset\tclade_display\tqc_status\n")

        for row in rows:
            sample = row["sample"]
            flu_type = row["type_blast"]
            subtype_ha = row["subtype_HA"]
            subtype_na = row["subtype_NA"]
            dataset_key = get_nextclade_key(flu_type, subtype_ha)
            sample_out = per_sample_dir / f"{sample}_nextclade.tsv"

            if not dataset_key:
                sample_out.write_text(
                    "seqName\tclade\tshort-clade\tlegacy-clade\tlegacy-clade-vic\tlegacy-clade-yam\tlineage\tqc.overallStatus\n"
                    f"{sample}\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n",
                    encoding="utf-8",
                )
                summary.write(f"{sample}\t{flu_type}\t{subtype_ha}\t{subtype_na}\tN/A\tâ€”\tN/A\n")
                continue

            segment4 = segment4_map.get(sample)
            if not segment4 or not segment4.exists() or segment4.stat().st_size == 0:
                sample_out.write_text(
                    "seqName\tclade\tshort-clade\tlegacy-clade\tlegacy-clade-vic\tlegacy-clade-yam\tlineage\tqc.overallStatus\n"
                    f"{sample}\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\n",
                    encoding="utf-8",
                )
                summary.write(f"{sample}\t{flu_type}\t{subtype_ha}\t{subtype_na}\t{dataset_key}\tâ€”\tN/A\n")
                continue

            dataset_dir = ensure_dataset(dataset_key, datasets_root, args.max_days, log_path)
            command = [
                "nextclade",
                "run",
                "--input-dataset",
                str(dataset_dir),
                "--output-tsv",
                str(sample_out),
                str(segment4),
            ]
            result = subprocess.run(command, capture_output=True, text=True)
            with log_path.open("a", encoding="utf-8") as handle:
                if result.stdout:
                    handle.write(result.stdout)
                    if not result.stdout.endswith("\n"):
                        handle.write("\n")
                if result.stderr:
                    handle.write(result.stderr)
                    if not result.stderr.endswith("\n"):
                        handle.write("\n")

            if result.returncode != 0:
                sample_out.write_text(
                    "seqName\tclade\tshort-clade\tlegacy-clade\tlegacy-clade-vic\tlegacy-clade-yam\tlineage\tqc.overallStatus\n"
                    f"{sample}\tERRO\tERRO\tERRO\tERRO\tERRO\tERRO\tERRO\n",
                    encoding="utf-8",
                )
                summary.write(f"{sample}\t{flu_type}\t{subtype_ha}\t{subtype_na}\t{dataset_key}\tâ€”\tERRO\n")
                continue

            parsed = parse_nextclade_summary(sample_out, flu_type)
            summary.write(
                f"{sample}\t{flu_type}\t{subtype_ha}\t{subtype_na}\t{dataset_key}\t{parsed['clade_display']}\t{parsed['qc_status']}\n"
            )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
