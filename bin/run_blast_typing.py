#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from pathlib import Path


BLAST_OUTFMT = "6 qseqid sseqid pident length evalue bitscore stitle"


def parse_blast_hit(stitle: str):
    stitle_lower = stitle.lower()
    if re.search(r"influenza a|type a|\(a/", stitle_lower):
        flu_type = "A"
    elif re.search(r"influenza b|type b|\(b/", stitle_lower):
        flu_type = "B"
    elif re.search(r"influenza c|\(c/", stitle_lower):
        flu_type = "C"
    elif re.search(r"influenza d|\(d/", stitle_lower):
        flu_type = "D"
    else:
        flu_type = "Not determined"

    subtype_match = re.search(r"H\d{1,2}N\d{1,2}", stitle, re.IGNORECASE)
    subtype = subtype_match.group(0).upper() if subtype_match else ""

    if flu_type == "B" and not subtype:
        if re.search(r"victoria|\bvic\b", stitle, re.IGNORECASE):
            subtype = "Victoria"
        elif re.search(r"yamagata|\byam\b", stitle, re.IGNORECASE):
            subtype = "Yamagata"
        else:
            subtype = "nd"

    return flu_type, subtype or "nd"


def parse_blast_results(blast_output: Path):
    if not blast_output.exists() or blast_output.stat().st_size == 0:
        return []

    results = []
    for line in blast_output.read_text(encoding="utf-8", errors="ignore").splitlines():
        fields = line.split("\t")
        if len(fields) < 7:
            continue
        try:
            results.append(
                {
                    "qseqid": fields[0],
                    "sseqid": fields[1],
                    "pident": float(fields[2]),
                    "length": int(float(fields[3])),
                    "evalue": float(fields[4]),
                    "bitscore": float(fields[5]),
                    "stitle": fields[6],
                }
            )
        except ValueError:
            continue
    return results


def select_best_hit(results):
    if not results:
        return None
    return sorted(
        results,
        key=lambda item: (
            -item["pident"],
            -item["bitscore"],
            -item["length"],
            item["evalue"],
        ),
    )[0]


def run_blast(query: Path, blast_db_prefix: str, output_path: Path, threads: int, log_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if not query.exists() or query.stat().st_size == 0:
        output_path.write_text("", encoding="utf-8")
        return None

    command = [
        "blastn",
        "-query",
        str(query),
        "-db",
        blast_db_prefix,
        "-out",
        str(output_path),
        "-outfmt",
        BLAST_OUTFMT,
        "-max_target_seqs",
        "5",
        "-max_hsps",
        "1",
        "-perc_identity",
        "80",
        "-evalue",
        "1e-10",
        "-num_threads",
        str(threads),
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stdout:
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(result.stdout)
            if not result.stdout.endswith("\n"):
                handle.write("\n")
    if result.stderr:
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(result.stderr)
            if not result.stderr.endswith("\n"):
                handle.write("\n")
    if result.returncode != 0:
        raise RuntimeError(
            f"blastn failed for query '{query.name}' against '{blast_db_prefix}' with exit code {result.returncode}"
        )
    return select_best_hit(parse_blast_results(output_path))


def sample_from_filename(path: Path):
    name = path.name
    match = re.match(r"(.+)_segment_([1-8])\.fasta$", name)
    if not match:
        return None, None
    return match.group(1), match.group(2)


def main():
    parser = argparse.ArgumentParser(description="Run BLAST typing from extracted influenza segments")
    parser.add_argument("--blast-db-prefix", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("segment_files", nargs="+")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("", encoding="utf-8")

    sample_segments = {}
    for item in args.segment_files:
        path = Path(item)
        sample_id, segment = sample_from_filename(path)
        if not sample_id or not segment:
            continue
        sample_segments.setdefault(sample_id, {})[segment] = path

    summary_path = output_dir / "blast_typing_summary.tsv"
    with summary_path.open("w", encoding="utf-8") as summary:
        summary.write("sample\ttype_blast\tsubtype_HA\tsubtype_NA\thit_HA\thit_NA\n")

        for sample_id in sorted(sample_segments):
            seg4 = sample_segments[sample_id].get("4")
            seg6 = sample_segments[sample_id].get("6")

            ha_output = output_dir / f"{sample_id}_HA_blast.tsv"
            na_output = output_dir / f"{sample_id}_NA_blast.tsv"

            ha_result = run_blast(seg4, args.blast_db_prefix, ha_output, args.threads, log_path) if seg4 else None
            na_result = run_blast(seg6, args.blast_db_prefix, na_output, args.threads, log_path) if seg6 else None

            flu_type = "Not determined"
            subtype_ha = "nd"
            hit_ha = "sem_resultado"
            if ha_result:
                flu_type, subtype_full = parse_blast_hit(ha_result["stitle"])
                subtype_ha_match = re.search(r"H\d{1,2}", subtype_full, re.IGNORECASE)
                subtype_ha = subtype_ha_match.group(0).upper() if subtype_ha_match else subtype_full
                hit_ha = ha_result["sseqid"]

            subtype_na = "nd"
            hit_na = "sem_resultado"
            if na_result:
                _, subtype_full_na = parse_blast_hit(na_result["stitle"])
                subtype_na_match = re.search(r"N\d{1,2}", subtype_full_na, re.IGNORECASE)
                subtype_na = subtype_na_match.group(0).upper() if subtype_na_match else subtype_full_na
                hit_na = na_result["sseqid"]

            summary.write(
                f"{sample_id}\t{flu_type}\t{subtype_ha}\t{subtype_na}\t{hit_ha}\t{hit_na}\n"
            )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
