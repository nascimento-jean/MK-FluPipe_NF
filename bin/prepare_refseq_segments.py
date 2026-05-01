#!/usr/bin/env python3

import argparse
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


REFSEQ_CACHE_PREFIXES_TO_PURGE = ("CY", "JN", "AY", "KF")

GENE_TO_SEGMENT = {
    "PB2": (1, "PB2"),
    "PB1": (2, "PB1"),
    "PB1-F2": (2, "PB1"),
    "PA": (3, "PA"),
    "PA-X": (3, "PA"),
    "HA": (4, "HA"),
    "NP": (5, "NP"),
    "NA": (6, "NA"),
    "NB": (6, "NA"),
    "M1": (7, "MP"),
    "M2": (7, "MP"),
    "BM2": (7, "MP"),
    "MP": (7, "MP"),
    "NS1": (8, "NS"),
    "NS2": (8, "NS"),
    "NEP": (8, "NS"),
}


def write_log(log_path: Path, text: str):
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(text)
        if text and not text.endswith("\n"):
            handle.write("\n")


def get_refseq_accs(flu_type: str, subtype_ha: str):
    subtype = (subtype_ha or "").upper()
    if flu_type == "A":
        if subtype.startswith("H1"):
            return ["NC_026438", "NC_026437", "NC_026436", "NC_026433", "NC_026435", "NC_026434", "NC_026432", "NC_026431"]
        if subtype.startswith("H3"):
            return ["NC_007373", "NC_007372", "NC_007371", "NC_007366", "NC_007370", "NC_007369", "NC_007368", "NC_007367"]
        if subtype.startswith("H5"):
            return ["NC_007364", "NC_007363", "NC_007362", "NC_007357", "NC_007361", "NC_007360", "NC_007359", "NC_007358"]
        return ["NC_007373", "NC_007372", "NC_007371", "NC_007366", "NC_007370", "NC_007369", "NC_007368", "NC_007367"]
    if flu_type == "B":
        return ["NC_002205", "NC_002204", "NC_002206", "NC_002207", "NC_002208", "NC_002209", "NC_002210", "NC_002211"]
    return []


def parse_requested_profiles(blast_summary: Path):
    requested = set()
    if not blast_summary.exists():
        return requested

    with blast_summary.open("r", encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: pos for pos, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            flu_type = fields[idx["type_blast"]] if "type_blast" in idx and idx["type_blast"] < len(fields) else ""
            subtype_ha = fields[idx["subtype_HA"]] if "subtype_HA" in idx and idx["subtype_HA"] < len(fields) else ""
            if not flu_type or flu_type == "Not determined":
                continue
            requested.add((flu_type, subtype_ha))
    return requested


def purge_old_cache(output_dir: Path, log_path: Path):
    removed = 0
    for prefix in REFSEQ_CACHE_PREFIXES_TO_PURGE:
        for item in output_dir.glob(f"{prefix}*"):
            if item.is_file():
                item.unlink(missing_ok=True)
                removed += 1
    if removed:
        write_log(log_path, f"Purged {removed} old non-NC_* cached file(s)")


def fetch_text(url: str):
    with urllib.request.urlopen(url, timeout=60) as response:
        return response.read().decode("utf-8", errors="ignore")


def fetch_refseq_fasta(accession: str):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        f"efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    )
    return fetch_text(url)


def fetch_refseq_gff(accession: str):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        f"efetch.fcgi?db=nuccore&id={accession}&rettype=gff3&retmode=text"
    )
    return fetch_text(url)


def ensure_refseq(accession: str, output_dir: Path, log_path: Path):
    fasta_path = output_dir / f"{accession}.fa"
    gff_path = output_dir / f"{accession}.gff3"

    if not fasta_path.exists() or fasta_path.stat().st_size == 0:
        write_log(log_path, f"Downloading RefSeq FASTA: {accession}")
        try:
            fasta_path.write_text(fetch_refseq_fasta(accession), encoding="utf-8")
        except Exception as exc:
            fasta_path.unlink(missing_ok=True)
            write_log(log_path, f"  -> failed to download FASTA for {accession}: {exc}")
            return None, None

    if not gff_path.exists() or gff_path.stat().st_size == 0:
        write_log(log_path, f"Downloading RefSeq GFF3: {accession}")
        try:
            gff_path.write_text(fetch_refseq_gff(accession), encoding="utf-8")
        except Exception as exc:
            gff_path.unlink(missing_ok=True)
            write_log(log_path, f"  -> failed to download GFF3 for {accession}: {exc}")
            return fasta_path, None

    if not gff_path.exists() or gff_path.stat().st_size == 0:
        return fasta_path, None

    gff_text = gff_path.read_text(encoding="utf-8", errors="ignore")
    if "\tCDS\t" not in gff_text:
        write_log(log_path, f"  -> GFF3 for {accession} has no CDS feature; removing for re-download")
        gff_path.unlink(missing_ok=True)
        return fasta_path, None

    return fasta_path, gff_path


def segment_from_gff3(gff_path: Path):
    if not gff_path.exists() or gff_path.stat().st_size == 0:
        return None

    for line in gff_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        attrs = {}
        for item in fields[8].split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attrs[key] = value
        gene = attrs.get("gene", "").replace(" ", "").upper()
        if gene in GENE_TO_SEGMENT:
            return GENE_TO_SEGMENT[gene]
    return None


def run_command(command: list[str], log_path: Path):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stdout:
        write_log(log_path, result.stdout)
    if result.stderr:
        write_log(log_path, result.stderr)
    return result.returncode == 0


def build_indexes(accession: str, fasta_path: Path, seq_mode: str, threads: int, log_path: Path):
    if seq_mode != "long" and shutil.which("bowtie2-build"):
        bt2_candidates = [fasta_path.with_suffix(fasta_path.suffix + ".1.bt2"), fasta_path.with_suffix(fasta_path.suffix + ".1.bt2l")]
        if not any(path.exists() for path in bt2_candidates):
            write_log(log_path, f"Building Bowtie2 index: {accession}")
            run_command(["bowtie2-build", "--threads", str(max(1, threads)), str(fasta_path), str(fasta_path)], log_path)

    if seq_mode == "long" and shutil.which("minimap2"):
        mmi_path = fasta_path.with_suffix(fasta_path.suffix + ".mmi")
        if not mmi_path.exists():
            write_log(log_path, f"Building minimap2 index: {accession}")
            run_command(["minimap2", "-d", str(mmi_path), str(fasta_path)], log_path)


def write_manifest(rows: list[tuple[str, int, str]], output_dir: Path):
    manifest = output_dir / "refseq_manifest.tsv"
    with manifest.open("w", encoding="utf-8") as handle:
        handle.write("accession\tsegment\tsegment_name\n")
        for accession, segment, segment_name in rows:
            handle.write(f"{accession}\t{segment}\t{segment_name}\n")


def main():
    parser = argparse.ArgumentParser(description="Prepare cached RefSeq FASTA/GFF3 assets for Step 10b")
    parser.add_argument("--blast-summary", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seq-mode", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--log-file", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.write_text("", encoding="utf-8")

    purge_old_cache(output_dir, log_path)

    requested_profiles = parse_requested_profiles(Path(args.blast_summary))
    if not requested_profiles:
        requested_profiles = {("A", "H3")}

    accessions = []
    seen = set()
    for flu_type, subtype_ha in sorted(requested_profiles):
        for accession in get_refseq_accs(flu_type, subtype_ha):
            if accession not in seen:
                accessions.append(accession)
                seen.add(accession)

    manifest_rows = []
    for accession in accessions:
        fasta_path, gff_path = ensure_refseq(accession, output_dir, log_path)
        if not fasta_path or not fasta_path.exists():
            continue
        if gff_path and gff_path.exists():
            segment_info = segment_from_gff3(gff_path)
            if segment_info:
                manifest_rows.append((accession, segment_info[0], segment_info[1]))
        build_indexes(accession, fasta_path, args.seq_mode, args.threads, log_path)

    write_manifest(manifest_rows, output_dir)
    write_log(log_path, f"Prepared RefSeq accessions: {len(accessions)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
