#!/usr/bin/env python3
"""Optional SnpEff annotation for MK-FluPipe full variant call outputs."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


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

SEGMENT_NAMES = {str(num): name for _, (num, name) in GENE_TO_SEGMENT.items()}

FIELDS = [
    "sample",
    "type",
    "subtype",
    "segment",
    "segment_name",
    "accession",
    "pos",
    "ref",
    "alt",
    "effect",
    "impact",
    "gene",
    "feature",
    "hgvs_c",
    "hgvs_p",
    "frequency",
    "depth",
    "status",
    "message",
]


def write_log(log_path: Path, text: str) -> None:
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(text)
        if text and not text.endswith("\n"):
            handle.write("\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_blast_summary(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    return {row.get("sample", ""): row for row in rows if row.get("sample")}


def subtype_from_blast(row: dict[str, str]) -> str:
    flu_type = row.get("type_blast", "")
    if flu_type == "B":
        return "B"
    subtype = f"{row.get('subtype_HA', '')}{row.get('subtype_NA', '')}".replace("-", "").replace("nd", "")
    return subtype or "-"


def get_refseq_accs(flu_type: str, subtype_ha: str) -> list[str]:
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


def segment_from_gff3(path: Path) -> tuple[str, str] | None:
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
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
            seg_num, seg_name = GENE_TO_SEGMENT[gene]
            return str(seg_num), seg_name
    return None


def build_refseq_manifest(refseq_dir: Path) -> dict[str, tuple[str, str, Path, Path]]:
    manifest = {}
    for gff in refseq_dir.glob("*.gff3"):
        accession = gff.stem
        fasta = refseq_dir / f"{accession}.fa"
        if not fasta.exists():
            continue
        seg = segment_from_gff3(gff)
        if seg:
            manifest[accession] = (seg[0], seg[1], fasta, gff)
    return manifest


def choose_refseq(
    refs: dict[str, tuple[str, str, Path, Path]],
    flu_type: str,
    subtype_ha: str,
    segment: str,
) -> tuple[str, str, Path, Path] | None:
    for accession in get_refseq_accs(flu_type, subtype_ha):
        ref = refs.get(accession)
        if ref and ref[0] == segment:
            return accession, ref[1], ref[2], ref[3]
    for accession, ref in sorted(refs.items()):
        if ref[0] == segment:
            return accession, ref[1], ref[2], ref[3]
    return None


def parse_ivar_name(path: Path) -> tuple[str, str, str]:
    base = path.name
    if base.endswith("_ivar.tsv"):
        base = base[: -len("_ivar.tsv")]
    parts = base.rsplit("_", 2)
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    return "", "", ""


def first_fasta_id(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                return line[1:].strip().split()[0]
    return path.stem


def normalize_bool(value: str) -> bool:
    return str(value).strip().upper() in {"TRUE", "T", "1", "YES", "PASS"}


def get_value(row: dict[str, str], *names: str) -> str:
    for name in names:
        if name in row and row[name] != "":
            return row[name]
    return ""


def ivar_to_vcf(ivar_tsv: Path, vcf_path: Path, chrom: str, sample: str) -> list[dict[str, str]]:
    rows = read_tsv(ivar_tsv)
    variants = []
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    with vcf_path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency from iVar">\n')
        handle.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth from iVar">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for idx, row in enumerate(rows, start=1):
            pass_value = get_value(row, "PASS")
            if pass_value and not normalize_bool(pass_value):
                continue
            pos = get_value(row, "POS", "Position", "position")
            ref = get_value(row, "REF", "Ref", "ref")
            alt = get_value(row, "ALT", "Alt", "alt")
            freq = get_value(row, "ALT_FREQ", "frequency", "freq") or "."
            depth = get_value(row, "TOTAL_DP", "depth", "DP") or "."
            if not pos or not ref or not alt:
                continue
            variant_id = f"{sample}_{idx}"
            info = f"AF={freq};DP={depth}"
            handle.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t.\tPASS\t{info}\n")
            variants.append({"pos": pos, "ref": ref, "alt": alt, "frequency": freq, "depth": depth})
    return variants


def run_command(command: list[str], log_path: Path) -> tuple[int, str]:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stdout:
        write_log(log_path, result.stdout)
    if result.stderr:
        write_log(log_path, result.stderr)
    return result.returncode, result.stderr.strip() or result.stdout.strip()


def prepare_snpeff_db(accession: str, fasta: Path, gff: Path, work_dir: Path, log_path: Path) -> tuple[str, Path] | None:
    genome = f"mkflupipe_{accession}"
    work_dir = work_dir.resolve()
    config = work_dir / "snpEff.config"
    data_dir = (work_dir / "data").resolve()
    genome_dir = data_dir / genome
    genome_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(fasta.resolve(), genome_dir / "sequences.fa")
    shutil.copy2(gff.resolve(), genome_dir / "genes.gff")
    config.write_text(f"data.dir = {data_dir}\n{genome}.genome : {accession}\n", encoding="utf-8")
    code, message = run_command(["snpEff", "build", "-gff3", "-noCheckCds", "-noCheckProtein", "-c", str(config), genome], log_path)
    if code != 0:
        write_log(log_path, f"SnpEff build failed for {accession}: {message}")
        return None
    return genome, config


def parse_info(info: str) -> dict[str, str]:
    parsed = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
    return parsed


def parse_ann(info: str) -> tuple[str, str, str, str, str, str]:
    parsed = parse_info(info)
    ann = parsed.get("ANN", "")
    if not ann:
        return "-", "-", "-", "-", "-", "-"
    first = ann.split(",", 1)[0]
    parts = first.split("|")
    effect = parts[1] if len(parts) > 1 and parts[1] else "-"
    impact = parts[2] if len(parts) > 2 and parts[2] else "-"
    gene = parts[3] if len(parts) > 3 and parts[3] else "-"
    feature = parts[6] if len(parts) > 6 and parts[6] else "-"
    hgvs_c = parts[9] if len(parts) > 9 and parts[9] else "-"
    hgvs_p = parts[10] if len(parts) > 10 and parts[10] else "-"
    return effect, impact, gene, feature, hgvs_c, hgvs_p


def parse_annotated_vcf(path: Path, sample: str, flu_type: str, subtype: str, segment: str, segment_name: str, accession: str) -> list[dict[str, str]]:
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            effect, impact, gene, feature, hgvs_c, hgvs_p = parse_ann(fields[7])
            rows.append(
                {
                    "sample": sample,
                    "type": flu_type,
                    "subtype": subtype,
                    "segment": segment,
                    "segment_name": segment_name,
                    "accession": accession,
                    "pos": fields[1],
                    "ref": fields[3],
                    "alt": fields[4],
                    "effect": effect,
                    "impact": impact,
                    "gene": gene,
                    "feature": feature,
                    "hgvs_c": hgvs_c,
                    "hgvs_p": hgvs_p,
                    "frequency": info.get("AF", "-"),
                    "depth": info.get("DP", "-"),
                    "status": "annotated",
                    "message": "-",
                }
            )
    return rows


def fallback_rows(variants: list[dict[str, str]], sample: str, flu_type: str, subtype: str, segment: str, segment_name: str, accession: str, status: str, message: str) -> list[dict[str, str]]:
    rows = []
    for var in variants:
        rows.append(
            {
                "sample": sample,
                "type": flu_type,
                "subtype": subtype,
                "segment": segment,
                "segment_name": segment_name,
                "accession": accession,
                "pos": var.get("pos", "-"),
                "ref": var.get("ref", "-"),
                "alt": var.get("alt", "-"),
                "effect": "-",
                "impact": "-",
                "gene": "-",
                "feature": "-",
                "hgvs_c": "-",
                "hgvs_p": "-",
                "frequency": var.get("frequency", "-"),
                "depth": var.get("depth", "-"),
                "status": status,
                "message": message,
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Run optional SnpEff annotation for full variant call outputs")
    parser.add_argument("--sample-dirs", nargs="+", required=True)
    parser.add_argument("--refseq-dir", required=True)
    parser.add_argument("--blast-summary", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("", encoding="utf-8")

    blast = parse_blast_summary(Path(args.blast_summary))
    refs = build_refseq_manifest(Path(args.refseq_dir))
    rows: list[dict[str, str]] = []
    snpeff_available = shutil.which("snpEff") is not None
    if not snpeff_available:
        write_log(log_path, "snpEff executable was not found; writing fallback rows with status=tool_missing")

    work_dir = output_dir / "_snpeff_work"
    for sample_dir_text in args.sample_dirs:
        sample_dir = Path(sample_dir_text)
        for ivar_tsv in sample_dir.glob("*_ivar.tsv"):
            sample, flu_type, segment = parse_ivar_name(ivar_tsv)
            if not sample or not segment:
                continue
            blast_row = blast.get(sample, {})
            subtype = subtype_from_blast(blast_row) if blast_row else "-"
            flu_type_for_ref = blast_row.get("type_blast", flu_type) if blast_row else flu_type
            subtype_ha_for_ref = blast_row.get("subtype_HA", "") if blast_row else ""
            selected_ref = choose_refseq(refs, flu_type_for_ref, subtype_ha_for_ref, segment)
            if selected_ref is None:
                write_log(log_path, f"No RefSeq FASTA/GFF3 found for segment {segment}; skipping {ivar_tsv.name}")
                continue
            accession, segment_name, fasta, gff = selected_ref
            chrom = first_fasta_id(fasta)
            vcf = output_dir / "vcf" / f"{ivar_tsv.stem}.vcf"
            variants = ivar_to_vcf(ivar_tsv, vcf, chrom, sample)
            if not variants:
                continue
            if not snpeff_available:
                rows.extend(fallback_rows(variants, sample, flu_type, subtype, segment, segment_name, accession, "tool_missing", "snpEff executable was not found"))
                continue
            prepared = prepare_snpeff_db(accession, fasta, gff, work_dir / accession, log_path)
            if not prepared:
                rows.extend(fallback_rows(variants, sample, flu_type, subtype, segment, segment_name, accession, "build_failed", "SnpEff database build failed"))
                continue
            genome, config = prepared
            annotated_vcf = output_dir / "vcf" / f"{ivar_tsv.stem}.snpeff.vcf"
            with annotated_vcf.open("w", encoding="utf-8") as handle:
                result = subprocess.run(["snpEff", "ann", "-noStats", "-c", str(config), genome, str(vcf)], stdout=handle, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                write_log(log_path, result.stderr)
                rows.extend(fallback_rows(variants, sample, flu_type, subtype, segment, segment_name, accession, "annotation_failed", result.stderr.strip()))
                continue
            rows.extend(parse_annotated_vcf(annotated_vcf, sample, flu_type, subtype, segment, segment_name, accession))

    write_tsv(output_dir / "snpeff_annotation.tsv", rows)
    write_log(log_path, f"SnpEff annotation rows: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
