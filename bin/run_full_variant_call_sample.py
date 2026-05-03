#!/usr/bin/env python3

import argparse
import csv
import shlex
import subprocess
import sys
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


def run_shell(command: str, log_path: Path, allow_fail: bool = False):
    result = subprocess.run(
        command,
        shell=True,
        executable="/usr/bin/bash",
        capture_output=True,
        text=True,
    )
    if result.stdout:
        write_log(log_path, result.stdout)
    if result.stderr:
        write_log(log_path, result.stderr)
    if result.returncode != 0 and not allow_fail:
        raise RuntimeError(f"Command failed ({result.returncode}): {command}")
    return result.returncode


def parse_ivar_to_protein(ivar_tsv: Path, seg_label: str, sample_id: str, flu_type: str, output_path: Path):
    with ivar_tsv.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        header_map = {name: idx for idx, name in enumerate(header or [])}
        rows = []
        for row in reader:
            if not row:
                continue
            pass_value = get_field(row, header_map, 13, "PASS")
            if str(pass_value).upper() != "TRUE":
                continue

            alt_freq = get_field(row, header_map, 10, "ALT_FREQ")
            total_dp = get_field(row, header_map, 11, "TOTAL_DP")
            gff_feature = get_field(row, header_map, 14, "GFF_FEATURE")
            ref_aa = get_field(row, header_map, 16, "REF_AA")
            alt_aa = get_field(row, header_map, 18, "ALT_AA")
            pos_aa = get_field(row, header_map, 19, "POS_AA")

            gene = str(gff_feature).split(":", 1)[0] if gff_feature and gff_feature != "NA" else "UTR"
            ref_aa = ref_aa or row[2] if len(row) > 2 else "."
            alt_aa = alt_aa or row[3] if len(row) > 3 else "."
            pos_aa = pos_aa or row[1] if len(row) > 1 else "."

            if not ref_aa or ref_aa == "NA" or not alt_aa or alt_aa == "NA":
                mutation = "."
                gene = "UTR"
            elif ref_aa == alt_aa:
                mutation = "."
            else:
                mutation = f"{ref_aa}{pos_aa}{alt_aa}"

            try:
                freq_text = f"{float(alt_freq) * 100:.2f}%"
            except Exception:
                freq_text = "0.00%"

            rows.append(
                [
                    sample_id,
                    flu_type,
                    seg_label,
                    gene,
                    pos_aa,
                    ref_aa,
                    alt_aa,
                    freq_text,
                    str(total_dp or "0"),
                    mutation,
                ]
            )

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "type", "segment", "gene", "aa_position", "ref_aa", "alt_aa", "frequency", "depth", "mutation"])
        writer.writerows(rows)
    return len(rows)


def get_field(row, header_map, fallback_index, *candidate_names):
    for name in candidate_names:
        if name in header_map and header_map[name] < len(row):
            return row[header_map[name]]
    if fallback_index < len(row):
        return row[fallback_index]
    return ""


def main():
    parser = argparse.ArgumentParser(description="Run Step 10b full variant calling for one sample")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--flu-type", required=True)
    parser.add_argument("--subtype-ha", default="")
    parser.add_argument("--subtype-na", default="")
    parser.add_argument("--seq-mode", required=True)
    parser.add_argument("--refseq-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--ivar-freq", type=float, required=True)
    parser.add_argument("--ivar-depth", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--reads", nargs="+", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.write_text("", encoding="utf-8")
    refseq_dir = Path(args.refseq_dir)

    if not args.flu_type or args.flu_type == "Not determined":
        write_log(log_path, f"Skipping {args.sample_id}: virus type unknown")
        write_manifest(output_dir, args.sample_id, args.flu_type, args.subtype_ha, args.subtype_na, 0, 0)
        return

    align_threads = max(1, args.threads - 1)
    sort_threads = max(1, args.threads - align_threads)
    segments_done = 0
    segments_skipped = 0

    for accession in get_refseq_accs(args.flu_type, args.subtype_ha):
        refseq_fa = refseq_dir / f"{accession}.fa"
        refseq_gff = refseq_dir / f"{accession}.gff3"
        if not refseq_fa.exists() or refseq_fa.stat().st_size == 0:
            write_log(log_path, f"Skipping accession {accession}: FASTA missing")
            segments_skipped += 1
            continue
        if not refseq_gff.exists() or refseq_gff.stat().st_size == 0:
            write_log(log_path, f"Skipping accession {accession}: GFF3 missing")
            segments_skipped += 1
            continue

        seg_info = segment_from_gff3(refseq_gff)
        if not seg_info:
            write_log(log_path, f"Skipping accession {accession}: unable to determine segment from GFF3")
            segments_skipped += 1
            continue

        seg_num, seg_name = seg_info
        seg_label = f"{args.sample_id}_{args.flu_type}_{seg_num}"
        bam_out = output_dir / f"{args.sample_id}_{seg_name}_refseq.bam"
        ivar_prefix = output_dir / f"{seg_label}_ivar"
        protein_tsv = output_dir / f"{seg_label}_protein_mutations.tsv"
        write_log(log_path, f"Processing {accession} -> segment {seg_num} ({seg_name})")

        if args.seq_mode == "long":
            align_cmd = (
                f"minimap2 -ax map-ont -t {align_threads} {shlex.quote(str(refseq_fa))} {shlex.quote(args.reads[0])} "
                f"| samtools sort -@ {sort_threads} -o {shlex.quote(str(bam_out))} -"
            )
        elif len(args.reads) > 1:
            align_cmd = (
                f"bowtie2 -x {shlex.quote(str(refseq_fa))} -1 {shlex.quote(args.reads[0])} -2 {shlex.quote(args.reads[1])} "
                f"-p {align_threads} --very-sensitive | samtools sort -@ {sort_threads} -o {shlex.quote(str(bam_out))} -"
            )
        else:
            align_cmd = (
                f"bowtie2 -x {shlex.quote(str(refseq_fa))} -U {shlex.quote(args.reads[0])} "
                f"-p {align_threads} --very-sensitive | samtools sort -@ {sort_threads} -o {shlex.quote(str(bam_out))} -"
            )

        run_shell(align_cmd, log_path, allow_fail=True)
        if not bam_out.exists() or bam_out.stat().st_size == 0:
            write_log(log_path, f"  -> empty BAM for {accession}")
            segments_skipped += 1
            continue

        run_shell(f"samtools index {shlex.quote(str(bam_out))}", log_path, allow_fail=True)

        mpileup_base = (
            f"samtools mpileup -A -d 0 -B -Q 0 --reference {shlex.quote(str(refseq_fa))} {shlex.quote(str(bam_out))}"
        )
        ivar_base = (
            f"ivar variants -p {shlex.quote(str(ivar_prefix))} -q 20 -t {args.ivar_freq} -m {args.ivar_depth} -r {shlex.quote(str(refseq_fa))}"
        )
        with_gff = f"{mpileup_base} | {ivar_base} -g {shlex.quote(str(refseq_gff))}"
        without_gff = f"{mpileup_base} | {ivar_base}"

        run_shell(with_gff, log_path, allow_fail=True)
        ivar_tsv = Path(f"{ivar_prefix}.tsv")
        if not ivar_tsv.exists() or ivar_tsv.stat().st_size == 0:
            write_log(log_path, f"  -> retrying iVar without GFF3 for {accession}")
            run_shell(without_gff, log_path, allow_fail=True)

        if not ivar_tsv.exists() or ivar_tsv.stat().st_size == 0:
            write_log(log_path, f"  -> no iVar TSV for {accession}")
            segments_skipped += 1
            continue

        n_rows = parse_ivar_to_protein(ivar_tsv, seg_label, args.sample_id, args.flu_type, protein_tsv)
        write_log(log_path, f"  -> {ivar_tsv.name}: {n_rows} protein-annotated row(s)")
        segments_done += 1

    write_manifest(output_dir, args.sample_id, args.flu_type, args.subtype_ha, args.subtype_na, segments_done, segments_skipped)


def write_manifest(output_dir: Path, sample_id: str, flu_type: str, subtype_ha: str, subtype_na: str, segments_done: int, segments_skipped: int):
    return


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
