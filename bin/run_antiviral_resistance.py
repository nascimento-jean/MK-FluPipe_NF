#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


HEADER = [
    "sample",
    "gene",
    "aa_position",
    "wt_who",
    "mut_who",
    "alt_observed",
    "frequency",
    "depth_total",
    "drug",
    "significance",
    "nomenclature",
]


def truthy(text: str) -> bool:
    return str(text).strip().upper() in {"TRUE", "T", "YES", "Y", "1", "PASS"}


def load_markers(db_path: Path):
    markers = []
    with db_path.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(
            (line for line in handle if line.strip() and not line.startswith("#")),
            delimiter="\t",
        ):
            markers.append(
                {
                    "gene": row["gene"].strip(),
                    "subtipo": row["subtipo"].strip(),
                    "codon_aa": int(row["codon_aa"]),
                    "wt_aa": row["wt_aa"].strip(),
                    "mut_aa": row["mut_aa"].strip(),
                    "drug": row["drug"].strip(),
                    "significance": row["significance"].strip(),
                }
            )
    return markers


def load_blast_summary(path: Path):
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(row)
    return rows


def load_fasta_sequence(path: Path):
    seq_parts = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line.upper())
    return "".join(seq_parts)


def infer_gene(tsv_path: Path):
    match = re.match(r"(.+)_(NA|PA|MP)_canonical_ivar\.tsv$", tsv_path.name)
    return match.group(2) if match else None


def infer_gene_from_medaka_dir(dir_path: Path):
    match = re.match(r"(.+)_(NA|PA|MP)_canonical_medaka$", dir_path.name)
    return match.group(2) if match else None


def get_field(row, names, default=""):
    for name in names:
        if name in row and row[name] not in (None, ""):
            return row[name]
    return default


def parse_ivar_rows(tsv_path: Path, min_freq: float):
    with tsv_path.open("r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pass_flag = get_field(row, ["PASS", "pass"], "")
            if pass_flag and not truthy(pass_flag):
                continue

            pos = get_field(row, ["POS", "Pos", "pos"], "")
            alt_freq = get_field(row, ["ALT_FREQ", "Alt_freq", "alt_freq"], "")
            total_dp = get_field(row, ["TOTAL_DP", "Total_DP", "total_dp"], "")
            alt_aa = get_field(row, ["ALT_AA", "Alt_AA", "alt_aa"], "").strip()

            try:
                pos_int = int(float(pos))
                freq_float = float(alt_freq)
            except ValueError:
                continue

            if freq_float < min_freq:
                continue

            if not alt_aa or alt_aa in {"?", "N", "N/D", "NA"}:
                continue

            yield {
                "pos": pos_int,
                "aa_pos": ((pos_int - 1) // 3) + 1,
                "alt_freq": freq_float,
                "total_dp": total_dp or "0",
                "alt_aa": alt_aa,
            }


GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_codon(codon: str):
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3 or any(base not in "ACGT" for base in codon):
        return ""
    return GENETIC_CODE.get(codon, "")


def parse_info_field(info_text: str):
    data = {}
    for entry in info_text.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            data[key] = value
        elif entry:
            data[entry] = True
    return data


def parse_medaka_vcf(vcf_path: Path, ref_sequence: str):
    with vcf_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            pos = int(fields[1])
            ref = fields[3].upper()
            alts = [alt.upper() for alt in fields[4].split(",") if alt and alt != "."]
            flt = fields[6]
            info = parse_info_field(fields[7])

            if flt not in {"PASS", "."}:
                continue

            depth = info.get("DP") or info.get("DPS") or "0"
            af_values = info.get("AF", "")
            af_list = af_values.split(",") if af_values else []

            for alt_idx, alt in enumerate(alts):
                if len(ref) != 1 or len(alt) != 1:
                    continue
                if pos < 1 or pos > len(ref_sequence):
                    continue

                aa_pos = ((pos - 1) // 3) + 1
                codon_start = ((aa_pos - 1) * 3)
                if codon_start + 3 > len(ref_sequence):
                    continue

                mutated = list(ref_sequence[codon_start:codon_start + 3])
                mutated[(pos - 1) % 3] = alt
                alt_aa = translate_codon("".join(mutated))
                if not alt_aa:
                    continue

                try:
                    alt_freq = float(af_list[alt_idx]) if alt_idx < len(af_list) and af_list[alt_idx] else None
                except ValueError:
                    alt_freq = None

                yield {
                    "pos": pos,
                    "aa_pos": aa_pos,
                    "alt_freq": alt_freq,
                    "total_dp": depth,
                    "alt_aa": alt_aa,
                }


def main():
    parser = argparse.ArgumentParser(description="Cross-reference canonical iVar calls with antiviral resistance markers")
    parser.add_argument("--db", required=True)
    parser.add_argument("--blast-summary", required=True)
    parser.add_argument("--canonical-refs-dir")
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--ivar-freq", type=float, default=0.03)
    parser.add_argument("sample_dirs", nargs="*")
    args = parser.parse_args()

    markers = load_markers(Path(args.db))
    blast_rows = load_blast_summary(Path(args.blast_summary))
    sample_dir_map = {Path(item).name: Path(item) for item in args.sample_dirs}
    canonical_refs_dir = Path(args.canonical_refs_dir) if args.canonical_refs_dir else None

    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    log_lines = []

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(HEADER)

        for row in blast_rows:
            sample = row.get("sample", "").strip()
            flu_type = row.get("type_blast", "").strip()
            subtype_na = row.get("subtype_NA", "").strip()

            if not sample:
                continue
            if flu_type == "B":
                log_lines.append(f"Skipping {sample}: Influenza B has no WHO canonical resistance mapping")
                continue
            if flu_type != "A":
                log_lines.append(f"Skipping {sample}: unsupported type '{flu_type or 'nd'}'")
                continue

            sample_dir = sample_dir_map.get(sample)
            if not sample_dir or not sample_dir.exists():
                log_lines.append(f"Skipping {sample}: no canonical iVar directory found")
                continue

            subtype_key = "N2" if subtype_na.upper().startswith("N2") else "N1"
            log_lines.append(f"Analyzing {sample} with subtype key {subtype_key}")

            findings = 0
            for ivar_tsv in sorted(sample_dir.glob("*_canonical_ivar.tsv")):
                gene = infer_gene(ivar_tsv)
                if not gene:
                    continue

                gene_bank = "M2" if gene == "MP" else gene
                relevant_markers = [
                    marker
                    for marker in markers
                    if marker["gene"] == gene_bank
                    and marker["subtipo"] in {"ALL", subtype_key}
                ]
                if not relevant_markers:
                    continue

                for call in parse_ivar_rows(ivar_tsv, args.ivar_freq):
                    for marker in relevant_markers:
                        if marker["codon_aa"] != call["aa_pos"]:
                            continue
                        if call["alt_aa"] != marker["mut_aa"]:
                            continue

                        writer.writerow(
                            [
                                sample,
                                gene_bank,
                                call["aa_pos"],
                                marker["wt_aa"],
                                marker["mut_aa"],
                                call["alt_aa"],
                                f"{call['alt_freq'] * 100:.2f}%",
                                call["total_dp"],
                                marker["drug"],
                                marker["significance"],
                                f"{marker['wt_aa']}{call['aa_pos']}{call['alt_aa']}",
                            ]
                        )
                        findings += 1
                        written += 1

            if canonical_refs_dir and canonical_refs_dir.exists():
                for medaka_dir in sorted(sample_dir.glob("*_canonical_medaka")):
                    gene = infer_gene_from_medaka_dir(medaka_dir)
                    if not gene:
                        continue

                    gene_bank = "M2" if gene == "MP" else gene
                    relevant_markers = [
                        marker
                        for marker in markers
                        if marker["gene"] == gene_bank
                        and marker["subtipo"] in {"ALL", subtype_key}
                    ]
                    if not relevant_markers:
                        continue

                    ref_key = f"{subtype_key}_{gene if gene != 'MP' else 'MP'}"
                    ref_path = canonical_refs_dir / f"{ref_key}.fa"
                    if not ref_path.exists():
                        continue
                    ref_sequence = load_fasta_sequence(ref_path)
                    if not ref_sequence:
                        continue

                    vcf_path = medaka_dir / "medaka.annotated.vcf"
                    if not vcf_path.exists():
                        continue

                    for call in parse_medaka_vcf(vcf_path, ref_sequence):
                        for marker in relevant_markers:
                            if marker["codon_aa"] != call["aa_pos"]:
                                continue
                            if call["alt_aa"] != marker["mut_aa"]:
                                continue

                            freq_value = (
                                f"{call['alt_freq'] * 100:.2f}%"
                                if call["alt_freq"] is not None
                                else "NA"
                            )
                            writer.writerow(
                                [
                                    sample,
                                    gene_bank,
                                    call["aa_pos"],
                                    marker["wt_aa"],
                                    marker["mut_aa"],
                                    call["alt_aa"],
                                    freq_value,
                                    call["total_dp"],
                                    marker["drug"],
                                    marker["significance"],
                                    f"{marker['wt_aa']}{call['aa_pos']}{call['alt_aa']}",
                                ]
                            )
                            findings += 1
                            written += 1

            if findings == 0:
                log_lines.append(f"  -> no resistance mutations detected at monitored positions")
            else:
                log_lines.append(f"  -> {findings} resistance marker(s) found")

    log_lines.append(f"Total rows written: {written}")
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
