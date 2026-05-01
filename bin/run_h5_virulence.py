#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


HEADER = [
    "sample",
    "gene",
    "aa_position",
    "wt_expected",
    "mut_reference",
    "observed_aa",
    "description",
    "nomenclature",
    "status",
]

H5_MARKERS = {
    ("PB2", 627): ("E", "K", "Mammalian adaptation - enhanced replication at 33C"),
    ("PB2", 701): ("D", "N", "Mammalian adaptation - enhanced replication"),
    ("PA", 97): ("I", "V", "Increased virulence in mammals"),
    ("PB1", 66): ("N", "S", "Increased virulence and pathogenicity"),
    ("NS1", 92): ("D", "E", "Interferon resistance"),
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


def load_blast_summary(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_fasta_sequence(path: Path):
    seq_parts = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line.upper())
    return "".join(seq_parts)


def translate_codon(codon: str):
    if len(codon) != 3 or any(base not in "ACGT" for base in codon):
        return ""
    return GENETIC_CODE.get(codon, "")


def find_gene_fasta(irma_dir: Path, gene: str):
    patterns = {
        "PB2": ["*PB2*.fasta"],
        "PB1": ["*PB1*.fasta"],
        "PA": ["*PA*.fasta"],
        "NS1": ["*NS*.fasta", "*NS1*.fasta"],
    }
    for pattern in patterns.get(gene, [f"*{gene}*.fasta"]):
        matches = sorted(irma_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def main():
    parser = argparse.ArgumentParser(description="Evaluate conditional H5 virulence markers from IRMA outputs")
    parser.add_argument("--blast-summary", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("irma_dirs", nargs="*")
    args = parser.parse_args()

    blast_rows = load_blast_summary(Path(args.blast_summary))
    irma_dir_map = {Path(item).name: Path(item) for item in args.irma_dirs}
    output_path = Path(args.output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    h5_samples = []
    for row in blast_rows:
        sample = (row.get("sample") or "").strip()
        subtype_ha = (row.get("subtype_HA") or "").strip().upper()
        if sample and subtype_ha.startswith("H5"):
            h5_samples.append(sample)

    log_lines = []
    rows_written = 0

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(HEADER)

        if not h5_samples:
            log_lines.append("No H5 samples detected. H5 virulence analysis not required.")
        else:
            log_lines.append(f"H5 samples detected: {' '.join(h5_samples)}")

        for sample in h5_samples:
            irma_dir = irma_dir_map.get(sample)
            if not irma_dir or not irma_dir.exists():
                log_lines.append(f"Skipping {sample}: IRMA directory not found")
                continue

            log_lines.append(f"Analyzing H5 virulence markers: {sample}")
            for (gene, aa_pos), (wt, mut, description) in H5_MARKERS.items():
                fasta_path = find_gene_fasta(irma_dir, gene)
                if not fasta_path:
                    continue

                seq = load_fasta_sequence(fasta_path)
                start = (aa_pos - 1) * 3
                codon = seq[start:start + 3]
                obs_aa = translate_codon(codon)
                if not obs_aa:
                    continue

                if obs_aa == wt:
                    status = "WT"
                    nomenclature = f"{wt}{aa_pos}{wt}"
                elif obs_aa == mut:
                    status = "MUT_DETECTED"
                    nomenclature = f"{wt}{aa_pos}{mut}"
                else:
                    status = "VARIANT"
                    nomenclature = f"{wt}{aa_pos}{obs_aa}"

                writer.writerow([
                    sample,
                    gene,
                    aa_pos,
                    wt,
                    mut,
                    obs_aa,
                    description,
                    nomenclature,
                    status,
                ])
                rows_written += 1

            log_lines.append(f"  -> H5 virulence: {sample} analyzed")

    log_lines.append(f"Total rows written: {rows_written}")
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
