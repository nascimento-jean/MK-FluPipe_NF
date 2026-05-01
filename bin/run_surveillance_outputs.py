#!/usr/bin/env python3

import argparse
import csv
import html
import json
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path


GENE_TO_SEGMENT = {
    "PB2": "1",
    "PB1": "2",
    "PA": "3",
    "HA": "4",
    "NP": "5",
    "NA": "6",
    "MP": "7",
    "NS": "8",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate final surveillance and GISAID outputs")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--irma-module", required=True)
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--gisaid-location", default="")
    parser.add_argument("--gisaid-year", default="")
    parser.add_argument("--dependencies", nargs="+", required=True)
    parser.add_argument("--consensus-fastas", nargs="+", required=True)
    parser.add_argument("--irma-dirs", nargs="+", required=True)
    return parser.parse_args()


def read_tsv(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_fasta(path: Path):
    header = None
    seq = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def normalize_sequence(seq: str):
    return "".join(base if base.upper() in {"A", "C", "T", "G"} else "N" for base in seq.upper())


def infer_segments_from_fasta(consensus_path: Path):
    segments = []
    for header, _ in parse_fasta(consensus_path):
        tail = header.split("_")[-1]
        if tail.isdigit():
            segments.append(tail)
            continue
        for gene, seg in GENE_TO_SEGMENT.items():
            if gene in header:
                segments.append(seg)
                break
    unique = sorted(set(segments), key=lambda x: int(x))
    return "|".join(unique) if unique else "Ã¢â‚¬â€"


def infer_type_gene_segment(segment_name: str):
    parts = segment_name.split("_")
    flu_type = parts[0] if parts else "Ã¢â‚¬â€"
    gene = parts[1] if len(parts) > 1 else "Ã¢â‚¬â€"
    seg_num = GENE_TO_SEGMENT.get(gene, "Ã¢â‚¬â€")
    return flu_type, gene, seg_num


def blast_classification(row):
    flu_type = row.get("type_blast", "")
    subtype_ha = row.get("subtype_HA", "nd")
    subtype_na = row.get("subtype_NA", "nd")
    if flu_type == "B":
        return subtype_ha
    if subtype_ha not in {"", "nd"} and subtype_na not in {"", "nd"}:
        return f"{subtype_ha}{subtype_na}"
    if subtype_ha not in {"", "nd"}:
        return subtype_ha
    if subtype_na not in {"", "nd"}:
        return subtype_na
    return "nd"


def pretty_header(name: str):
    return name.replace("_", " ").title() if name else name


def badge(value: str):
    lv = str(value).strip().lower()
    if lv in {"pass", "good", "ok", "wt"}:
        css = "bg-success"
    elif lv in {"warn", "mediocre", "verify"}:
        css = "bg-warning text-dark"
    elif lv in {"fail", "bad", "mut_detected"}:
        css = "bg-danger"
    elif "blast+nextclade" in lv:
        css = "bg-primary"
    elif lv == "variant":
        return '<span class="badge rounded-pill" style="background:#fd7e14;color:#fff">{}</span>'.format(html.escape(str(value)))
    else:
        css = "bg-secondary"
    return '<span class="badge {} rounded-pill">{}</span>'.format(css, html.escape(str(value)))


def render_cell(header: str, value):
    text = "-" if str(value).strip() == "nd" else str(value)
    header_lower = header.lower()
    if header_lower == "type":
        if text == "A":
            return '<span class="badge rounded-pill" style="background:#66CDAA;color:#1a4a3a">A</span>'
        if text == "B":
            return '<span class="badge rounded-pill" style="background:#FFDEAD;color:#5a3e1b">B</span>'
    if any(key in header_lower for key in ("qc_assembly", "qc_nextclade", "classification_source", "coinfection_status", "status")):
        return badge(text)
    if header_lower == "nomenclature" and text not in {"N/D", "?", ""}:
        clean = text.replace("?", "")
        if len(clean) >= 3 and clean[0] == clean[-1]:
            return '<span style="color:#198754;font-weight:600">{}</span>'.format(html.escape(text))
        return '<span style="color:#dc3545;font-weight:700">{}</span>'.format(html.escape(text))
    if header_lower == "mutation":
        if text in {".", "", "-"}:
            return '<span style="color:#adb5bd;font-size:.78rem">syn</span>'
        if "*" in text:
            return '<span style="color:#dc3545;font-weight:700;background:#ffeaea;padding:1px 4px;border-radius:3px">{}</span>'.format(html.escape(text))
        return '<span style="color:#6f42c1;font-weight:700">{}</span>'.format(html.escape(text))
    if header_lower in {"frequency", "freq", "alt_freq"}:
        try:
            frequency_value = float(text.replace("%", "").strip())
            if frequency_value >= 80:
                color = "#198754"
            elif frequency_value >= 30:
                color = "#fd7e14"
            else:
                color = "#dc3545"
            return '<span style="color:{};font-weight:600">{}</span>'.format(color, html.escape(text))
        except Exception:
            pass
    return html.escape(text)


def render_download_button(table_id: str):
    return (
        '<div class="text-end mb-1">'
        '<button onclick="dlExcel(\'{}\')" class="btn btn-sm btn-outline-success" style="font-size:.75rem">'
        '&#11015; Download Excel</button></div>'
    ).format(table_id)


def render_html_table(table_id: str, rows, headers, caption=None):
    if not headers:
        return '<p class="text-muted fst-italic">No data.</p>'
    head_html = "".join(f"<th>{html.escape(pretty_header(h))}</th>" for h in headers)
    body = []
    for row in rows:
        cells = "".join(f"<td>{render_cell(headers[i], row.get(headers[i], ''))}</td>" for i in range(len(headers)))
        body.append(f"<tr>{cells}</tr>")
    caption_html = ""
    if caption:
        caption_html = '<caption style="caption-side:top;font-weight:600;color:#495057;font-size:.9rem;padding:.4rem 0">{}</caption>'.format(
            html.escape(caption)
        )
    return (
        render_download_button(table_id)
        + '<div class="table-responsive mt-1">'
        + f'<table id="{table_id}" class="table table-sm table-striped table-hover align-middle" style="font-size:.82rem">'
        + caption_html
        + f'<thead class="table-dark"><tr>{head_html}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'
    )


def summarize_protein_mutations(fullvar_rows):
    gene_order = {
        "PB2": 1,
        "PB1": 2,
        "PB1-F2": 3,
        "PA": 2,
        "PA-X": 3,
        "HA": 2,
        "NP": 2,
        "NA": 2,
        "NB": 3,
        "M1": 2,
        "M2": 3,
        "BM2": 3,
        "NS1": 2,
        "NEP": 3,
        "NS2": 3,
    }
    color_map = {
        "nonsynonymous": "#d62728",
        "insertion": "#2ca02c",
        "deletion": "#006400",
    }

    def parse_frequency(value):
        try:
            return float(str(value).replace("%", "").strip())
        except Exception:
            return None

    def mutation_class(gene, alt_aa, mutation):
        alt = str(alt_aa)
        if alt.startswith("+"):
            return "insertion"
        if alt.startswith("-"):
            return "deletion"
        if gene == "UTR":
            return "utr"
        if str(mutation) == ".":
            return "synonymous"
        return "nonsynonymous"

    def mutation_label(gene, aa_pos, ref_aa, alt_aa, mutation, mclass):
        if mclass == "insertion":
            return f"{gene}:{ref_aa}{aa_pos}ins{str(alt_aa).lstrip('+')}"
        if mclass == "deletion":
            return f"{gene}:{ref_aa}{aa_pos}del{str(alt_aa).lstrip('-')}"
        if mclass == "nonsynonymous":
            return f"{gene}:{mutation}"
        return None

    def segment_number(segment_id):
        tail = str(segment_id).split("_")[-1]
        try:
            return int(tail)
        except Exception:
            return 999

    grouped = defaultdict(lambda: defaultdict(list))
    segment_sort = {}
    for row in fullvar_rows:
        frequency = parse_frequency(row.get("frequency", ""))
        try:
            depth = float(row.get("depth", "0"))
        except Exception:
            depth = 0
        gene = str(row.get("gene", ""))
        mutation = row.get("mutation", "")
        alt_aa = row.get("alt_aa", "")
        if frequency is None or depth < 10 or frequency <= 75:
            continue
        mclass = mutation_class(gene, alt_aa, mutation)
        if mclass not in {"nonsynonymous", "insertion", "deletion"}:
            continue
        if gene == "UTR":
            continue
        label = mutation_label(gene, row.get("aa_position", ""), row.get("ref_aa", ""), alt_aa, mutation, mclass)
        if not label:
            continue
        span = '<span style="color:{};font-weight:600">{}</span>'.format(color_map.get(mclass, "#333"), html.escape(label))
        key = (row.get("sample", ""), row.get("type", ""), row.get("segment", ""))
        grouped[key][gene].append(span)
        segment_sort[key] = segment_number(row.get("segment", ""))

    summary_rows = []
    for key in sorted(grouped.keys(), key=lambda item: (item[0], segment_sort.get(item, 999))):
        sample, virus_type, segment = key
        gene_map = grouped[key]
        sorted_genes = sorted(gene_map.keys(), key=lambda gene: (gene_order.get(gene, 50), gene))
        mutation_blocks = []
        for gene in sorted_genes:
            mutation_blocks.append(f"<b>{html.escape(gene)}</b>:[{', '.join(gene_map[gene])}]")
        summary_rows.append(
            {
                "sample": sample,
                "type": virus_type,
                "segment": segment,
                "n_mutations": str(sum(len(values) for values in gene_map.values())),
                "mutations_by_segment": " | ".join(mutation_blocks),
            }
        )
    return summary_rows


def render_protein_table(table_id: str, rows):
    headers = ["sample", "type", "segment", "n_mutations", "mutations_by_segment"]
    if not rows:
        return '<p class="text-muted fst-italic">No data.</p>'
    head_html = "".join(f"<th>{html.escape(pretty_header(h))}</th>" for h in headers)
    table_rows = []
    for row in rows:
        cells = []
        for header in headers:
            value = row.get(header, "")
            if header == "mutations_by_segment":
                cells.append(
                    '<td style="white-space:normal;word-break:break-word;min-width:300px;text-align:left !important">{}</td>'.format(value)
                )
            else:
                cells.append(f"<td>{render_cell(header, value)}</td>")
        table_rows.append("<tr>" + "".join(cells) + "</tr>")
    legend = (
        '<div style="margin:.4rem 0 .6rem;font-size:.8rem">'
        '<b>Legend:</b> '
        '<span style="color:#d62728;font-weight:600">Non-synonymous</span> | '
        '<span style="color:#2ca02c;font-weight:600">Insertion</span> | '
        '<span style="color:#006400;font-weight:600">Deletion</span>'
        ' &nbsp;&mdash;&nbsp; Filters: depth &ge; 10&times;, frequency &gt; 75%'
        '</div>'
    )
    return (
        legend
        + render_download_button(table_id)
        + '<div class="table-responsive mt-1">'
        + f'<table id="{table_id}" class="table table-sm table-striped table-hover align-middle" style="font-size:.82rem">'
        + f'<thead class="table-dark"><tr>{head_html}</tr></thead><tbody>{"".join(table_rows)}</tbody></table></div>'
    )


def build_gisaid_name(summary_row, location: str, year: str):
    sample = summary_row["sample"]
    flu_type = summary_row["type"]
    subtype_ha = summary_row["subtype_HA"]
    subtype_na = summary_row["subtype_NA"]
    blast_class = summary_row["blast_classification"]
    loc = location or "Brasil"
    if flu_type == "A":
        if subtype_ha != "nd" and subtype_na != "nd":
            subtype = f"{subtype_ha}{subtype_na}"
            return f"A/{loc}/{sample}/{year}({subtype})", subtype
        return f"A/{loc}/{sample}/{year}", blast_class or "H?N?"
    if flu_type == "B":
        label = blast_class or "B"
        if "victoria" in label.lower():
            return f"B/{loc}/{sample}/{year}(Victoria)", "B"
        if "yamagata" in label.lower():
            return f"B/{loc}/{sample}/{year}(Yamagata)", "B"
        return f"B/{loc}/{sample}/{year}", "B"
    return f"{flu_type or 'A'}/{loc}/{sample}/{year}", blast_class or "unknown"


def rewrite_gisaid_sequences(summary_rows, consensus_map, out_fasta: Path):
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w", encoding="utf-8") as handle:
        for row in summary_rows:
            sample = row["sample"]
            consensus_path = consensus_map.get(sample)
            if not consensus_path or not consensus_path.exists():
                continue
            for header, seq in parse_fasta(consensus_path):
                seg_num = ""
                tail = header.split("_")[-1]
                if tail.isdigit():
                    seg_num = tail
                else:
                    for gene, seg in GENE_TO_SEGMENT.items():
                        if gene in header:
                            seg_num = seg
                            break
                if not seg_num:
                    continue
                handle.write(f">{sample}_{row['type']}_{seg_num}\n{normalize_sequence(seq)}\n")


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    dep_map = {Path(item).name: Path(item) for item in args.dependencies}
    consensus_paths = [Path(item) for item in args.consensus_fastas]
    irma_dirs = [Path(item) for item in args.irma_dirs]

    consensus_map = {path.stem: path for path in consensus_paths}

    blast_rows = read_tsv(dep_map.get("blast_typing_summary.tsv", Path("__missing__")))
    nextclade_rows = read_tsv(dep_map.get("nextclade_summary.tsv", Path("__missing__")))
    assembly_qc_rows = read_tsv(dep_map.get("assembly_qc_report.tsv", Path("__missing__")))
    depth_rows = read_tsv(dep_map.get("depth_summary.tsv", Path("__missing__")))
    coinfection_rows = read_tsv(dep_map.get("coinfection_report.tsv", Path("__missing__")))
    antiviral_rows = read_tsv(dep_map.get("antiviral_resistance.tsv", Path("__missing__")))
    h5_rows = read_tsv(dep_map.get("h5_virulence_markers.tsv", Path("__missing__")))
    fullvar_rows = read_tsv(dep_map.get("all_samples_protein_mutations.tsv", Path("__missing__")))

    nextclade_map = {row["sample"]: row for row in nextclade_rows}
    assembly_qc_map = {row["sample"]: row for row in assembly_qc_rows}

    typing_rows = []
    for row in blast_rows:
        sample = row["sample"]
        next_row = nextclade_map.get(sample, {})
        qc_row = assembly_qc_map.get(sample, {})
        clade_display = next_row.get("clade_display", "Ã¢â‚¬â€") or "Ã¢â‚¬â€"
        qc_nextclade = next_row.get("qc_status", "Ã¢â‚¬â€") or "Ã¢â‚¬â€"
        source = "BLAST+Nextclade" if clade_display not in {"", "Ã¢â‚¬â€"} else "BLAST"

        typing_rows.append(
            {
                "sample": sample,
                "type": row.get("type_blast", "Ã¢â‚¬â€"),
                "subtype_HA": row.get("subtype_HA", "nd"),
                "subtype_NA": row.get("subtype_NA", "nd"),
                "blast_classification": blast_classification(row),
                "nextclade_clade": clade_display,
                "qc_nextclade": qc_nextclade,
                "classification_source": source,
                "hit_blast_HA": row.get("hit_HA", "sem_resultado"),
                "hit_blast_NA": row.get("hit_NA", "sem_resultado"),
                "segments": infer_segments_from_fasta(consensus_map[sample]) if sample in consensus_map else "Ã¢â‚¬â€",
                "qc_assembly": qc_row.get("qc_assembly", "Ã¢â‚¬â€"),
                "qc_detail": qc_row.get("qc_detail", "Ã¢â‚¬â€"),
            }
        )

    typing_path = out_dir / "typing_results.tsv"
    typing_fields = [
        "sample",
        "type",
        "subtype_HA",
        "subtype_NA",
        "blast_classification",
        "nextclade_clade",
        "qc_nextclade",
        "classification_source",
        "hit_blast_HA",
        "hit_blast_NA",
        "segments",
        "qc_assembly",
        "qc_detail",
    ]
    write_tsv(typing_path, typing_rows, typing_fields)

    coverage_rows = []
    for row in depth_rows:
        flu_type, gene, seg_num = infer_type_gene_segment(row.get("segment", ""))
        coverage_rows.append(
            {
                "sample": row.get("sample", ""),
                "type": flu_type,
                "gene": gene,
                "segment": seg_num,
                "cov_mean": row.get("cov_mean", ""),
                "cov_min": row.get("cov_min", ""),
                "cov_max": row.get("cov_max", ""),
                "positions_covered": row.get("positions_covered", ""),
                "ref_length": row.get("ref_length", ""),
            }
        )
    coverage_path = out_dir / "coverage_per_segment.tsv"
    write_tsv(
        coverage_path,
        coverage_rows,
        ["sample", "type", "gene", "segment", "cov_mean", "cov_min", "cov_max", "positions_covered", "ref_length"],
    )

    run_summary_rows = []
    today = date.today().isoformat()
    for row in typing_rows:
        merged = dict(row)
        merged["analysis_date"] = today
        merged["irma_module"] = args.irma_module
        merged["pipeline_version"] = args.pipeline_version
        run_summary_rows.append(merged)

    run_summary_path = out_dir / "run_summary.tsv"
    write_tsv(
        run_summary_path,
        run_summary_rows,
        typing_fields + ["analysis_date", "irma_module", "pipeline_version"],
    )

    run_summary_json = out_dir / "run_summary.json"
    run_summary_json.write_text(
        json.dumps(
            {
                "pipeline": "MK Flu-Pipe Nextflow",
                "version": args.pipeline_version,
                "analysis_date": today,
                "irma_module": args.irma_module,
                "total_samples": len(run_summary_rows),
                "samples": run_summary_rows,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    multisample_path = out_dir / "multisample_consensus.fasta"
    with multisample_path.open("w", encoding="utf-8") as handle:
        for consensus_path in consensus_paths:
            if consensus_path.exists() and consensus_path.stat().st_size > 0:
                for header, seq in parse_fasta(consensus_path):
                    handle.write(f">{header}\n{normalize_sequence(seq)}\n")

    html_path = out_dir / "surveillance_report.html"
    protein_summary_rows = summarize_protein_mutations(fullvar_rows)
    total_samples = len(typing_rows)
    total_pass = sum(1 for row in typing_rows if row.get("qc_assembly") == "PASS")
    total_warn_fail = sum(1 for row in typing_rows if row.get("qc_assembly") in {"WARN", "FAIL"})
    total_coinfection_warn = sum(1 for row in coinfection_rows if row.get("coinfection_status") == "WARN")
    total_antiviral = len(antiviral_rows)
    total_protein_mutations = sum(1 for row in fullvar_rows if row.get("mutation") not in {".", "", "mutation"})

    css = (
        "body{font-family:'Segoe UI',system-ui,sans-serif;background:#f8f9fa}"
        ".sc{border-left:4px solid;border-radius:.5rem}"
        ".sc.bl{border-color:#0d6efd}.sc.gr{border-color:#198754}"
        ".sc.or{border-color:#fd7e14}.sc.rd{border-color:#dc3545}"
        ".sec{font-size:.95rem;font-weight:700;color:#1e293b;text-transform:uppercase;letter-spacing:.04em;"
        "border-bottom:2px solid #dee2e6;padding-bottom:.3rem;margin:1.4rem 0 .6rem}"
        "thead.table-dark th{font-size:.78rem;white-space:nowrap;text-align:center !important}"
        "td{text-align:center !important}"
        ".badge{font-size:.72rem}"
    )

    cards = (
        '<div class="row g-3 mb-4">'
        '<div class="col-6 col-md-2"><div class="card sc bl p-3"><div class="text-muted" style="font-size:.75rem">Samples</div><div class="fs-3 fw-bold text-primary">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc gr p-3"><div class="text-muted" style="font-size:.75rem">PASS</div><div class="fs-3 fw-bold text-success">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc or p-3"><div class="text-muted" style="font-size:.75rem">WARN / FAIL</div><div class="fs-3 fw-bold text-warning">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc rd p-3"><div class="text-muted" style="font-size:.75rem">Co-infect. (alert)</div><div class="fs-3 fw-bold text-danger">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc or p-3"><div class="text-muted" style="font-size:.75rem">Resist. markers</div><div class="fs-3 fw-bold text-warning">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc bl p-3"><div class="text-muted" style="font-size:.75rem">Protein mutations</div><div class="fs-3 fw-bold" style="color:#6f42c1">{}</div></div></div>'
        '</div>'
        '<div class="row g-3 mb-4"><div class="col-6 col-md-2"><div class="card sc bl p-3"><div class="text-muted" style="font-size:.75rem">Data</div><div style="font-size:.85rem;font-weight:600;color:#0d6efd">{}</div></div></div></div>'
    ).format(
        total_samples,
        total_pass,
        total_warn_fail,
        total_coinfection_warn,
        total_antiviral,
        total_protein_mutations,
        date.today().isoformat(),
    )

    html_parts = [
        "<!DOCTYPE html><html lang='en'><head>",
        "<meta charset='UTF-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        f"<title>MK Flu-Pipe v{html.escape(args.pipeline_version)} - Influenza Genomic Surveillance</title>",
        "<link rel='stylesheet' href='https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.3/css/bootstrap.min.css'>",
        "<link rel='stylesheet' href='https://cdnjs.cloudflare.com/ajax/libs/datatables/1.10.21/css/dataTables.bootstrap5.min.css'>",
        "<script src='https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js'></script>",
        f"<style>{css}</style></head><body>",
        "<div class='container-fluid py-3 px-4'>",
        "<div class='mb-4 text-center'>",
        f"<h1 class='mb-1 fw-bold'>MK Flu-Pipe <span class='badge bg-primary fs-4'>v{html.escape(args.pipeline_version)}</span></h1>",
        f"<p class='text-muted mb-0'>LACEN/AL &nbsp;|&nbsp; Run: <strong>{html.escape(datetime.now().strftime('%Y-%m-%d %H:%M'))}</strong> &nbsp;|&nbsp; IRMA Module: <strong>{html.escape(args.irma_module)}</strong></p>",
        "</div>",
        cards,
        "<div class='sec'>Typing and quality control</div>",
        render_html_table("tbl_typ", typing_rows, typing_fields),
        "<div class='sec'>Segment Coverage</div>",
        render_html_table(
            "tbl_depth",
            coverage_rows,
            ["sample", "type", "gene", "segment", "cov_mean", "cov_min", "cov_max", "positions_covered", "ref_length"],
        ),
        "<div class='sec'>Co-infection and subtype mixing</div>",
        render_html_table("tbl_coinf", coinfection_rows, ["sample", "coinfection_status", "details"]),
        "<div class='sec'>Antiviral resistance (FluSurver/WHO)</div>",
        render_html_table(
            "tbl_antiv",
            antiviral_rows,
            list(antiviral_rows[0].keys()) if antiviral_rows else ["sample", "gene", "aa_position", "wt_who", "mut_who", "alt_observed", "frequency", "depth_total", "drug", "significance", "nomenclature"],
        ),
    ]
    if h5_rows:
        html_parts.extend(
            [
                "<div class='sec'>H5 virulence markers</div>",
                render_html_table("tbl_h5", h5_rows, list(h5_rows[0].keys())),
            ]
        )
    if protein_summary_rows:
        html_parts.extend(
            [
                "<div class='sec'>Protein Mutations</div>",
                render_protein_table("tbl_fvc", protein_summary_rows),
            ]
        )

    html_parts.extend(
        [
            "<hr class='mt-4'>",
            f"<footer class='text-muted' style='font-size:.75rem'>MK Flu-Pipe v{html.escape(args.pipeline_version)} &mdash; Jean Phellipe M. Nascimento &mdash; LACEN/AL &mdash; {html.escape(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))}</footer>",
            "</div>",
            "<script src='https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js'></script>",
            "<script src='https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.3/js/bootstrap.bundle.min.js'></script>",
            "<script src='https://cdnjs.cloudflare.com/ajax/libs/datatables/1.10.21/js/jquery.dataTables.min.js'></script>",
            "<script src='https://cdnjs.cloudflare.com/ajax/libs/datatables/1.10.21/js/dataTables.bootstrap5.min.js'></script>",
            "<script>$(function(){['tbl_typ','tbl_depth','tbl_coinf','tbl_antiv','tbl_h5','tbl_fvc'].forEach(function(id){if($('#'+id).length)$('#'+id).DataTable({pageLength:25,language:{search:'Search:',info:'_START_-_END_ de _TOTAL_',paginate:{previous:'Ant',next:'Prox'}}});});});</script>",
            "<script>function dlExcel(tid){var wb=XLSX.utils.book_new();var ws=XLSX.utils.table_to_sheet(document.getElementById(tid));XLSX.utils.book_append_sheet(wb,ws,tid);XLSX.writeFile(wb,'MKFluPipe_'+tid+'_'+new Date().toISOString().slice(0,10)+'.xlsx');}</script>",
            "</body></html>",
        ]
    )
    html_path.write_text("".join(html_parts), encoding="utf-8")

    readme_path = out_dir / "README_outputs.txt"
    readme_text = (
        "MK Flu-Pipe Nextflow -- Surveillance_Outputs\n"
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
        "FILES\n"
        "  typing_results.tsv              Typing consolidated (BLAST + Nextclade + QC)\n"
        "  coverage_per_segment.tsv        Coverage per segment and sample\n"
        "  run_summary.tsv                 Typing + run metadata\n"
        "  run_summary.json                Same summary in JSON\n"
        "  multisample_consensus.fasta     Multi-sample FASTA\n"
        "  surveillance_report.html        Static HTML report\n"
        "  coinfection/coinfection_report.tsv          Co-infection analysis per sample\n"
        "  antiviral_resistance/antiviral_resistance.tsv Antiviral resistance mutations\n"
        "  h5_virulence/h5_virulence_markers.tsv       H5 virulence markers (if detected)\n"
        "  GISAID_ready/                   FASTA + CSV template EpiFlu (if location provided)\n"
        "  README_outputs.txt              This file\n"
    )
    readme_path.write_text(readme_text, encoding="utf-8")

    if args.gisaid_location.strip():
        gisaid_dir = out_dir / "GISAID_ready"
        gisaid_dir.mkdir(parents=True, exist_ok=True)
        gisaid_fasta = gisaid_dir / "gisaid_sequences.fasta"
        rewrite_gisaid_sequences(run_summary_rows, consensus_map, gisaid_fasta)
        gisaid_csv = gisaid_dir / "gisaid_metadata.csv"
        with gisaid_csv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "Isolate_Id",
                    "Isolate_Name",
                    "Subtype",
                    "Lineage",
                    "Passage_History",
                    "Location",
                    "Host",
                    "Originating_Lab",
                    "Submitting_Lab",
                    "Authors",
                    "Submission_Date",
                    "Collection_Date",
                    "Sequence_Length",
                ]
            )
            for row in run_summary_rows:
                virus_name, subtype = build_gisaid_name(row, args.gisaid_location.strip(), args.gisaid_year or str(date.today().year))
                writer.writerow(
                    [
                        row["sample"],
                        virus_name,
                        subtype,
                        row["nextclade_clade"],
                        "Original",
                        args.gisaid_location.strip(),
                        "Human",
                        "LACEN-AL",
                        "LACEN-AL",
                        "Jean Nascimento et al.",
                        today,
                        "",
                        "",
                    ]
                )

    log_lines = [
        f"Typing rows: {len(typing_rows)}",
        f"Coverage rows: {len(coverage_rows)}",
        f"Coinfection rows: {len(coinfection_rows)}",
        f"Antiviral rows: {len(antiviral_rows)}",
        f"H5 rows: {len(h5_rows)}",
        f"Protein mutation rows: {len(fullvar_rows)}",
    ]
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
