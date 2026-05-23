#!/usr/bin/env python3

import argparse
import csv
import html
import json
import shutil
from collections import Counter, defaultdict
from itertools import zip_longest
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
    "M": "7",
}

MISSING = {"", "N/A", "ERRO", "nd", "no_result", "sem_resultado", "sem resultado", "-", "—", "â€”", "Ã¢â‚¬â€"}


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


def text(value):
    return str(value or "").strip()


def is_missing(value):
    return text(value) in MISSING


def normalize_dash(value):
    return "-" if is_missing(value) else text(value)


def normalize_hit(value):
    return "-" if is_missing(value) else text(value)


def dedupe_rows(rows, key):
    seen = {}
    for row in rows:
        row_key = text(row.get(key))
        if row_key:
            seen[row_key] = row
    return list(seen.values())


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
    return "|".join(unique) if unique else "-"


def infer_type_gene_segment(segment_name: str):
    parts = segment_name.split("_")
    flu_type = parts[0] if parts else "-"
    gene = parts[1] if len(parts) > 1 else "-"
    seg_num = GENE_TO_SEGMENT.get(gene, "-")
    return flu_type, gene, seg_num


def blast_classification(row):
    flu_type = text(row.get("type_blast", ""))
    subtype_ha = text(row.get("subtype_HA", "nd"))
    subtype_na = text(row.get("subtype_NA", "nd"))
    if flu_type == "B":
        return normalize_dash(subtype_ha)
    if subtype_ha not in {"", "nd"} and subtype_na not in {"", "nd"}:
        return f"{subtype_ha}{subtype_na}"
    if subtype_ha not in {"", "nd"}:
        return subtype_ha
    if subtype_na not in {"", "nd"}:
        return subtype_na
    return "-"


def pretty_header(name: str):
    mapping = {
        "subtype_HA": "Subtype Ha",
        "subtype_NA": "Subtype Na",
        "blast_classification": "Blast Classification",
        "nextclade_clade": "Nextclade Clade",
        "qc_nextclade": "Qc Nextclade",
        "classification_source": "Classification Source",
        "hit_blast_HA": "Hit Blast Ha",
        "hit_blast_NA": "Hit Blast Na",
        "qc_assembly": "Qc Assembly",
        "qc_detail": "Qc Detail",
        "coinfection_status": "Coinfection Status",
        "aa_position": "Aa Position",
        "wt_who": "Wt Who",
        "mut_who": "Mut Who",
        "alt_observed": "Alt Observed",
        "depth_total": "Depth Total",
        "n_mutations": "N Mutations",
        "mutations_by_segment": "Mutations By Segment",
        "q30_pct": "Q30 Pct",
        "retention_pct": "Retention Pct",
        "input_reads": "Input Reads",
        "output_reads": "Output Reads",
        "mean_length": "Mean Length",
        "sample_id": "Sample",
    }
    return mapping.get(name, name.replace("_", " ").title())


def badge(value: str):
    lv = text(value).lower()
    if lv in {"pass", "good", "ok", "wt"}:
        css = "bg-success"
    elif lv in {"warn", "mediocre", "verify"}:
        css = "bg-warning text-dark"
    elif lv in {"fail", "bad", "mut_detected"}:
        css = "bg-danger"
    elif "blast+nextclade" in lv:
        css = "bg-primary"
    elif lv == "variant":
        return '<span class="badge rounded-pill" style="background:#fd7e14;color:#fff">{}</span>'.format(html.escape(text(value)))
    else:
        css = "bg-secondary"
    return '<span class="badge {} rounded-pill">{}</span>'.format(css, html.escape(text(value)))


def render_cell(header: str, value):
    raw = text(value)
    display = "-" if raw == "nd" else raw
    header_lower = header.lower()
    if header_lower == "type":
        if display == "A":
            return '<span class="badge rounded-pill" style="background:#66CDAA;color:#1a4a3a">A</span>'
        if display == "B":
            return '<span class="badge rounded-pill" style="background:#FFDEAD;color:#5a3e1b">B</span>'
    if any(
        key in header_lower
        for key in ("qc_assembly", "qc_nextclade", "classification_source", "coinfection_status", "status")
    ):
        return badge(display)
    if header_lower == "nomenclature" and display not in {"N/D", "?", "", "-"}:
        clean = display.replace("?", "")
        if len(clean) >= 3 and clean[0] == clean[-1]:
            return '<span style="color:#198754;font-weight:600">{}</span>'.format(html.escape(display))
        return '<span style="color:#dc3545;font-weight:700">{}</span>'.format(html.escape(display))
    if header_lower == "mutation":
        if display in {".", "", "-"}:
            return '<span style="color:#adb5bd;font-size:.78rem">syn</span>'
        if "*" in display:
            return '<span style="color:#dc3545;font-weight:700;background:#ffeaea;padding:1px 4px;border-radius:3px">{}</span>'.format(html.escape(display))
        return '<span style="color:#6f42c1;font-weight:700">{}</span>'.format(html.escape(display))
    if header_lower in {"frequency", "retention_pct", "q30_pct"}:
        try:
            freq = float(display.replace("%", ""))
            color = "#198754" if freq >= 80 else "#fd7e14" if freq >= 30 else "#dc3545"
            return '<span style="color:{};font-weight:600">{}</span>'.format(color, html.escape(display))
        except Exception:
            pass
    if header_lower in {"qc_detail", "details"}:
        return '<div style="white-space:pre-line;text-align:left">{}</div>'.format(html.escape(display))
    return html.escape(display)


def render_download_button(table_id: str):
    return (
        '<div class="text-end mb-1">'
        '<button onclick="dlExcel(\'{}\')" class="btn btn-sm btn-outline-success" style="font-size:.75rem">'
        '&#11015; Download Excel</button></div>'
    ).format(table_id)


def render_empty_panel(message):
    return '<p class="text-muted fst-italic mb-0">{}</p>'.format(html.escape(message))


def render_html_table(table_id: str, rows, headers, caption=None):
    if not rows:
        return render_empty_panel("No data available for this section.")
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


def render_protein_table(table_id: str, rows):
    if not rows:
        return render_empty_panel("No protein mutation summary available.")
    headers = ["sample", "type", "segment", "n_mutations", "mutations_by_segment"]
    head_html = "".join(f"<th>{html.escape(pretty_header(h))}</th>" for h in headers)
    body = []
    for row in rows:
        body.append(
            "<tr>"
            f"<td>{html.escape(row['sample'])}</td>"
            f"<td>{render_cell('type', row['type'])}</td>"
            f"<td>{html.escape(row['segment'])}</td>"
            f"<td>{html.escape(row['n_mutations'])}</td>"
            f"<td style='text-align:left!important;color:#dc3545;font-weight:500'>{row['mutations_by_segment']}</td>"
            "</tr>"
        )
    return (
        render_download_button(table_id)
        + '<div class="table-responsive mt-1">'
        + f'<table id="{table_id}" class="table table-sm table-striped table-hover align-middle" style="font-size:.82rem">'
        + f'<thead class="table-dark"><tr>{head_html}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'
    )


def parse_fastp_json(path: Path):
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    summary = payload.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    filtering = payload.get("filtering_result", {})
    input_reads = before.get("total_reads") or 0
    output_reads = after.get("total_reads") or 0
    retention_pct = (output_reads / input_reads * 100) if input_reads else 0.0
    q30 = after.get("q30_rate")
    q30_pct = q30 * 100 if isinstance(q30, (int, float)) and q30 <= 1 else (q30 or 0)
    sample = path.name.replace(".fastp.json", "")
    return {
        "sample": sample,
        "step": "fastp",
        "input_reads": str(input_reads),
        "output_reads": str(output_reads),
        "retention_pct": f"{retention_pct:.2f}%",
        "q30_pct": f"{float(q30_pct):.2f}%",
        "mean_length": str(after.get("read1_mean_length") or after.get("mean_length") or ""),
        "notes": "passed_filter_reads={}".format(filtering.get("passed_filter_reads", output_reads)),
    }


def parse_filtlong_stats(path: Path):
    rows = read_tsv(path)
    if not rows:
        return None
    row = rows[0]
    sample = row.get("sample") or path.name.replace(".filtlong.stats.tsv", "")
    retention = row.get("pct_reads_retained") or row.get("retention_pct") or row.get("pct_reads_kept") or row.get("fraction_reads_retained") or ""
    if retention and not text(retention).endswith("%"):
        try:
            retention = f"{float(retention):.2f}%"
        except Exception:
            pass
    return {
        "sample": sample,
        "step": "filtlong",
        "input_reads": str(row.get("input_reads") or row.get("input_reads_total") or row.get("total_reads") or ""),
        "output_reads": str(row.get("kept_reads") or row.get("output_reads") or row.get("reads_retained") or ""),
        "retention_pct": str(retention),
        "q30_pct": "-",
        "mean_length": str(row.get("mean_length") or row.get("mean_read_length") or ""),
        "notes": "bases_retained={}".format(row.get("kept_bases") or row.get("bases_retained") or ""),
    }


def parse_host_depletion_stats(path: Path):
    rows = read_tsv(path)
    if not rows:
        return None
    row = rows[0]
    return {
        "sample": text(row.get("sample_id")),
        "input_reads": int(float(row.get("input_reads") or 0)),
        "output_reads": int(float(row.get("output_reads") or 0)),
        "retention_pct": f"{float(row.get('read_retention_pct') or 0):.2f}%",
        "input_mean_len": row.get("input_mean_len", ""),
        "output_mean_len": row.get("output_mean_len", ""),
        "layout": row.get("layout", ""),
        "seq_type": row.get("seq_type", ""),
    }


def build_preprocessing_rows(dependency_paths):
    rows = []
    for path in dependency_paths:
        if not path.exists() or path.is_dir():
            continue
        name = path.name
        parsed = None
        if name.endswith(".fastp.json"):
            parsed = parse_fastp_json(path)
        elif name.endswith(".filtlong.stats.tsv"):
            parsed = parse_filtlong_stats(path)
        if parsed:
            rows.append(parsed)
    rows.sort(key=lambda item: (item["sample"], item["step"]))
    return rows


def build_host_depletion_rows(dependency_paths):
    rows = []
    for path in dependency_paths:
        if not path.exists() or path.is_dir():
            continue
        if path.name.endswith(".host_depletion.stats.tsv"):
            parsed = parse_host_depletion_stats(path)
            if parsed:
                rows.append(parsed)
    rows.sort(key=lambda item: item["sample"])
    return dedupe_rows(rows, "sample")




def format_blast_hit(value):
    raw = text(value)
    if is_missing(raw):
        return "-"
    parts = raw.split("|")
    accession = ""
    organism = ""
    if len(parts) >= 5:
        accession = parts[3].strip()
        organism = parts[4].strip()
    elif len(parts) >= 2:
        accession = parts[-2].strip()
        organism = parts[-1].strip()
    if accession and organism:
        return f"Accession: {accession} | {organism}"
    return raw


def format_qc_detail(value):
    raw = text(value)
    if is_missing(raw):
        return "-"
    pieces = []
    for token in [part.strip() for part in raw.split("|") if part.strip()]:
        if token.startswith("segs:"):
            pieces.append(f"Segments: {token.split(':', 1)[1]}")
        elif token.startswith("low_cov:"):
            values = [item for item in token.split(":", 1)[1].split() if item]
            if values:
                pieces.append("Low Coverage: " + ", ".join(values))
        else:
            label, sep, content = token.partition(":")
            if sep:
                pretty = label.replace("_", " ").title()
                values = [item for item in content.split() if item]
                pieces.append(f"{pretty}: {', '.join(values) if values else content}")
            else:
                pieces.append(token)
    return "\n".join(pieces) if pieces else raw


def format_coinfection_detail(value):
    raw = text(value)
    if is_missing(raw):
        return "-"
    segs_analyzed = ""
    segs_flagged = ""
    alerts = ""
    for token in [part.strip() for part in raw.split("|") if part.strip()]:
        if token.startswith("segs_analyzed:"):
            segs_analyzed = token.split(":", 1)[1].strip()
        elif token.startswith("segs_flagged:"):
            segs_flagged = token.split(":", 1)[1].strip()
        elif token == "no_alert":
            alerts = ""
        else:
            alerts = token
    prefix = f"{segs_analyzed} segments analyzed" if segs_analyzed else "Segments analyzed"
    if not alerts or segs_flagged in {"", "0"}:
        return f"{prefix} and no alert"
    entries = []
    for item in [part for part in alerts.split() if part]:
        item = item.replace("pos,", " positions, ")
        entries.append(item)
    if len(entries) == 1:
        joined = entries[0]
    else:
        joined = ", ".join(entries[:-1]) + f" and {entries[-1]}"
    return f"{prefix} and alerts found in: {joined}"


def build_gisaid_name(row, location: str, year: str):
    subtype_ha = row.get("subtype_HA", "nd")
    subtype_na = row.get("subtype_NA", "nd")
    flu_type = row.get("type", "-")
    sample = row.get("sample", "sample")
    subtype = "B" if flu_type == "B" else f"{subtype_ha}{subtype_na}".replace("nd", "")
    virus_name = f"A/{location}/{sample}/{year}" if flu_type == "A" else f"B/{location}/{sample}/{year}"
    return virus_name, subtype or flu_type


def ensure_local_copy(source: Path, target: Path):
    if not source or not source.exists():
        return None
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(source, target)
    else:
        shutil.copy2(source, target)
    return target


def summarize_protein_mutations(fullvar_rows):
    grouped = defaultdict(list)
    for row in fullvar_rows:
        sample = row.get("sample", "")
        segment = row.get("segment", "")
        mutation = row.get("mutation", ".")
        gene = row.get("gene", "")
        if not sample or not segment or mutation in {".", "", "mutation"}:
            continue
        grouped[(sample, segment)].append((gene, mutation))

    rows = []
    for (sample, segment), items in sorted(grouped.items()):
        flu_type = segment.split("_")[1] if "_" in segment and len(segment.split("_")) > 1 else "-"
        display = []
        by_gene = defaultdict(list)
        for gene, mutation in items:
            by_gene[gene].append(mutation)
        for gene in sorted(by_gene):
            joined = ", ".join(f"{gene}:{m}" for m in by_gene[gene])
            display.append(f"<strong>{html.escape(gene)}</strong>: {html.escape(joined)}")
        rows.append(
            {
                "sample": sample,
                "type": flu_type,
                "segment": segment,
                "n_mutations": str(len(items)),
                "mutations_by_segment": " | ".join(display),
            }
        )
    return rows


def ordered_counts(rows, key, order):
    counter = Counter(text(row.get(key, "")) for row in rows)
    labels = [item for item in order if counter.get(item, 0) > 0]
    values = [counter[item] for item in labels]
    return labels, values


def main():
    args = parse_args()

    out_dir = Path(args.output_dir)
    log_path = Path(args.log_file)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    dependency_paths = [Path(p) for p in args.dependencies]
    dep_map = {path.name: path for path in dependency_paths}
    consensus_paths = [Path(p) for p in args.consensus_fastas]
    irma_paths = [Path(p) for p in args.irma_dirs]
    consensus_map = {path.stem: path for path in consensus_paths}

    blast_rows = dedupe_rows(read_tsv(dep_map.get("blast_typing_summary.tsv", Path("__missing__"))), "sample")
    nextclade_rows = dedupe_rows(read_tsv(dep_map.get("nextclade_summary.tsv", Path("__missing__"))), "sample")
    assembly_qc_rows = dedupe_rows(read_tsv(dep_map.get("assembly_qc_report.tsv", Path("__missing__"))), "sample")
    depth_rows = read_tsv(dep_map.get("depth_summary.tsv", Path("__missing__")))
    coinfection_rows = dedupe_rows(read_tsv(dep_map.get("coinfection_report.tsv", Path("__missing__"))), "sample")
    for row in coinfection_rows:
        row["details"] = format_coinfection_detail(row.get("details", "-"))
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
        sample_type = normalize_dash(row.get("type_blast", "-"))
        clade_display = normalize_dash(next_row.get("clade_display", "-"))
        if sample_type.lower() == "not determined":
            clade_display = "-"
        qc_nextclade = text(next_row.get("qc_status", "N/A")) or "N/A"
        source = "BLAST+Nextclade" if clade_display != "-" else "BLAST"

        typing_rows.append(
            {
                "sample": sample,
                "type": sample_type,
                "subtype_HA": normalize_dash(row.get("subtype_HA", "nd")),
                "subtype_NA": normalize_dash(row.get("subtype_NA", "nd")),
                "blast_classification": blast_classification(row),
                "nextclade_clade": clade_display,
                "qc_nextclade": qc_nextclade,
                "classification_source": source,
                "hit_blast_HA": format_blast_hit(row.get("hit_HA", "-")),
                "hit_blast_NA": format_blast_hit(row.get("hit_NA", "-")),
                "segments": infer_segments_from_fasta(consensus_map[sample]) if sample in consensus_map else "-",
                "qc_assembly": normalize_dash(qc_row.get("qc_assembly", "-")),
                "qc_detail": format_qc_detail(qc_row.get("qc_detail", "-")),
            }
        )

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
    write_tsv(out_dir / "typing_results.tsv", typing_rows, typing_fields)

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
    write_tsv(
        out_dir / "coverage_per_segment.tsv",
        coverage_rows,
        ["sample", "type", "gene", "segment", "cov_mean", "cov_min", "cov_max", "positions_covered", "ref_length"],
    )

    preprocessing_rows = build_preprocessing_rows(dependency_paths)
    preprocessing_fields = ["sample", "step", "input_reads", "output_reads", "retention_pct", "q30_pct", "mean_length", "notes"]
    write_tsv(out_dir / "preprocessing_summary.tsv", preprocessing_rows, preprocessing_fields)

    host_depletion_rows = build_host_depletion_rows(dependency_paths)
    host_fields = ["sample", "input_reads", "output_reads", "retention_pct", "input_mean_len", "output_mean_len", "layout", "seq_type"]
    write_tsv(out_dir / "host_depletion_summary.tsv", host_depletion_rows, host_fields)

    run_summary_rows = []
    today = date.today().isoformat()
    for row in typing_rows:
        merged = dict(row)
        merged["analysis_date"] = today
        merged["irma_module"] = args.irma_module
        merged["pipeline_version"] = args.pipeline_version
        run_summary_rows.append(merged)
    write_tsv(out_dir / "run_summary.tsv", run_summary_rows, typing_fields + ["analysis_date", "irma_module", "pipeline_version"])
    (out_dir / "run_summary.json").write_text(
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

    local_coinfection = ensure_local_copy(dep_map.get("coinfection_report.tsv"), out_dir / "coinfection" / "coinfection_report.tsv")
    local_antiviral = ensure_local_copy(dep_map.get("antiviral_resistance.tsv"), out_dir / "antiviral_resistance" / "antiviral_resistance.tsv")
    local_h5 = ensure_local_copy(dep_map.get("h5_virulence_markers.tsv"), out_dir / "h5_virulence" / "h5_virulence_markers.tsv")
    local_fullvar = ensure_local_copy(dep_map.get("all_samples_protein_mutations.tsv"), out_dir / "full_variant_calls" / "all_samples_protein_mutations.tsv")
    local_multiqc = ensure_local_copy(dep_map.get("multiqc_report.html"), out_dir / "multiqc_report.html")
    if local_coinfection:
        write_tsv(local_coinfection, coinfection_rows, ["sample", "coinfection_status", "details"])

    protein_summary_rows = summarize_protein_mutations(fullvar_rows)
    total_samples = len(typing_rows)
    total_pass = sum(1 for row in typing_rows if row.get("qc_assembly") == "PASS")
    total_warn_fail = sum(1 for row in typing_rows if row.get("qc_assembly") in {"WARN", "FAIL"})
    total_coinfection_warn = sum(1 for row in coinfection_rows if row.get("coinfection_status") == "WARN")
    total_protein_mutations = len(protein_summary_rows)
    total_h5 = len({row.get("sample", "") for row in h5_rows if row.get("sample")})

    type_counts = Counter(row.get("type", "Unknown") for row in typing_rows)
    subtype_chart_labels = [
        ("Not determined" if row.get("blast_classification") == "-" else row.get("blast_classification", "Unknown"))
        for row in typing_rows
        if row.get("type") == "A"
    ]
    subtype_counts = Counter(subtype_chart_labels)
    assembly_labels, assembly_values = ordered_counts(typing_rows, "qc_assembly", ["PASS", "WARN", "FAIL"])
    coinf_labels, coinf_values = ordered_counts(coinfection_rows, "coinfection_status", ["OK", "WARN"])

    segment_cov = defaultdict(list)
    for row in coverage_rows:
        try:
            segment_cov[row["segment"]].append(float(row["cov_mean"]))
        except Exception:
            continue
    ordered_cov_items = sorted(segment_cov.items(), key=lambda item: int(item[0]) if text(item[0]).isdigit() else 999)
    coverage_chart = {
        "labels": [seg for seg, _ in ordered_cov_items],
        "values": [round(sum(values) / len(values), 2) for _, values in ordered_cov_items],
    }

    host_chart = {
        "labels": [row["sample"] for row in host_depletion_rows],
        "before": [row["input_reads"] for row in host_depletion_rows],
        "after": [row["output_reads"] for row in host_depletion_rows],
    }

    steps_present = sorted({row["step"] for row in preprocessing_rows})
    if steps_present == ["fastp"]:
        preprocessing_title = "Fastp Summary"
        preprocessing_caption = "fastp summary"
    elif steps_present == ["filtlong"]:
        preprocessing_title = "Filtlong Summary"
        preprocessing_caption = "Filtlong summary"
    else:
        preprocessing_title = "Preprocessing Summary"
        preprocessing_caption = "Preprocessing summary"

    download_links = [
        ("Typing results", "typing_results.tsv"),
        ("Coverage per segment", "coverage_per_segment.tsv"),
        ("Preprocessing summary", "preprocessing_summary.tsv"),
        ("Host depletion summary", "host_depletion_summary.tsv"),
        ("Run summary (TSV)", "run_summary.tsv"),
        ("Run summary (JSON)", "run_summary.json"),
        ("Multisample consensus FASTA", "multisample_consensus.fasta"),
        ("Coinfection report", "coinfection/coinfection_report.tsv" if local_coinfection else ""),
        ("Antiviral resistance", "antiviral_resistance/antiviral_resistance.tsv" if local_antiviral else ""),
        ("H5 virulence markers", "h5_virulence/h5_virulence_markers.tsv" if local_h5 else ""),
        ("Protein mutations", "full_variant_calls/all_samples_protein_mutations.tsv" if local_fullvar else ""),
        ("MultiQC full report", "multiqc_report.html" if local_multiqc else ""),
    ]

    if args.gisaid_location.strip():
        gisaid_dir = out_dir / "GISAID_ready"
        gisaid_dir.mkdir(parents=True, exist_ok=True)
        gisaid_fasta = gisaid_dir / "gisaid_sequences.fasta"
        with gisaid_fasta.open("w", encoding="utf-8") as handle:
            for consensus_path in consensus_paths:
                if consensus_path.exists() and consensus_path.stat().st_size > 0:
                    for header, seq in parse_fasta(consensus_path):
                        handle.write(f">{header}\n{normalize_sequence(seq)}\n")

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
        download_links.extend(
            [
                ("GISAID FASTA", "GISAID_ready/gisaid_sequences.fasta"),
                ("GISAID metadata", "GISAID_ready/gisaid_metadata.csv"),
            ]
        )

    cards = (
        '<div class="row g-3 mb-4">'
        '<div class="col-6 col-md-2"><div class="card sc bl p-3"><div class="text-muted" style="font-size:.75rem">Samples</div><div class="fs-3 fw-bold text-primary">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc gr p-3"><div class="text-muted" style="font-size:.75rem">Assembly PASS</div><div class="fs-3 fw-bold text-success">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc or p-3"><div class="text-muted" style="font-size:.75rem">Assembly WARN/FAIL</div><div class="fs-3 fw-bold text-warning">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc rd p-3"><div class="text-muted" style="font-size:.75rem">Coinfection alerts</div><div class="fs-3 fw-bold text-danger">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc pu p-3"><div class="text-muted" style="font-size:.75rem">Protein mutation rows</div><div class="fs-3 fw-bold text-purple">{}</div></div></div>'
        '<div class="col-6 col-md-2"><div class="card sc te p-3"><div class="text-muted" style="font-size:.75rem">H5 flagged rows</div><div class="fs-3 fw-bold text-info">{}</div></div></div>'
        "</div>"
    ).format(total_samples, total_pass, total_warn_fail, total_coinfection_warn, total_protein_mutations, total_h5)

    charts = {
        "types": {"labels": list(type_counts.keys()), "values": list(type_counts.values())},
        "subtypes": {"labels": list(subtype_counts.keys()), "values": list(subtype_counts.values())},
        "assemblyQc": {"labels": assembly_labels, "values": assembly_values},
        "coinfection": {"labels": coinf_labels, "values": coinf_values},
        "coverage": coverage_chart,
        "hostDepletion": host_chart,
    }

    css = (
        "body{font-family:'Segoe UI',system-ui,sans-serif;background:#f6f8fb;color:#1f2937}"
        ".hero{background:linear-gradient(135deg,#0f172a,#1d4ed8);color:#fff;border-radius:1rem;padding:1.5rem 1.8rem;box-shadow:0 20px 40px rgba(15,23,42,.18)}"
        ".hero p{opacity:.92}.sc{border:0;border-left:5px solid;border-radius:.85rem;box-shadow:0 10px 24px rgba(15,23,42,.08)}"
        ".sc.bl{border-color:#0d6efd}.sc.gr{border-color:#198754}.sc.or{border-color:#fd7e14}.sc.rd{border-color:#dc3545}.sc.pu{border-color:#6f42c1}.sc.te{border-color:#0dcaf0}"
        ".text-purple{color:#6f42c1}.section{font-size:1rem;font-weight:800;color:#0f172a;text-transform:uppercase;letter-spacing:.04em;border-bottom:2px solid #dbe4f0;padding-bottom:.4rem;margin:1.2rem 0 .8rem}"
        ".panel{background:#fff;border-radius:1rem;box-shadow:0 10px 28px rgba(15,23,42,.08);padding:1rem 1rem 1.2rem}.nav-pills .nav-link{border-radius:999px;padding:.55rem 1rem;font-weight:700}"
        ".nav-pills .nav-link.active{background:#1d4ed8}thead.table-dark th{font-size:.78rem;white-space:nowrap;text-align:center!important}"
        "td{text-align:center!important;vertical-align:middle}.badge{font-size:.72rem}.chart-card{background:#fff;border-radius:1rem;box-shadow:0 8px 24px rgba(15,23,42,.08);padding:1rem;height:100%}"
        ".chart-card canvas{max-height:360px!important}.download-list a{display:flex;justify-content:space-between;align-items:center;background:#fff;border:1px solid #dbe4f0;border-radius:.8rem;padding:.75rem 1rem;text-decoration:none;color:#0f172a;margin-bottom:.6rem}"
        ".download-list a:hover{border-color:#1d4ed8;box-shadow:0 8px 20px rgba(29,78,216,.12)}"
    )

    html_parts = [
        "<!DOCTYPE html><html lang='en'><head>",
        "<meta charset='UTF-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        f"<title>MK Flu-Pipe v{html.escape(args.pipeline_version)} - Influenza Surveillance Dashboard</title>",
        "<link rel='stylesheet' href='https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.3/css/bootstrap.min.css'>",
        "<link rel='stylesheet' href='https://cdnjs.cloudflare.com/ajax/libs/datatables/1.10.21/css/dataTables.bootstrap5.min.css'>",
        "<script src='https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js'></script>",
        "<script src='https://cdn.jsdelivr.net/npm/chart.js'></script>",
        f"<style>{css}</style></head><body>",
        "<div class='container-fluid py-4 px-4'>",
        "<div class='hero mb-4 text-center'>",
        "<h1 class='display-5 fw-bold mb-2'>MK Flu-Pipe Nextflow Dashboard</h1>",
        f"<p class='mb-1'>Influenza genomic surveillance summary | IRMA module <strong>{html.escape(args.irma_module)}</strong> | Pipeline <strong>v{html.escape(args.pipeline_version)}</strong></p>",
        f"<p class='mb-0'>Generated on <strong>{html.escape(datetime.now().strftime('%Y-%m-%d %H:%M'))}</strong></p>",
        "</div>",
        cards,
        "<ul class='nav nav-pills mb-4 justify-content-center' id='dashboard-tabs' role='tablist'>",
        "<li class='nav-item'><button class='nav-link active' data-bs-toggle='pill' data-bs-target='#overview' type='button'>Overview</button></li>",
        "<li class='nav-item'><button class='nav-link' data-bs-toggle='pill' data-bs-target='#qc' type='button'>QC Dashboard</button></li>",
        "<li class='nav-item'><button class='nav-link' data-bs-toggle='pill' data-bs-target='#typing' type='button'>Typing</button></li>",
        "<li class='nav-item'><button class='nav-link' data-bs-toggle='pill' data-bs-target='#alerts' type='button'>Resistance & H5</button></li>",
        "<li class='nav-item'><button class='nav-link' data-bs-toggle='pill' data-bs-target='#coinfection' type='button'>Coinfection</button></li>",
        "<li class='nav-item'><button class='nav-link' data-bs-toggle='pill' data-bs-target='#mutations' type='button'>Protein Mutations</button></li>",
        "<li class='nav-item'><button class='nav-link' data-bs-toggle='pill' data-bs-target='#downloads' type='button'>Downloads</button></li>",
        "</ul>",
        "<div class='tab-content'>",
        "<div class='tab-pane fade show active' id='overview'>",
        "<div class='row g-3 mb-4'>",
        "<div class='col-lg-4'><div class='chart-card'><h6 class='fw-bold mb-3'>Type distribution</h6><canvas id='chartTypes'></canvas></div></div>",
        "<div class='col-lg-4'><div class='chart-card'><h6 class='fw-bold mb-3'>Subtype distribution</h6><canvas id='chartSubtypes'></canvas></div></div>",
        "<div class='col-lg-4'><div class='chart-card'><h6 class='fw-bold mb-3'>Assembly QC status</h6><canvas id='chartAssembly'></canvas></div></div>",
        "</div>",
        "<div class='row g-3 mb-4'>",
        "<div class='col-lg-6'><div class='chart-card'><h6 class='fw-bold mb-3'>Average coverage by segment</h6><canvas id='chartCoverage'></canvas></div></div>",
        "<div class='col-lg-6'><div class='chart-card'><h6 class='fw-bold mb-3'>Coinfection alerts</h6><canvas id='chartCoinfection'></canvas></div></div>",
        "</div>",
        "<div class='section'>Run summary</div><div class='panel'>",
        render_html_table("tbl_typ_summary", typing_rows, typing_fields, "Integrated typing and assembly overview"),
        "</div></div>",
        "<div class='tab-pane fade' id='qc'>",
        f"<div class='section'>{html.escape(preprocessing_title)}</div><div class='panel'>",
        render_html_table("tbl_preprocess", preprocessing_rows, preprocessing_fields, preprocessing_caption),
        "</div>",
        "<div class='row g-3 mb-4'>",
        "<div class='col-lg-8'><div class='chart-card'><h6 class='fw-bold mb-3'>Host depletion read counts</h6><canvas id='chartHostDepletion'></canvas></div></div>",
        "<div class='col-lg-4'><div class='panel h-100'><div class='section mt-0'>QC files</div>",
        ("<p class='mb-2'><a class='btn btn-outline-primary btn-sm' href='multiqc_report.html' target='_blank'>Open full MultiQC report</a></p>" if local_multiqc else "<p class='text-muted'>MultiQC report not available.</p>"),
        "<p class='text-muted mb-2' style='font-size:.88rem'>QC charts summarize host depletion and preprocessing metrics. Use MultiQC for detailed FastQC and fastp inspection.</p>",
        "</div></div></div>",
        "<div class='section'>Host depletion summary</div><div class='panel'>",
        render_html_table("tbl_hostdep", host_depletion_rows, host_fields, "Host depletion summary"),
        "</div></div>",
        "<div class='tab-pane fade' id='typing'>",
        "<div class='section'>Typing and classification</div><div class='panel'>",
        render_html_table("tbl_typing", typing_rows, typing_fields),
        "</div>",
        "<div class='section'>Coverage per segment</div><div class='panel'>",
        render_html_table("tbl_depth", coverage_rows, ["sample", "type", "gene", "segment", "cov_mean", "cov_min", "cov_max", "positions_covered", "ref_length"]),
        "</div></div>",
        "<div class='tab-pane fade' id='alerts'>",
        "<div class='section'>Antiviral resistance (FluSurver/WHO)</div><div class='panel'>",
        render_html_table("tbl_antiv", antiviral_rows, list(antiviral_rows[0].keys())) if antiviral_rows else render_empty_panel("No antiviral resistance markers were detected in this run."),
        "</div>",
        "<div class='section'>H5 virulence markers</div><div class='panel'>",
        render_html_table("tbl_h5", h5_rows, list(h5_rows[0].keys())) if h5_rows else render_empty_panel("No H5-associated virulence markers were reported in this run."),
        "</div></div>",
        "<div class='tab-pane fade' id='coinfection'>",
        "<div class='section'>Coinfection and subtype mixing</div><div class='panel'>",
        render_html_table("tbl_coinf", coinfection_rows, ["sample", "coinfection_status", "details"]),
        "</div></div>",
        "<div class='tab-pane fade' id='mutations'>",
        "<div class='section'>Protein mutations by segment</div><div class='panel'>",
        render_protein_table("tbl_fvc", protein_summary_rows),
        "</div></div>",
        "<div class='tab-pane fade' id='downloads'>",
        "<div class='section'>Download center</div><div class='download-list'>",
    ]
    for label, href in download_links:
        if href:
            html_parts.append(f"<a href='{html.escape(href)}' target='_blank'><span>{html.escape(label)}</span><span class='badge bg-light text-dark'>open</span></a>")
    html_parts.extend(
        [
            "</div></div></div>",
            "<hr class='mt-4'>",
            f"<footer class='text-muted text-center' style='font-size:.78rem'>MK Flu-Pipe v{html.escape(args.pipeline_version)} | Dashboard generated {html.escape(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))}</footer>",
            "</div>",
            "<script src='https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js'></script>",
            "<script src='https://cdnjs.cloudflare.com/ajax/libs/bootstrap/5.3.3/js/bootstrap.bundle.min.js'></script>",
            "<script src='https://cdnjs.cloudflare.com/ajax/libs/datatables/1.10.21/js/jquery.dataTables.min.js'></script>",
            "<script src='https://cdnjs.cloudflare.com/ajax/libs/datatables/1.10.21/js/dataTables.bootstrap5.min.js'></script>",
            "<script>function dt(id, opts){if($('#'+id).length)$('#'+id).DataTable(Object.assign({pageLength:25,order:[],language:{search:'Search:',info:'_START_-_END_ of _TOTAL_',paginate:{previous:'Prev',next:'Next'}}},opts||{}));}$(function(){dt('tbl_typ_summary');dt('tbl_preprocess');dt('tbl_hostdep');dt('tbl_typing');dt('tbl_depth');dt('tbl_antiv');dt('tbl_h5');dt('tbl_fvc');dt('tbl_coinf',{paging:false,info:false});});</script>",
            "<script>function dlExcel(tid){var wb=XLSX.utils.book_new();var ws=XLSX.utils.table_to_sheet(document.getElementById(tid));XLSX.utils.book_append_sheet(wb,ws,tid);XLSX.writeFile(wb,'MKFluPipe_'+tid+'_'+new Date().toISOString().slice(0,10)+'.xlsx');}</script>",
            f"<script>const charts={json.dumps(charts, ensure_ascii=False)};",
            """
function mkPie(id, labels, values, colors){
  const el=document.getElementById(id); if(!el) return;
  new Chart(el,{type:'doughnut',data:{labels:labels,datasets:[{data:values,backgroundColor:colors}]},options:{plugins:{legend:{position:'bottom'}},maintainAspectRatio:false}});
}
function mkBar(id, labels, values, color){
  const el=document.getElementById(id); if(!el) return;
  new Chart(el,{type:'bar',data:{labels:labels,datasets:[{data:values,backgroundColor:color,borderRadius:6}]},options:{scales:{y:{beginAtZero:true}},plugins:{legend:{display:false}},maintainAspectRatio:false}});
}
function mkGroupedBar(id, labels, beforeVals, afterVals){
  const el=document.getElementById(id); if(!el) return;
  new Chart(el,{
    type:'bar',
    data:{
      labels:labels,
      datasets:[
        {label:'Before host depletion',data:beforeVals,backgroundColor:'#6c757d',borderRadius:6},
        {label:'After host depletion',data:afterVals,backgroundColor:'#0d6efd',borderRadius:6}
      ]
    },
    options:{responsive:true,scales:{y:{beginAtZero:true}},plugins:{legend:{position:'bottom'}},maintainAspectRatio:false}
  });
}
mkPie('chartTypes', charts.types.labels, charts.types.values, ['#66CDAA','#FFDEAD','#adb5bd','#5b8def']);
mkBar('chartSubtypes', charts.subtypes.labels, charts.subtypes.values, '#5b8def');
const assemblyColors = charts.assemblyQc.labels.map(label => label === 'PASS' ? '#198754' : label === 'WARN' ? '#fd7e14' : '#dc3545');
mkPie('chartAssembly', charts.assemblyQc.labels, charts.assemblyQc.values, assemblyColors);
mkBar('chartCoverage', charts.coverage.labels, charts.coverage.values, '#198754');
const coinfColors = charts.coinfection.labels.map(label => label === 'OK' ? '#198754' : '#fd7e14');
mkPie('chartCoinfection', charts.coinfection.labels, charts.coinfection.values, coinfColors);
mkGroupedBar('chartHostDepletion', charts.hostDepletion.labels, charts.hostDepletion.before, charts.hostDepletion.after);
</script>""",
            "</body></html>",
        ]
    )

    (out_dir / "surveillance_report.html").write_text("".join(html_parts), encoding="utf-8")

    readme_text = (
        "MK Flu-Pipe Nextflow -- Surveillance_Outputs\n"
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
        "FILES\n"
        "  typing_results.tsv                         Integrated BLAST + Nextclade + assembly QC typing table\n"
        "  coverage_per_segment.tsv                   Coverage per segment and sample\n"
        "  preprocessing_summary.tsv                  fastp / Filtlong preprocessing summary\n"
        "  host_depletion_summary.tsv                 Read counts before and after host depletion\n"
        "  run_summary.tsv                            Typing + run metadata\n"
        "  run_summary.json                           Same summary in JSON\n"
        "  multisample_consensus.fasta                Multi-sample FASTA\n"
        "  surveillance_report.html                   Interactive dashboard report\n"
        "  multiqc_report.html                        Full MultiQC report (copied when available)\n"
        "  coinfection/coinfection_report.tsv         Co-infection analysis per sample\n"
        "  antiviral_resistance/antiviral_resistance.tsv Antiviral resistance mutations\n"
        "  h5_virulence/h5_virulence_markers.tsv      H5 virulence markers (if detected)\n"
        "  full_variant_calls/all_samples_protein_mutations.tsv Protein mutation calls\n"
        "  GISAID_ready/                              FASTA + CSV template EpiFlu (if location provided)\n"
        "  README_outputs.txt                         This file\n"
    )
    (out_dir / "README_outputs.txt").write_text(readme_text, encoding="utf-8")

    log_lines = [
        f"Typing rows: {len(typing_rows)}",
        f"Coverage rows: {len(coverage_rows)}",
        f"Preprocessing rows: {len(preprocessing_rows)}",
        f"Host depletion rows: {len(host_depletion_rows)}",
        f"Coinfection rows: {len(coinfection_rows)}",
        f"Antiviral rows: {len(antiviral_rows)}",
        f"H5 rows: {len(h5_rows)}",
        f"Protein mutation rows: {len(fullvar_rows)}",
        f"MultiQC report copied: {'yes' if local_multiqc else 'no'}",
        f"IRMA directories staged: {len(irma_paths)}",
    ]
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
