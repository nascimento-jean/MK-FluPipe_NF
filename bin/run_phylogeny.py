#!/usr/bin/env python3
"""Build optional HA and NA influenza phylogenies with Augur."""

from __future__ import annotations

import argparse
import colorsys
import csv
import html
import json
import re
import shutil
import subprocess
from collections import defaultdict
from datetime import date
from pathlib import Path


SUMMARY_FIELDS = [
    "group",
    "type",
    "segment",
    "subtype",
    "pipeline_sequences",
    "context_sequences",
    "total_sequences",
    "status",
    "message",
    "auspice_json",
    "tree_html",
]


USER_SEQUENCE_COLOR = "#8B0000"


class TreeNode:
    def __init__(self, name="", length=0.0):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None
        self.x = 0.0
        self.y = 0.0


def parse_args():
    parser = argparse.ArgumentParser(description="Create Augur trees for influenza HA and NA segments")
    parser.add_argument("--blast-summary", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--context-fasta", required=True)
    parser.add_argument("--context-metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-sequences", type=int, default=3)
    parser.add_argument("segment_files", nargs="+")
    return parser.parse_args()


def read_delimited(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        first_line = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in first_line else ","
        return [
            {key.strip(): str(value or "").strip() for key, value in row.items()}
            for row in csv.DictReader(handle, delimiter=delimiter)
        ]


def parse_fasta(path: Path):
    records = {}
    if not path.exists() or path.stat().st_size == 0:
        return records
    header = None
    sequence = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header:
                records[header] = "".join(sequence).upper()
            header = line[1:].split()[0]
            sequence = []
        elif header:
            sequence.append(line)
    if header:
        records[header] = "".join(sequence).upper()
    return records


def validate_context(rows, records):
    if not rows and not records:
        return
    errors = []
    if not rows:
        errors.append("context FASTA was supplied without context metadata rows")
    if rows and not records:
        errors.append("context metadata was supplied without context FASTA records")
    allowed_missing = {"", "-", "ND", "N/A", "NOT DETERMINED"}
    for line_number, row in enumerate(rows, start=2):
        strain = row.get("strain", row.get("sample_name", "")).strip()
        collection_date = row.get("collection_date", row.get("date", "")).strip()
        flu_type = row.get("type", "").strip().upper()
        segment = row.get("segment", "").strip().upper()
        if not strain:
            errors.append(f"line {line_number}: strain is required")
        elif strain not in records:
            errors.append(f"line {line_number}: strain '{strain}' is absent from context FASTA")
        try:
            date.fromisoformat(collection_date)
        except ValueError:
            errors.append(f"line {line_number}: collection_date '{collection_date}' must use YYYY-MM-DD")
        if flu_type not in {"A", "B"}:
            errors.append(f"line {line_number}: type must be A or B")
        if segment not in {"HA", "NA"}:
            errors.append(f"line {line_number}: segment must be HA or NA")
        subtype = (
            row.get("subtype_HA", row.get("subtype_ha", ""))
            if segment == "HA"
            else row.get("subtype_NA", row.get("subtype_na", ""))
        ).strip().upper()
        if flu_type == "A" and subtype in allowed_missing:
            errors.append(f"line {line_number}: influenza A {segment} records require a subtype")
    metadata_strains = {
        row.get("strain", row.get("sample_name", "")).strip() for row in rows
    }
    missing_metadata = sorted(set(records) - metadata_strains)
    if missing_metadata:
        errors.append("context FASTA records missing metadata: " + ", ".join(missing_metadata))
    if errors:
        raise SystemExit("Invalid phylogeny context: " + "; ".join(errors))


def normalize_value(value: str):
    return str(value or "").strip()


def group_for_record(flu_type: str, segment: str, subtype_ha: str, subtype_na: str):
    flu_type = normalize_value(flu_type).upper()
    segment = normalize_value(segment).upper()
    subtype_ha = normalize_value(subtype_ha).upper()
    subtype_na = normalize_value(subtype_na).upper()
    invalid = {"", "-", "ND", "N/A", "NOT DETERMINED"}

    if segment not in {"HA", "NA"} or flu_type not in {"A", "B"}:
        return None
    if flu_type == "B":
        return f"B_{segment}", flu_type, segment, "-"
    subtype = subtype_ha if segment == "HA" else subtype_na
    if subtype in invalid:
        return None
    return f"A_{subtype}_{segment}", flu_type, segment, subtype


def safe_identifier(identifier: str):
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", identifier.strip())
    return cleaned.strip("_") or "unnamed_sequence"


def segment_number_to_name(path: Path):
    match = re.search(r"_segment_([46])\.fasta$", path.name)
    if not match:
        return None
    return "HA" if match.group(1) == "4" else "NA"


def display_group(metadata, source):
    if source == "MK-FluPipe":
        return "User Sequences"
    state = str(metadata.get("state", "") or "").replace("_", " ").strip()
    state = re.sub(r"\s+", " ", state)
    return state or "State not available"


def hex_to_rgb(color: str):
    color = color.lstrip("#")
    return tuple(int(color[index:index + 2], 16) / 255.0 for index in (0, 2, 4))


def rgb_to_lab(red: float, green: float, blue: float):
    def pivot_rgb(value):
        return ((value + 0.055) / 1.055) ** 2.4 if value > 0.04045 else value / 12.92

    def pivot_xyz(value):
        return value ** (1 / 3) if value > 0.008856 else (7.787 * value) + (16 / 116)

    red, green, blue = [pivot_rgb(value) for value in (red, green, blue)]
    x = (red * 0.4124 + green * 0.3576 + blue * 0.1805) / 0.95047
    y = (red * 0.2126 + green * 0.7152 + blue * 0.0722) / 1.00000
    z = (red * 0.0193 + green * 0.1192 + blue * 0.9505) / 1.08883
    x, y, z = [pivot_xyz(value) for value in (x, y, z)]
    return (116 * y - 16, 500 * (x - y), 200 * (y - z))


def color_distance(color_a: str, color_b: str):
    lab_a = rgb_to_lab(*hex_to_rgb(color_a))
    lab_b = rgb_to_lab(*hex_to_rgb(color_b))
    return sum((left - right) ** 2 for left, right in zip(lab_a, lab_b)) ** 0.5


def generate_color_candidates():
    candidates = []
    for lightness in (0.50, 0.38, 0.64, 0.28, 0.72):
        for saturation in (0.95, 0.80, 0.65):
            for step in range(48):
                hue = (step / 48.0 + 0.07) % 1.0
                red, green, blue = colorsys.hls_to_rgb(hue, lightness, saturation)
                color = "#{:02X}{:02X}{:02X}".format(
                    round(red * 255),
                    round(green * 255),
                    round(blue * 255),
                )
                if color.upper() != USER_SEQUENCE_COLOR:
                    candidates.append(color)
    return list(dict.fromkeys(candidates))


def choose_display_group_colors(groups):
    available = generate_color_candidates()
    chosen = [USER_SEQUENCE_COLOR]
    color_map = {}
    for group in groups:
        best_color = max(
            available,
            key=lambda candidate: min(color_distance(candidate, used) for used in chosen),
        )
        color_map[group] = best_color
        chosen.append(best_color)
        available.remove(best_color)
    return color_map


def append_record(groups, group_info, identifier, sequence, metadata, source):
    if not group_info or not sequence:
        return
    group_name, flu_type, segment, subtype = group_info
    group = groups[group_name]
    tip = safe_identifier(identifier)
    existing = {record["strain"] for record in group}
    if tip in existing:
        tip = safe_identifier(f"{tip}_{source}")
    group.append(
        {
            "strain": tip,
            "sequence": sequence,
            "date": metadata.get("collection_date", metadata.get("date", "")),
            "country": metadata.get("country", ""),
            "state": metadata.get("state", ""),
            "city": metadata.get("city", ""),
            "source": source,
            "display_group": display_group(metadata, source),
            "type": flu_type,
            "segment": segment,
            "subtype": subtype,
        }
    )


def write_group_inputs(group_dir: Path, records):
    group_dir.mkdir(parents=True, exist_ok=True)
    with (group_dir / "sequences.fasta").open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record['strain']}\n{record['sequence']}\n")
    fields = ["strain", "date", "country", "state", "city", "source", "display_group", "type", "segment", "subtype"]
    with (group_dir / "metadata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows({key: record.get(key, "") for key in fields} for record in records)


def write_display_group_colors(path: Path, records):
    groups = sorted(
        {
            record.get("display_group", "").strip()
            for record in records
            if record.get("display_group", "").strip()
            and record.get("display_group", "").strip() != "User Sequences"
        }
    )
    color_map = choose_display_group_colors(groups)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"display_group\tUser Sequences\t{USER_SEQUENCE_COLOR}\n")
        for group in groups:
            handle.write(f"display_group\t{group}\t{color_map[group]}\n")


def read_color_map(path: Path):
    color_map = {}
    if not path.exists():
        return color_map
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 3 and parts[0] == "display_group":
            color_map[parts[1]] = parts[2]
    return color_map


def parse_newick(path: Path):
    text = path.read_text(encoding="utf-8", errors="ignore").strip()
    index = 0

    def skip_ws():
        nonlocal index
        while index < len(text) and text[index].isspace():
            index += 1

    def parse_label():
        nonlocal index
        skip_ws()
        if index < len(text) and text[index] in {"'", '"'}:
            quote = text[index]
            index += 1
            start = index
            while index < len(text) and text[index] != quote:
                index += 1
            label = text[start:index]
            if index < len(text):
                index += 1
            return label
        start = index
        while index < len(text) and text[index] not in ":,();":
            index += 1
        return text[start:index].strip()

    def parse_length():
        nonlocal index
        skip_ws()
        if index >= len(text) or text[index] != ":":
            return 0.0
        index += 1
        start = index
        while index < len(text) and text[index] not in ",();":
            index += 1
        try:
            return max(0.0, float(text[start:index].strip()))
        except ValueError:
            return 0.0

    def parse_node():
        nonlocal index
        skip_ws()
        node = TreeNode()
        if index < len(text) and text[index] == "(":
            index += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                skip_ws()
                if index < len(text) and text[index] == ",":
                    index += 1
                    continue
                if index < len(text) and text[index] == ")":
                    index += 1
                break
            node.name = parse_label()
            node.length = parse_length()
        else:
            node.name = parse_label()
            node.length = parse_length()
        return node

    return parse_node()


def iter_nodes(node):
    yield node
    for child in node.children:
        yield from iter_nodes(child)


def tree_leaves(node):
    if not node.children:
        return [node]
    leaves = []
    for child in node.children:
        leaves.extend(tree_leaves(child))
    return leaves


def assign_tree_layout(root):
    leaves = tree_leaves(root)
    for order, leaf in enumerate(leaves):
        leaf.y = order

    def assign_x(node, parent_x=0.0):
        node.x = parent_x + node.length
        for child in node.children:
            assign_x(child, node.x)

    def assign_y(node):
        if node.children:
            for child in node.children:
                assign_y(child)
            node.y = sum(child.y for child in node.children) / len(node.children)

    assign_x(root, 0.0)
    assign_y(root)
    return leaves


def build_svg_tree(root, metadata_by_strain, color_map):
    leaves = assign_tree_layout(root)
    max_x = max((node.x for node in iter_nodes(root)), default=1.0) or 1.0
    width = 1180
    left = 28
    right = 260
    top = 36
    row_height = 26
    plot_width = width - left - right
    height = max(220, top * 2 + row_height * max(1, len(leaves) - 1))

    def sx(value):
        return left + (value / max_x * plot_width)

    def sy(value):
        return top + value * row_height

    elements = [
        f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img'>",
        "<style>.branch{stroke:#a8b0b0;stroke-width:3;fill:none;stroke-linecap:round}.tip-label{font:13px Arial,sans-serif;fill:#111827}.axis{stroke:#e5e7eb;stroke-width:1}</style>",
    ]

    for frac in [i / 5 for i in range(6)]:
        x = left + frac * plot_width
        elements.append(f"<line class='axis' x1='{x:.2f}' x2='{x:.2f}' y1='{top / 2:.2f}' y2='{height - top / 2:.2f}'/>")

    def draw(node):
        if not node.children:
            return
        child_ys = [sy(child.y) for child in node.children]
        x = sx(node.x)
        elements.append(f"<line class='branch' x1='{x:.2f}' x2='{x:.2f}' y1='{min(child_ys):.2f}' y2='{max(child_ys):.2f}'/>")
        for child in node.children:
            cx = sx(child.x)
            cy = sy(child.y)
            elements.append(f"<line class='branch' x1='{x:.2f}' x2='{cx:.2f}' y1='{cy:.2f}' y2='{cy:.2f}'/>")
            draw(child)

    draw(root)

    for leaf in leaves:
        meta = metadata_by_strain.get(leaf.name, {})
        group = meta.get("display_group", "State not available")
        color = color_map.get(group, "#9CA3AF")
        x = sx(leaf.x)
        y = sy(leaf.y)
        label = html.escape(leaf.name)
        title = html.escape(f"{leaf.name} | {group} | {meta.get('date', '')}")
        elements.append(f"<circle cx='{x:.2f}' cy='{y:.2f}' r='5.5' fill='{color}' stroke='#374151' stroke-width='.55'><title>{title}</title></circle>")
        elements.append(f"<text class='tip-label' x='{x + 10:.2f}' y='{y + 4:.2f}'>{label}</text>")

    elements.append("</svg>")
    return "\n".join(elements)


def write_tree_html(group_dir: Path, group_name: str, records):
    tree_path = group_dir / "tree.nwk"
    if not tree_path.exists() or tree_path.stat().st_size == 0:
        return ""
    metadata_by_strain = {record["strain"]: record for record in records}
    color_map = read_color_map(group_dir / "colors.tsv")
    root = parse_newick(tree_path)
    svg = build_svg_tree(root, metadata_by_strain, color_map)
    legend_items = []
    for group, color in sorted(color_map.items(), key=lambda item: (item[0] != "User Sequences", item[0])):
        legend_items.append(
            "<span class='legend-item'><span class='swatch' style='background:{}'></span>{}</span>".format(
                html.escape(color),
                html.escape(group),
            )
        )
    html_path = group_dir / f"{group_name}.html"
    html_path.write_text(
        """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>MK Flu-Pipe {group} phylogeny</title>
  <style>
    body{{font-family:Arial,sans-serif;margin:0;background:#f8fafc;color:#111827}}
    header{{padding:1rem 1.25rem;background:linear-gradient(135deg,#0f172a,#1d4ed8);color:#fff}}
    main{{padding:1rem 1.25rem}}
    .panel{{background:#fff;border-radius:14px;box-shadow:0 10px 25px rgba(15,23,42,.10);padding:1rem;overflow:auto}}
    .legend{{display:flex;flex-wrap:wrap;gap:.5rem 1rem;margin:.75rem 0 1rem}}
    .legend-item{{display:inline-flex;align-items:center;gap:.35rem;font-size:.9rem}}
    .swatch{{display:inline-block;width:14px;height:14px;border-radius:50%;border:1px solid rgba(17,24,39,.3)}}
    .links a{{display:inline-block;margin-right:.75rem;color:#1d4ed8;font-weight:700}}
    svg{{min-width:900px;width:100%;height:auto}}
  </style>
</head>
<body>
  <header>
    <h1>MK Flu-Pipe {group} phylogeny</h1>
    <p>Static offline tree generated from the Augur Newick output. Use the Auspice JSON for the fully interactive Nextstrain view.</p>
  </header>
  <main>
    <div class="links">
      <a href="{group}.json" target="_blank">Open Auspice JSON</a>
      <a href="tree.nwk" target="_blank">Open Newick tree</a>
    </div>
    <div class="legend">{legend}</div>
    <div class="panel">{svg}</div>
  </main>
</body>
</html>
""".format(group=html.escape(group_name), legend="".join(legend_items), svg=svg),
        encoding="utf-8",
    )
    return f"{group_name}/{group_name}.html"


def run_command(command, log_handle):
    log_handle.write("$ " + " ".join(command) + "\n")
    result = subprocess.run(command, capture_output=True, text=True)
    if result.stdout:
        log_handle.write(result.stdout)
        if not result.stdout.endswith("\n"):
            log_handle.write("\n")
    if result.stderr:
        log_handle.write(result.stderr)
        if not result.stderr.endswith("\n"):
            log_handle.write("\n")
    return result.returncode == 0


def build_tree(group_dir: Path, group_name: str, records, threads: int, log_handle):
    aligned = group_dir / "aligned.fasta"
    raw_tree = group_dir / "tree_raw.nwk"
    refined_tree = group_dir / "tree.nwk"
    node_data = group_dir / "branch_lengths.json"
    auspice_json = group_dir / f"{group_name}.json"
    auspice_config = group_dir / "auspice_config.json"
    colors_file = group_dir / "colors.tsv"
    commands = [
        ["augur", "align", "--sequences", str(group_dir / "sequences.fasta"), "--output", str(aligned), "--nthreads", str(threads), "--fill-gaps"],
        ["augur", "tree", "--alignment", str(aligned), "--method", "iqtree", "--nthreads", str(threads), "--output", str(raw_tree)],
    ]
    for command in commands:
        if not run_command(command, log_handle):
            return "FAILED", "Augur alignment or tree inference failed", "", ""

    time_command = [
        "augur", "refine", "--tree", str(raw_tree), "--alignment", str(aligned),
        "--metadata", str(group_dir / "metadata.tsv"), "--metadata-id-columns", "strain",
        "--timetree", "--output-tree", str(refined_tree), "--output-node-data", str(node_data),
    ]
    status = "PASS"
    message = "Time-scaled Augur tree generated"
    if not run_command(time_command, log_handle):
        fallback_command = [
            "augur", "refine", "--tree", str(raw_tree), "--alignment", str(aligned),
            "--metadata", str(group_dir / "metadata.tsv"), "--metadata-id-columns", "strain",
            "--output-tree", str(refined_tree), "--output-node-data", str(node_data),
        ]
        if not run_command(fallback_command, log_handle):
            return "FAILED", "Augur refinement failed", "", ""
        status = "WARN"
        message = "Tree generated without time scaling; inspect temporal signal"

    auspice_config.write_text(
        json.dumps(
            {
                "title": f"MK Flu-Pipe {group_name} phylogeny",
                "colorings": [
                    {
                        "key": "display_group",
                        "title": "State / Source",
                        "type": "categorical",
                    }
                ],
                "display_defaults": {"color_by": "display_group"},
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    write_display_group_colors(colors_file, records)
    export_command = [
        "augur", "export", "v2", "--tree", str(refined_tree), "--node-data", str(node_data),
        "--metadata", str(group_dir / "metadata.tsv"), "--metadata-id-columns", "strain",
        "--auspice-config", str(auspice_config),
        "--metadata-columns", "source", "display_group", "type", "segment", "subtype", "country", "state", "city", "date",
        "--colors", str(colors_file),
        "--output", str(auspice_json),
    ]
    if not run_command(export_command, log_handle):
        return "FAILED", "Auspice JSON export failed", "", ""
    tree_html = write_tree_html(group_dir, group_name, records)
    return status, message, f"{group_name}/{group_name}.json", tree_html


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = Path(args.log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    blast_rows = {row["sample"]: row for row in read_delimited(Path(args.blast_summary))}
    sample_metadata = {row["sample_name"]: row for row in read_delimited(Path(args.metadata))}
    groups = defaultdict(list)

    for segment_file in [Path(item) for item in args.segment_files]:
        segment = segment_number_to_name(segment_file)
        if not segment:
            continue
        match = re.match(r"(.+)_segment_[46]\.fasta$", segment_file.name)
        if not match:
            continue
        sample = match.group(1)
        typing = blast_rows.get(sample, {})
        group_info = group_for_record(
            typing.get("type_blast", ""),
            segment,
            typing.get("subtype_HA", ""),
            typing.get("subtype_NA", ""),
        )
        for _, sequence in parse_fasta(segment_file).items():
            append_record(groups, group_info, sample, sequence, sample_metadata.get(sample, {}), "MK-FluPipe")

    context_records = parse_fasta(Path(args.context_fasta))
    context_rows = read_delimited(Path(args.context_metadata))
    validate_context(context_rows, context_records)
    for row in context_rows:
        strain = row.get("strain", row.get("sample_name", ""))
        if not strain or strain not in context_records:
            continue
        group_info = group_for_record(
            row.get("type", ""),
            row.get("segment", ""),
            row.get("subtype_HA", row.get("subtype_ha", "")),
            row.get("subtype_NA", row.get("subtype_na", "")),
        )
        append_record(groups, group_info, strain, context_records[strain], row, row.get("source", "Context") or "Context")

    summaries = []
    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write(f"MK Flu-Pipe phylogeny started {date.today().isoformat()}\n")
        for group_name in sorted(groups):
            records = groups[group_name]
            pipeline_count = sum(record["source"] == "MK-FluPipe" for record in records)
            context_count = len(records) - pipeline_count
            group_dir = out_dir / group_name
            write_group_inputs(group_dir, records)
            type_value = records[0]["type"]
            segment_value = records[0]["segment"]
            subtype_value = records[0]["subtype"]
            if len(records) < args.min_sequences:
                status = "SKIPPED"
                message = f"At least {args.min_sequences} sequences are required"
                auspice_json = ""
                tree_html = ""
            else:
                log_handle.write(f"\n[{group_name}] {len(records)} sequences\n")
                status, message, auspice_json, tree_html = build_tree(group_dir, group_name, records, args.threads, log_handle)
            summaries.append(
                {
                    "group": group_name,
                    "type": type_value,
                    "segment": segment_value,
                    "subtype": subtype_value,
                    "pipeline_sequences": pipeline_count,
                    "context_sequences": context_count,
                    "total_sequences": len(records),
                    "status": status,
                    "message": message,
                    "auspice_json": auspice_json,
                    "tree_html": tree_html,
                }
            )

    with (out_dir / "phylogeny_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)


if __name__ == "__main__":
    main()
