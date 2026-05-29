process RUN_MERGE_IRMA_STATUS {
    tag 'merge-irma-status'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}", pattern: 'pipeline_status/*', mode: 'copy', overwrite: true

    input:
    path status_files, stageAs: 'input??/*'

    output:
    path('pipeline_status/irma_status.tsv'), emit: report
    path('pipeline_status/irma_failures.tsv'), emit: failures

    script:
    """
    mkdir -p pipeline_status

    python3 - <<'PY'
import csv
from pathlib import Path

files = sorted(Path('.').glob('input*/*.tsv'))
rows = []
fieldnames = [
    'sample_id', 'layout', 'seq_type', 'status', 'reason',
    'read1', 'read2', 'irma_module', 'irma_dir', 'consensus_fasta', 'log_file'
]

for path in files:
    with path.open('r', encoding='utf-8', errors='ignore', newline='') as handle:
        rows.extend(csv.DictReader(handle, delimiter='\\t'))

rows.sort(key=lambda row: row.get('sample_id', ''))

report_path = Path('pipeline_status/irma_status.tsv')
with report_path.open('w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\\t')
    writer.writeheader()
    writer.writerows(rows)

failure_rows = [row for row in rows if row.get('status') != 'success']
failures_path = Path('pipeline_status/irma_failures.tsv')
with failures_path.open('w', encoding='utf-8', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\\t')
    writer.writeheader()
    writer.writerows(failure_rows)
PY
    """
}
