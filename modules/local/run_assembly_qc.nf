process RUN_ASSEMBLY_QC {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/qc_reports/assembly_qc/reports", pattern: '*.assembly_qc.tsv', mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(irma_dir), path(consensus_fasta)

    output:
    path("*.assembly_qc.tsv"), emit: qc_rows

    script:
    """
    python "${projectDir}/bin/assembly_qc.py" \
        --sample-id "${meta.id}" \
        --consensus-fasta "${consensus_fasta}" \
        --irma-dir "${irma_dir}" \
        --min-coverage "${params.min_coverage}" \
        --max-n-pct "${params.max_n_pct}" \
        --min-segments "${params.min_segments}" \
        --output-tsv "${meta.id}.assembly_qc.tsv"
    """
}
