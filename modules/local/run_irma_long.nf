process RUN_IRMA_LONG {
    tag { meta.id }
    label 'irma_tools'
    publishDir "${params.output_dir}", pattern: 'irma_runs_long/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}", pattern: 'irma_failures_long/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}", pattern: 'assembly_final/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}", pattern: 'pipeline_status/*', mode: 'copy', overwrite: true

    input:
    val ready
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("irma_runs_long/${meta.id}"), emit: irma_dirs, optional: true
    tuple val(meta), path("assembly_final/${meta.id}.fasta"), emit: consensus_fastas, optional: true
    tuple val(meta), path("reports/*"), emit: reports
    tuple val(meta), path("${meta.id}_irma_manifest.tsv"), emit: manifests, optional: true
    tuple val(meta), path("pipeline_status/${meta.id}.irma.tsv"), emit: statuses

    script:
    def irmaModule = params.irma_module as String

    """
    mkdir -p tmp_irma irma_runs_long irma_failures_long assembly_final reports pipeline_status

    tmp_dir="tmp_irma/${meta.id}"
    success_dir="irma_runs_long/${meta.id}"
    failure_dir="irma_failures_long/${meta.id}"
    status="success"
    reason="-"

    if ! IRMA "${irmaModule}" \
        "${reads[0]}" \
        "\$tmp_dir" \
        > "reports/${meta.id}.irma.log" 2>&1; then
        status="failed"
        reason="irma_exit_nonzero"
    elif ! ls "\$tmp_dir/amended_consensus/"*.fa >/dev/null 2>&1; then
        status="failed"
        reason="amended_consensus_missing_or_empty"
    else
        mv "\$tmp_dir" "\$success_dir"
        cat "\$success_dir/amended_consensus/"*.fa > "assembly_final/${meta.id}.fasta"
        sed -i '/^>/!s/[^ACTGactg]/N/g' "assembly_final/${meta.id}.fasta"

        {
            printf 'sample_id\\tlayout\\tseq_type\\tread1\\tread2\\tirma_module\\tirma_dir\\tconsensus_fasta\\tlog_file\\n'
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
                "${meta.id}" \
                "${meta.layout}" \
                "${meta.seq_type}" \
                "${reads[0]}" \
                "-" \
                "${irmaModule}" \
                "\$success_dir" \
                "assembly_final/${meta.id}.fasta" \
                "reports/${meta.id}.irma.log"
        } > "${meta.id}_irma_manifest.tsv"
    fi

    if [ "\$status" = "failed" ]; then
        rm -f "assembly_final/${meta.id}.fasta" "${meta.id}_irma_manifest.tsv"
        if [ -d "\$tmp_dir" ]; then
            mv "\$tmp_dir" "\$failure_dir"
        else
            mkdir -p "\$failure_dir"
        fi
        cp "reports/${meta.id}.irma.log" "\$failure_dir/" 2>/dev/null || true
    fi

    {
        printf 'sample_id\\tlayout\\tseq_type\\tstatus\\treason\\tread1\\tread2\\tirma_module\\tirma_dir\\tconsensus_fasta\\tlog_file\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
            "${meta.id}" \
            "${meta.layout}" \
            "${meta.seq_type}" \
            "\$status" \
            "\$reason" \
            "${reads[0]}" \
            "-" \
            "${irmaModule}" \
            "\$([ "\$status" = "success" ] && printf '%s' "\$success_dir" || printf '%s' "\$failure_dir")" \
            "\$([ "\$status" = "success" ] && printf '%s' "assembly_final/${meta.id}.fasta" || printf '%s' '-')" \
            "reports/${meta.id}.irma.log"
    } > "pipeline_status/${meta.id}.irma.tsv"
    """
}
