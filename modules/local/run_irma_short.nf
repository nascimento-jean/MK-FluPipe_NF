process RUN_IRMA_SHORT {
    tag { meta.id }
    label 'irma_tools'
    publishDir "${params.output_dir}", pattern: 'irma_runs_short/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}", pattern: 'assembly_final/*', mode: 'copy', overwrite: true

    input:
    val ready
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("irma_runs_short/${meta.id}"), emit: irma_dirs
    tuple val(meta), path("assembly_final/${meta.id}.fasta"), emit: consensus_fastas
    tuple val(meta), path("reports/*"), emit: reports
    tuple val(meta), path("${meta.id}_irma_manifest.tsv"), emit: manifests

    script:
    def irmaModule = params.irma_module as String

    if( meta.layout == 'paired' ) {
        """
        mkdir -p irma_runs_short assembly_final reports

        IRMA "${irmaModule}" \
            "${reads[0]}" \
            "${reads[1]}" \
            "irma_runs_short/${meta.id}" \
            > "reports/${meta.id}.irma.log" 2>&1

        if ls "irma_runs_short/${meta.id}/amended_consensus/"*.fa >/dev/null 2>&1; then
            cat "irma_runs_short/${meta.id}/amended_consensus/"*.fa > "assembly_final/${meta.id}.fasta"
            sed -i '/^>/!s/[^ACTGactg]/N/g' "assembly_final/${meta.id}.fasta"
        else
            echo "IRMA finished but amended_consensus is missing or empty for ${meta.id}" >&2
            exit 1
        fi

        {
            printf 'sample_id\\tlayout\\tseq_type\\tread1\\tread2\\tirma_module\\tirma_dir\\tconsensus_fasta\\tlog_file\\n'
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
                "${meta.id}" \
                "${meta.layout}" \
                "${meta.seq_type}" \
                "${reads[0]}" \
                "${reads[1]}" \
                "${irmaModule}" \
                "irma_runs_short/${meta.id}" \
                "assembly_final/${meta.id}.fasta" \
                "reports/${meta.id}.irma.log"
        } > "${meta.id}_irma_manifest.tsv"
        """
    }
    else {
        """
        mkdir -p irma_runs_short assembly_final reports

        IRMA "${irmaModule}" \
            "${reads[0]}" \
            "irma_runs_short/${meta.id}" \
            > "reports/${meta.id}.irma.log" 2>&1

        if ls "irma_runs_short/${meta.id}/amended_consensus/"*.fa >/dev/null 2>&1; then
            cat "irma_runs_short/${meta.id}/amended_consensus/"*.fa > "assembly_final/${meta.id}.fasta"
            sed -i '/^>/!s/[^ACTGactg]/N/g' "assembly_final/${meta.id}.fasta"
        else
            echo "IRMA finished but amended_consensus is missing or empty for ${meta.id}" >&2
            exit 1
        fi

        {
            printf 'sample_id\\tlayout\\tseq_type\\tread1\\tread2\\tirma_module\\tirma_dir\\tconsensus_fasta\\tlog_file\\n'
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
                "${meta.id}" \
                "${meta.layout}" \
                "${meta.seq_type}" \
                "${reads[0]}" \
                "-" \
                "${irmaModule}" \
                "irma_runs_short/${meta.id}" \
                "assembly_final/${meta.id}.fasta" \
                "reports/${meta.id}.irma.log"
        } > "${meta.id}_irma_manifest.tsv"
        """
    }
}
