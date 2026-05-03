process RUN_FASTQC {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/qc_reports/fastqc_raw", mode: 'copy', overwrite: true

    input:
    val ready
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_fastqc"), emit: report_dirs

    script:
    def readsArg = reads.collect { "\"${it}\"" }.join(' ')

    """
    mkdir -p "${meta.id}_fastqc"

    fastqc \
        --threads ${task.cpus} \
        --outdir "${meta.id}_fastqc" \
        ${readsArg}
    """
}
