process RUN_MULTIQC {
    tag 'multiqc'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/qc_reports/multiqc", pattern: 'MultiQC/*', mode: 'copy', overwrite: true

    input:
    val ready
    path qc_artifacts, stageAs: 'input??/*'

    output:
    path('MultiQC/multiqc_report.html'), emit: report
    path('MultiQC/multiqc_data'), emit: data_dir
    path('reports/multiqc.log'), emit: log

    script:
    """
    mkdir -p MultiQC reports

    multiqc \
        --force \
        --outdir MultiQC \
        input* \
        > reports/multiqc.log 2>&1
    """
}