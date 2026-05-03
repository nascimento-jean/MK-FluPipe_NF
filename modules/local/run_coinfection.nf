process RUN_COINFECTION {
    tag 'coinfection'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'coinfection/*', mode: 'copy', overwrite: true

    input:
    path irma_dirs, stageAs: 'input??/*'

    output:
    path('coinfection/coinfection_report.tsv'), emit: report
    path('reports/coinfection.log'), emit: log

    script:
    """
    mkdir -p coinfection reports

    python "${projectDir}/bin/run_coinfection.py" \
        --output-tsv coinfection/coinfection_report.tsv \
        --log-file reports/coinfection.log \
        --minority-freq ${params.minority_freq} \
        --coinfection-pct ${params.coinfection_pct} \
        input*/*
    """
}