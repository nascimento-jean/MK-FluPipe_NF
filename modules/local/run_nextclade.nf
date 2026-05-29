process RUN_NEXTCLADE {
    tag 'nextclade'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'nextclade_results/*', mode: 'copy', overwrite: true

    input:
    val ready
    path blast_summary
    path segment4_files, stageAs: 'input??/*'

    output:
    path('nextclade_results/*'), emit: reports
    path('nextclade_results/nextclade_summary.tsv'), emit: summary

    script:
    def nextcladeDatasetsDir = params.nextclade_datasets_dir as String
    def nextcladeMaxDays = params.nextclade_max_days as Integer

    """
    mkdir -p nextclade_results

    python "${projectDir}/bin/run_nextclade_batch.py"         --blast-summary "${blast_summary}"         --datasets-root "${nextcladeDatasetsDir}"         --max-days ${nextcladeMaxDays}         --output-dir nextclade_results         --log-file nextclade_results/nextclade.log         input*/*
    """

    stub:
    """
    mkdir -p nextclade_results
    printf 'sample\\tclade_display\\tqc_status\\n' > nextclade_results/nextclade_summary.tsv
    awk -F '\\t' 'NR>1 { print \$1 "\\t3C.2a1b.2a.2\\tgood" }' "${blast_summary}" >> nextclade_results/nextclade_summary.tsv
    echo 'stub Nextclade log' > nextclade_results/nextclade.log
    """
}
