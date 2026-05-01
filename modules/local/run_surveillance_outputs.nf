process RUN_SURVEILLANCE_OUTPUTS {
    tag 'surveillance-outputs'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/surveillance_outputs", pattern: 'reports/*', mode: 'copy', overwrite: true

    input:
    path dependency_files
    path consensus_fastas
    path irma_dirs
    val irma_module
    val pipeline_version
    val gisaid_location
    val gisaid_year

    output:
    path('Surveillance_Outputs'), emit: outputs
    path('reports/surveillance_outputs.log'), emit: log

    script:
    def stagedDeps = dependency_files.collect { "\"${it.getFileName().toString()}\"" }.join(' ')
    def stagedConsensus = consensus_fastas.collect { "\"${it.getFileName().toString()}\"" }.join(' ')
    def stagedIrmaDirs = irma_dirs.collect { "\"${it.getFileName().toString()}\"" }.join(' ')

    """
    mkdir -p Surveillance_Outputs reports

    python "${projectDir}/bin/run_surveillance_outputs.py" \
        --output-dir Surveillance_Outputs \
        --log-file reports/surveillance_outputs.log \
        --irma-module "${irma_module}" \
        --pipeline-version "${pipeline_version}" \
        --gisaid-location "${gisaid_location}" \
        --gisaid-year "${gisaid_year}" \
        --dependencies ${stagedDeps} \
        --consensus-fastas ${stagedConsensus} \
        --irma-dirs ${stagedIrmaDirs}
    """
}
