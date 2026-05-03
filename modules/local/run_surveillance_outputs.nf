process RUN_SURVEILLANCE_OUTPUTS {
    tag 'surveillance-outputs'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}", mode: 'copy', overwrite: true

    input:
    path dependency_files, stageAs: 'deps??/*'
    path consensus_fastas, stageAs: 'cons??/*'
    path irma_dirs, stageAs: 'irma??/*'
    val irma_module
    val pipeline_version
    val gisaid_location
    val gisaid_year

    output:
    path('Surveillance_Outputs'), emit: outputs
    path('reports/surveillance_outputs.log'), emit: log

    script:
    """
    mkdir -p Surveillance_Outputs reports

    python "${projectDir}/bin/run_surveillance_outputs.py"         --output-dir Surveillance_Outputs         --log-file reports/surveillance_outputs.log         --irma-module "${irma_module}"         --pipeline-version "${pipeline_version}"         --gisaid-location "${gisaid_location}"         --gisaid-year "${gisaid_year}"         --dependencies deps*/*         --consensus-fastas cons*/*         --irma-dirs irma*/*
    """
}
