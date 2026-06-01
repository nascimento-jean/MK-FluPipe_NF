process RUN_H5_AVIAN_ANALYSIS {
    tag 'h5-avian-analysis'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'h5_avian/*', mode: 'copy', overwrite: true

    input:
    path blast_summary
    path segment_files, stageAs: 'input??/*'

    output:
    path('h5_avian/h5_avian_summary.tsv'), emit: summary
    path('h5_avian/*'), emit: reports
    path('reports/h5_avian_analysis.log'), emit: log

    script:
    """
    mkdir -p h5_avian reports

    python3 "${projectDir}/bin/run_h5_avian_analysis.py" \
        --blast-summary "${blast_summary}" \
        --output-dir h5_avian \
        --log-file reports/h5_avian_analysis.log \
        input*/*
    """

    stub:
    """
    mkdir -p h5_avian reports
    printf 'sample\ttype\tsubtype_HA\tsubtype_NA\tsegments_available\tselected_for_h5_avian\tgenoflu_status\tgenoflu_genotype\tgenoflu_message\tflumut_status\tflumut_markers_detected\tflumut_message\n' > h5_avian/h5_avian_summary.tsv
    printf 'sample\ttype\tsubtype_HA\tsubtype_NA\tsegments_available\tselected_for_h5_avian\tgenoflu_status\tgenoflu_genotype\tgenoflu_message\n' > h5_avian/genoflu_summary.tsv
    printf 'sample\tmarker\n' > h5_avian/flumut_markers.tsv
    printf 'sample\tmutation\n' > h5_avian/flumut_mutations.tsv
    printf 'sample\treference\n' > h5_avian/flumut_literature.tsv
    echo 'stub H5 avian analysis' > reports/h5_avian_analysis.log
    """
}
