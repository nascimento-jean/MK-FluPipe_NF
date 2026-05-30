process RUN_H5_VIRULENCE {
    tag 'h5-virulence'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'h5_virulence/*', mode: 'copy', overwrite: true

    input:
    path blast_summary
    path irma_dirs, stageAs: 'input??/*'

    output:
    path('h5_virulence/h5_virulence_markers.tsv'), emit: report
    path('reports/h5_virulence.log'), emit: log

    script:
    """
    mkdir -p h5_virulence reports

    python "${projectDir}/bin/run_h5_virulence.py"         --blast-summary "${blast_summary}"         --output-tsv h5_virulence/h5_virulence_markers.tsv         --log-file reports/h5_virulence.log         input*/*
    """

    stub:
    """
    mkdir -p h5_virulence reports
    printf 'sample\tsegment\tposition\tmarker\tobserved\tnote\n' > h5_virulence/h5_virulence_markers.tsv
    echo 'stub H5 virulence report' > reports/h5_virulence.log
    """
}
