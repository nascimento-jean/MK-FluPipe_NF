process RUN_ANTIVIRAL_RESISTANCE {
    tag 'antiviral-resistance'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'antiviral_resistance/*', mode: 'copy', overwrite: true

    input:
    path antiviral_db
    path canonical_refs_dir
    path blast_summary
    path sample_dirs, stageAs: 'input??/*'

    output:
    path('antiviral_resistance/antiviral_resistance.tsv'), emit: report
    path('reports/antiviral_resistance.log'), emit: log

    script:
    """
    mkdir -p antiviral_resistance reports

    python "${projectDir}/bin/run_antiviral_resistance.py"         --db "${antiviral_db}"         --canonical-refs-dir "${canonical_refs_dir}"         --blast-summary "${blast_summary}"         --output-tsv antiviral_resistance/antiviral_resistance.tsv         --log-file reports/antiviral_resistance.log         --ivar-freq ${params.ivar_freq}         input*/*
    """

    stub:
    """
    mkdir -p antiviral_resistance reports
    printf 'sample\tgene\taa_position\twt_who\tmut_who\talt_observed\tfrequency\tdepth_total\tdrug\tsignificance\tnomenclature\n' > antiviral_resistance/antiviral_resistance.tsv
    echo 'stub short-read antiviral resistance report' > reports/antiviral_resistance.log
    """
}
