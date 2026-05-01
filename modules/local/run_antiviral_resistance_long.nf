process RUN_ANTIVIRAL_RESISTANCE_LONG {
    tag 'antiviral-resistance-long'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'antiviral_resistance/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/antiviral_resistance", pattern: 'reports/*', mode: 'copy', overwrite: true

    input:
    path antiviral_db
    path canonical_refs_dir
    path blast_summary
    path sample_dirs

    output:
    path('antiviral_resistance/antiviral_resistance.tsv'), emit: report
    path('reports/antiviral_resistance.log'), emit: log

    script:
    def stagedSampleDirs = sample_dirs.collect { "\"${it.getFileName().toString()}\"" }.join(' ')

    """
    mkdir -p antiviral_resistance reports

    python "${projectDir}/bin/run_antiviral_resistance.py" \
        --db "${antiviral_db}" \
        --canonical-refs-dir "${canonical_refs_dir}" \
        --blast-summary "${blast_summary}" \
        --output-tsv antiviral_resistance/antiviral_resistance.tsv \
        --log-file reports/antiviral_resistance.log \
        --ivar-freq ${params.ivar_freq} \
        ${stagedSampleDirs}
    """
}
