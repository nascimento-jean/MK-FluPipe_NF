process RUN_BLAST_TYPING {
    tag 'blast-typing'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'blast_results/*', mode: 'copy', overwrite: true

    input:
    val ready
    path blast_db_dir
    path segment_files

    output:
    path('blast_results/*'), emit: reports
    path('blast_results/blast_typing_summary.tsv'), emit: summary

    script:
    def stagedSegmentFiles = segment_files.collect { "\"${it.getFileName().toString()}\"" }.join(' ')

    """
    mkdir -p blast_results

    python "${projectDir}/bin/run_blast_typing.py" \
        --blast-db-prefix "${blast_db_dir}/influenza_blast_db" \
        --output-dir blast_results \
        --log-file blast_results/blast_typing.log \
        --threads ${task.cpus} \
        ${stagedSegmentFiles}
    """
}
