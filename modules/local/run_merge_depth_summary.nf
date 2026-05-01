process RUN_MERGE_DEPTH_SUMMARY {
    tag 'merge-depth-summary'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'depth_summary.tsv', mode: 'copy', overwrite: true

    input:
    path depth_stats

    output:
    path("depth_summary.tsv"), emit: summary

    script:
    def stagedDepthStats = depth_stats.collect { "\"${it.getFileName().toString()}\"" }.join(' ')

    """
    printf 'sample\tsegment\tcov_mean\tcov_min\tcov_max\tpositions_covered\tref_length\n' > depth_summary.tsv
    for row in ${stagedDepthStats}; do
        tail -n +2 "\$row" >> depth_summary.tsv
    done
    """
}
