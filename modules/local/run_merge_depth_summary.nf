process RUN_MERGE_DEPTH_SUMMARY {
    tag 'merge-depth-summary'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'depth_summary.tsv', mode: 'copy', overwrite: true

    input:
    path depth_stats, stageAs: 'input??/*'

    output:
    path("depth_summary.tsv"), emit: summary

    script:
    """
    printf 'sample	segment	cov_mean	cov_min	cov_max	positions_covered	ref_length
' > depth_summary.tsv
    for row in input*/*; do
        tail -n +2 "\$row" >> depth_summary.tsv
    done
    """
}
