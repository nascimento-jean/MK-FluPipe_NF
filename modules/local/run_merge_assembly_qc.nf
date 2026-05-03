process RUN_MERGE_ASSEMBLY_QC {
    tag 'merge-assembly-qc'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'assembly_qc_report.tsv', mode: 'copy', overwrite: true

    input:
    path qc_rows, stageAs: 'input??/*'

    output:
    path("assembly_qc_report.tsv"), emit: report

    script:
    """
    printf 'sample	qc_assembly	qc_detail
' > assembly_qc_report.tsv
    for row in input*/*; do
        tail -n +2 "\$row" >> assembly_qc_report.tsv
    done
    """
}
