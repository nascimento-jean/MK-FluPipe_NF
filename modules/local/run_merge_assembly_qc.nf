process RUN_MERGE_ASSEMBLY_QC {
    tag 'merge-assembly-qc'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'assembly_qc_report.tsv', mode: 'copy', overwrite: true

    input:
    path qc_rows

    output:
    path("assembly_qc_report.tsv"), emit: report

    script:
    def stagedQcRows = qc_rows.collect { "\"${it.getFileName().toString()}\"" }.join(' ')

    """
    printf 'sample\tqc_assembly\tqc_detail\n' > assembly_qc_report.tsv
    for row in ${stagedQcRows}; do
        tail -n +2 "\$row" >> assembly_qc_report.tsv
    done
    """
}
