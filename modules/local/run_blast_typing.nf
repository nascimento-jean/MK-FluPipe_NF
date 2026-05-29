process RUN_BLAST_TYPING {
    tag 'blast-typing'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final", pattern: 'blast_results/*', mode: 'copy', overwrite: true

    input:
    val ready
    path blast_db_dir
    path segment_files, stageAs: 'input??/*'

    output:
    path('blast_results/*'), emit: reports
    path('blast_results/blast_typing_summary.tsv'), emit: summary

    script:
    """
    mkdir -p blast_results

    python "${projectDir}/bin/run_blast_typing.py"         --blast-db-prefix "${blast_db_dir}/influenza_blast_db"         --output-dir blast_results         --log-file blast_results/blast_typing.log         --threads ${task.cpus}         input*/*
    """

    stub:
    """
    mkdir -p blast_results
    printf 'sample\\ttype_blast\\tsubtype_HA\\tsubtype_NA\\thit_HA\\thit_NA\\n' > blast_results/blast_typing_summary.tsv
    for fasta in input*/*_segment_4.fasta; do
        [ -f "\$fasta" ] || continue
        sample=\$(basename "\$fasta" _segment_4.fasta)
        printf '%s\\tA\\tH3\\tN2\\tgi|1|gb|LK054740|Influenza\\tgi|2|gb|KY925263|Influenza\\n' "\$sample" >> blast_results/blast_typing_summary.tsv
    done
    echo 'stub BLAST typing log' > blast_results/blast_typing.log
    """
}
