process RUN_EXTRACT_SEGMENTS {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final/segments", pattern: 'single_segment*/*', mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(consensus_fasta)

    output:
    path("single_segment*/*.fasta"), emit: segment_files

    script:
    """
    python "${projectDir}/bin/extract_segments.py" \
        --input-fasta "${consensus_fasta}" \
        --sample-id "${meta.id}" \
        --output-dir .
    """

    stub:
    """
    for seg in 1 2 3 4 5 6 7 8; do
        mkdir -p "single_segment\${seg}"
        printf '>%s_%s\\nACGTACGTACGT\\n' "${meta.id}" "\${seg}" > "single_segment\${seg}/${meta.id}_segment_\${seg}.fasta"
    done
    """
}
