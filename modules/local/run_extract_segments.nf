process RUN_EXTRACT_SEGMENTS {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final/segments", pattern: 'single_segment*/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/segment_extraction", pattern: '*_segments_manifest.tsv', mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(consensus_fasta)

    output:
    path("single_segment*/*.fasta"), emit: segment_files
    path("*_segments_manifest.tsv"), emit: manifests

    script:
    """
    python "${projectDir}/bin/extract_segments.py" \
        --input-fasta "${consensus_fasta}" \
        --sample-id "${meta.id}" \
        --output-dir .
    """
}
