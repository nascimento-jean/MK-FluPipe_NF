process PREPARE_REFSEQ_SEGMENTS {
    tag 'refseq-segments'
    label 'mk_flu_tools'

    input:
    path blast_summary
    val seq_mode

    output:
    path 'refseq_segments', emit: refs_dir
    path 'refseq_segments.log', emit: log

    script:
    def refseqSegmentsDir = params.refseq_segments_dir as String

    """
    mkdir -p "${refseqSegmentsDir}"

    python "${projectDir}/bin/prepare_refseq_segments.py" \
        --blast-summary "${blast_summary}" \
        --output-dir "${refseqSegmentsDir}" \
        --seq-mode "${seq_mode}" \
        --threads ${task.cpus} \
        --log-file refseq_segments.log

    mkdir -p refseq_segments
    for ref in "${refseqSegmentsDir}"/*; do
        [ -f "\$ref" ] || continue
        ln -sf "\$ref" "refseq_segments/\$(basename "\$ref")"
    done
    """

    stub:
    """
    mkdir -p refseq_segments
    for ref in A_HA A_NA A_MP B_HA B_NA B_MP; do
        printf '>%s\\nACGTACGTACGTACGT\\n' "\$ref" > "refseq_segments/\${ref}.fa"
    done
    echo 'stub refseq segments' > refseq_segments.log
    """
}
