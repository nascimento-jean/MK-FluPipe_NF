process RUN_MERGE_SEGMENTS {
    tag 'merge-segments'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/assembly_final/segments", pattern: 'segment_*.fasta', mode: 'copy', overwrite: true

    input:
    path segment_files, stageAs: 'input??/*'

    output:
    path('segment_*.fasta'), emit: merged_segments

    script:
    """
    for seg in 1 2 3 4 5 6 7 8; do
        : > "segment_\${seg}.fasta"
        for fasta in input*/*; do
            case "\$fasta" in
                *"_segment_\${seg}.fasta")
                    if [ -s "\$fasta" ]; then
                        cat "\$fasta" >> "segment_\${seg}.fasta"
                    fi
                    ;;
            esac
        done
    done
    """
}