process RUN_SAMTOOLS_DEPTH {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}", pattern: 'depth_per_position/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/samtools_depth", pattern: '*.depth_stats.tsv', mode: 'copy', overwrite: true

    input:
    val ready
    tuple val(meta), path(irma_dir)

    output:
    path("*.depth_stats.tsv"), emit: depth_stats

    script:
    """
    mkdir -p "depth_per_position/${meta.id}"

    bam_found=0

    for bam in "${irma_dir}"/*.bam; do
        [ -f "\$bam" ] || continue
        seg_name=\$(basename "\$bam" .bam)
        depth_file="depth_per_position/${meta.id}/\${seg_name}.depth.tsv"
        samtools depth -a -d 0 "\$bam" > "\$depth_file"
        if [ -s "\$depth_file" ]; then
            bam_found=\$((bam_found + 1))
        fi
    done

    for bam in "${irma_dir}"/intermediate/*.bam "${irma_dir}"/secondary/*.bam; do
        [ -f "\$bam" ] || continue
        seg_name=\$(basename "\$bam" .bam)
        depth_file="depth_per_position/${meta.id}/\${seg_name}_secondary.depth.tsv"
        samtools depth -a -d 0 "\$bam" > "\$depth_file"
    done

    {
        printf 'sample\tsegment\tcov_mean\tcov_min\tcov_max\tpositions_covered\tref_length\n'
        for depth_file in depth_per_position/${meta.id}/*.depth.tsv; do
            [ -f "\$depth_file" ] || continue
            seg_label=\$(basename "\$depth_file" .depth.tsv)
            stats=\$(awk '
                BEGIN { sum=0; n=0; min=999999; max=0; covered=0 }
                {
                    depth=\$3+0
                    sum += depth
                    n++
                    if (depth < min) min = depth
                    if (depth > max) max = depth
                    if (depth > 0) covered++
                }
                END {
                    if (n > 0)
                        printf "%.1f\t%d\t%d\t%d\t%d", sum/n, min, max, covered, n
                    else
                        printf "0\t0\t0\t0\t0"
                }' "\$depth_file")
            printf '%s\t%s\t%s\n' "${meta.id}" "\$seg_label" "\$stats"
        done
    } > "${meta.id}.depth_stats.tsv"
    """

    stub:
    """
    mkdir -p "depth_per_position/${meta.id}"
    printf 'sample\\tsegment\\tcov_mean\\tcov_min\\tcov_max\\tpositions_covered\\tref_length\\n' > "${meta.id}.depth_stats.tsv"
    printf '%s\\tA_HA\\t100.0\\t50\\t150\\t12\\t12\\n' "${meta.id}" >> "${meta.id}.depth_stats.tsv"
    printf 'A_HA\\t1\\t100\\n' > "depth_per_position/${meta.id}/A_HA.depth.tsv"
    """
}
