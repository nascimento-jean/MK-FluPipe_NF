process RUN_MEDAKA_VARIANTS {
    tag { meta.id }
    label 'medaka_tools'
    publishDir "${params.output_dir}", pattern: 'variant_calls/*', mode: 'copy', overwrite: true

    input:
    val ready
    tuple val(meta), path(irma_dir)

    output:
    path("variant_calls/${meta.id}"), emit: sample_dirs
    path("reports/${meta.id}.medaka.log"), emit: reports

    script:
    def sortThreads = Math.max(1, task.cpus as int)
    """
    mkdir -p "variant_calls/${meta.id}" reports
    : > "reports/${meta.id}.medaka.log"

    medaka_cmd=""
    if command -v medaka_haploid_variant >/dev/null 2>&1; then
        medaka_cmd="medaka_haploid_variant"
    elif command -v medaka_variant >/dev/null 2>&1; then
        medaka_cmd="medaka_variant"
    else
        printf 'Neither medaka_haploid_variant nor medaka_variant was found in PATH\n' >> "reports/${meta.id}.medaka.log"
        exit 127
    fi

    called=0
    bam_seen=0
    for bam in "${irma_dir}"/*.bam; do
        [ -f "\$bam" ] || continue
        bam_seen=\$((bam_seen + 1))
        seg_name=\$(basename "\$bam" .bam)
        ref_fa="${irma_dir}/\${seg_name}.fasta"

        if [ ! -f "\$ref_fa" ]; then
            seg_idx=""
            case "\$seg_name" in
                *_PB2) seg_idx="1" ;;
                *_PB1) seg_idx="2" ;;
                *_PA)  seg_idx="3" ;;
                *_HA*|*_HA) seg_idx="4" ;;
                *_NP)  seg_idx="5" ;;
                *_NA*|*_NA) seg_idx="6" ;;
                *_MP)  seg_idx="7" ;;
                *_NS)  seg_idx="8" ;;
            esac

            if [ -n "\$seg_idx" ]; then
                ref_fa="${irma_dir}/amended_consensus/${meta.id}_\${seg_idx}.fa"
            fi
        fi

        if [ ! -f "\$ref_fa" ]; then
            printf 'Skipping %s: reference FASTA not found for BAM %s\n' "\$seg_name" "\$bam" >> "reports/${meta.id}.medaka.log"
            continue
        fi

        sorted_bam="variant_calls/${meta.id}/\${seg_name}.sorted.bam"
        samtools sort -@ ${sortThreads} -o "\$sorted_bam" "\$bam" >> "reports/${meta.id}.medaka.log" 2>&1
        samtools index "\$sorted_bam" >> "reports/${meta.id}.medaka.log" 2>&1

        medaka_out="variant_calls/${meta.id}/\${seg_name}_medaka"
        mkdir -p "\$medaka_out"

        "\$medaka_cmd" \
            -i "\$sorted_bam" \
            -r "\$ref_fa" \
            -o "\$medaka_out" \
            -t ${task.cpus} \
            >> "reports/${meta.id}.medaka.log" 2>&1

        [[ -f "\$medaka_out/medaka.annotated.vcf" ]] && called=\$((called + 1))
    done

    if [[ "\$bam_seen" -eq 0 ]]; then
        printf 'No BAM files found for %s in %s\n' "${meta.id}" "${irma_dir}" >> "reports/${meta.id}.medaka.log"
    fi
    """
}
