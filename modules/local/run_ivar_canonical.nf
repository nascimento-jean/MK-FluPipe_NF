process RUN_IVAR_CANONICAL {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}", pattern: 'variant_calls/*', mode: 'copy', overwrite: true

    input:
    val ready
    path canonical_refs_dir
    tuple val(meta), path(reads), val(flu_type), val(subtype_ha), val(subtype_na)

    output:
    path("variant_calls/${meta.id}"), emit: sample_dirs
    path("reports/${meta.id}.ivar.log"), emit: reports

    script:
    def sortThreads = Math.max(1, (task.cpus as int) - 1)
    def alignThreads = Math.max(1, (task.cpus as int) - sortThreads)
    """
    mkdir -p "variant_calls/${meta.id}" reports

    if [[ "${flu_type}" == "B" || "${flu_type}" == "Not determined" || -z "${subtype_na}" || "${subtype_na}" == "nd" ]]; then
        printf 'Skipping canonical iVar for %s (type=%s, NA=%s)\n' "${meta.id}" "${flu_type}" "${subtype_na}" > "reports/${meta.id}.ivar.log"
        exit 0
    fi

    subtype_key="N1"
    case "${subtype_na}" in
        N2*|n2*) subtype_key="N2" ;;
    esac
    called=0

    for gene in NA PA MP; do
        ref_fa="${canonical_refs_dir}/\${subtype_key}_\${gene}.fa"
        [ -f "\$ref_fa" ] || continue

        prefix="variant_calls/${meta.id}/${meta.id}_\${gene}_canonical"
        canon_bam="\${prefix}.bam"
        ivar_prefix="\${prefix}_ivar"

        if [[ ${meta.layout} == "paired" && ${reads.size()} -gt 1 ]]; then
            minimap2 -ax sr -t ${alignThreads} "\$ref_fa" "${reads[0]}" "${reads[1]}" 2>> "reports/${meta.id}.ivar.log" \
                | samtools sort -@ ${sortThreads} -o "\$canon_bam" - >> "reports/${meta.id}.ivar.log" 2>&1
        else
            minimap2 -ax sr -t ${alignThreads} "\$ref_fa" "${reads[0]}" 2>> "reports/${meta.id}.ivar.log" \
                | samtools sort -@ ${sortThreads} -o "\$canon_bam" - >> "reports/${meta.id}.ivar.log" 2>&1
        fi

        [ -s "\$canon_bam" ] || continue
        samtools index "\$canon_bam" >> "reports/${meta.id}.ivar.log" 2>&1

        samtools mpileup \
            -A -d 0 -B -Q 0 \
            --reference "\$ref_fa" \
            "\$canon_bam" 2>> "reports/${meta.id}.ivar.log" \
        | ivar variants \
            -p "\$ivar_prefix" \
            -q 20 \
            -t ${params.ivar_freq} \
            -m ${params.ivar_depth} \
            -r "\$ref_fa" \
            >> "reports/${meta.id}.ivar.log" 2>&1

        if [[ -s "\${ivar_prefix}.tsv" ]]; then
            called=\$((called + 1))
        fi
    done
    """

    stub:
    """
    mkdir -p "variant_calls/${meta.id}" reports
    printf 'REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\n' > "variant_calls/${meta.id}/${meta.id}_NA_canonical_ivar.tsv"
    echo 'stub canonical iVar calls for ${meta.id}' > "reports/${meta.id}.ivar.log"
    """
}
