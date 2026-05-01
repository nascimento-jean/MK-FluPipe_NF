process RUN_MEDAKA_CANONICAL {
    tag { meta.id }
    label 'medaka_tools'
    publishDir "${params.output_dir}", pattern: 'variant_calls_canonical_long/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/medaka_canonical", pattern: 'reports/*', mode: 'copy', overwrite: true

    input:
    val ready
    path canonical_refs_dir
    tuple val(meta), path(reads), val(flu_type), val(subtype_ha), val(subtype_na)

    output:
    path("variant_calls_canonical_long/${meta.id}"), emit: sample_dirs
    path("reports/${meta.id}.medaka_canonical.log"), emit: reports

    script:
    def sortThreads = Math.max(1, (task.cpus as int) - 1)
    def alignThreads = Math.max(1, (task.cpus as int) - sortThreads)
    """
    mkdir -p "variant_calls_canonical_long/${meta.id}" reports
    : > "reports/${meta.id}.medaka_canonical.log"

    if [[ "${flu_type}" == "B" || "${flu_type}" == "Not determined" || -z "${subtype_na}" || "${subtype_na}" == "nd" ]]; then
        printf 'sample_id\tflu_type\tsubtype_HA\tsubtype_NA\tgenes_called\n' > "variant_calls_canonical_long/${meta.id}/${meta.id}_medaka_canonical_manifest.tsv"
        printf '%s\t%s\t%s\t%s\t0\n' "${meta.id}" "${flu_type}" "${subtype_ha}" "${subtype_na}" >> "variant_calls_canonical_long/${meta.id}/${meta.id}_medaka_canonical_manifest.tsv"
        printf 'Skipping canonical Medaka for %s (type=%s, NA=%s)\n' "${meta.id}" "${flu_type}" "${subtype_na}" > "reports/${meta.id}.medaka_canonical.log"
        exit 0
    fi

    medaka_cmd=""
    if command -v medaka_haploid_variant >/dev/null 2>&1; then
        medaka_cmd="medaka_haploid_variant"
    elif command -v medaka_variant >/dev/null 2>&1; then
        medaka_cmd="medaka_variant"
    else
        printf 'Neither medaka_haploid_variant nor medaka_variant was found in PATH\n' >> "reports/${meta.id}.medaka_canonical.log"
        exit 127
    fi

    subtype_key="N1"
    case "${subtype_na}" in
        N2*|n2*) subtype_key="N2" ;;
    esac

    called=0
    for gene in NA PA MP; do
        ref_fa="${canonical_refs_dir}/\${subtype_key}_\${gene}.fa"
        [ -f "\$ref_fa" ] || continue

        prefix="variant_calls_canonical_long/${meta.id}/${meta.id}_\${gene}_canonical"
        canon_bam="\${prefix}.bam"
        medaka_out="\${prefix}_medaka"

        minimap2 -ax map-ont -t ${alignThreads} "\$ref_fa" "${reads[0]}" 2>> "reports/${meta.id}.medaka_canonical.log" \
            | samtools sort -@ ${sortThreads} -o "\$canon_bam" - >> "reports/${meta.id}.medaka_canonical.log" 2>&1

        [ -s "\$canon_bam" ] || continue
        samtools index "\$canon_bam" >> "reports/${meta.id}.medaka_canonical.log" 2>&1

        mkdir -p "\$medaka_out"
        "\$medaka_cmd" \
            -i "\$canon_bam" \
            -r "\$ref_fa" \
            -o "\$medaka_out" \
            -t ${task.cpus} \
            >> "reports/${meta.id}.medaka_canonical.log" 2>&1

        [[ -f "\$medaka_out/medaka.annotated.vcf" ]] && called=\$((called + 1))
    done

    printf 'sample_id\tflu_type\tsubtype_HA\tsubtype_NA\tgenes_called\n' > "variant_calls_canonical_long/${meta.id}/${meta.id}_medaka_canonical_manifest.tsv"
    printf '%s\t%s\t%s\t%s\t%s\n' "${meta.id}" "${flu_type}" "${subtype_ha}" "${subtype_na}" "\$called" >> "variant_calls_canonical_long/${meta.id}/${meta.id}_medaka_canonical_manifest.tsv"
    """
}
