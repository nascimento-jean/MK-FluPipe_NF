process RUN_FULLVARCALL_SAMPLE {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}", pattern: 'full_variant_calls/*', mode: 'copy', overwrite: true

    input:
    val ready
    path refseq_segments_dir
    val seq_mode
    tuple val(meta), path(reads), val(flu_type), val(subtype_ha), val(subtype_na)

    output:
    path("full_variant_calls/${meta.id}"), emit: sample_dirs
    path("reports/${meta.id}.fullvarcall.log"), emit: reports

    script:
    def readsArgs = reads.collect { "\"${it.getFileName().toString()}\"" }.join(' ')

    """
    mkdir -p "full_variant_calls/${meta.id}" reports

    python "${projectDir}/bin/run_full_variant_call_sample.py" \
        --sample-id "${meta.id}" \
        --flu-type "${flu_type}" \
        --subtype-ha "${subtype_ha}" \
        --subtype-na "${subtype_na}" \
        --seq-mode "${seq_mode}" \
        --refseq-dir "${refseq_segments_dir}" \
        --output-dir "full_variant_calls/${meta.id}" \
        --log-file "reports/${meta.id}.fullvarcall.log" \
        --ivar-freq ${params.ivar_freq} \
        --ivar-depth ${params.ivar_depth} \
        --threads ${task.cpus} \
        --reads ${readsArgs}
    """

    stub:
    """
    mkdir -p "full_variant_calls/${meta.id}" reports
    printf 'sample\ttype\tsegment\tgene\taa_position\tref_aa\talt_aa\tfrequency\tdepth\tmutation\n' > "full_variant_calls/${meta.id}/${meta.id}_protein_mutations.tsv"
    printf '${meta.id}\t${flu_type}\tHA\tHA\t1\tM\tI\t0.5000\t100\tM1I\n' >> "full_variant_calls/${meta.id}/${meta.id}_protein_mutations.tsv"
    echo 'stub full variant calls for ${meta.id}' > "reports/${meta.id}.fullvarcall.log"
    """
}
