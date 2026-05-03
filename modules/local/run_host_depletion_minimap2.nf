process RUN_HOST_DEPLETION_MINIMAP2 {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/depleted_reads/minimap2", pattern: 'depleted/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/host_depletion_minimap2", pattern: 'reports/*.host_depletion.stats.tsv', mode: 'copy', overwrite: true

    input:
    val ready
    path human_fasta
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('depleted/*'), emit: depleted_reads
    tuple val(meta), path('reports/*'), emit: reports
    tuple val(meta), path("${meta.id}_host_depletion_manifest.tsv"), emit: manifests

    script:
    """
    mkdir -p depleted reports human
    cp -L "${human_fasta}" human/

    minimap2 \
        -ax map-ont \
        -t ${task.cpus} \
        "human/${human_fasta.getName()}" \
        "${reads[0]}" \
        2> "reports/${meta.id}.minimap2_host_depletion.log" \
    | samtools view -f 4 -bS - \
    | samtools fastq - \
    | gzip -c > "depleted/${meta.id}_depleted.fastq.gz"

    # Original pipeline behavior: if depletion yields no usable reads, fall back to filtered/raw reads.
    [ -s "depleted/${meta.id}_depleted.fastq.gz" ] || cp "${reads[0]}" "depleted/${meta.id}_depleted.fastq.gz"

    python3 "${projectDir}/bin/read_set_stats.py" \
        --sample-id "${meta.id}" \
        --layout "${meta.layout}" \
        --seq-type "${meta.seq_type}" \
        --input-fastq "${reads[0]}" \
        --output-fastq "depleted/${meta.id}_depleted.fastq.gz" \
        --output-tsv "reports/${meta.id}.host_depletion.stats.tsv"

    {
        printf 'sample_id\tlayout\tseq_type\tread1\tread2\tlog_file\tstats_tsv\n'
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${meta.id}" \
            "${meta.layout}" \
            "${meta.seq_type}" \
            "depleted/${meta.id}_depleted.fastq.gz" \
            "-" \
            "reports/${meta.id}.minimap2_host_depletion.log" \
            "reports/${meta.id}.host_depletion.stats.tsv"
    } > "${meta.id}_host_depletion_manifest.tsv"
    """
}
