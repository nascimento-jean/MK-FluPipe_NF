process RUN_HOST_DEPLETION_BOWTIE2 {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/depleted_reads/bowtie2", pattern: 'depleted/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/host_depletion_bowtie2", pattern: 'reports/*.host_depletion.stats.tsv', mode: 'copy', overwrite: true

    input:
    val ready
    path human_fasta
    path human_index_files
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('depleted/*'), emit: depleted_reads
    tuple val(meta), path('reports/*'), emit: reports
    tuple val(meta), path("${meta.id}_host_depletion_manifest.tsv"), emit: manifests

    script:
    def indexPrefix = "human/${params.human_index_prefix as String}"

    if( meta.layout == 'paired' ) {
        """
        mkdir -p depleted reports human
        cp -L "${human_fasta}" human/
        cp -L ${human_index_files.collect { "\"${it}\"" }.join(' ')} human/

        bowtie2 \
            -x "${indexPrefix}" \
            -1 "${reads[0]}" \
            -2 "${reads[1]}" \
            --un-conc-gz "depleted/${meta.id}_depleted_R%.fastq.gz" \
            -p ${task.cpus} \
            --very-sensitive \
            --no-unal \
            -S /dev/null \
            2> "reports/${meta.id}.bowtie2_host_depletion.log" || true

        [ -f "depleted/${meta.id}_depleted_R_1.fastq.gz" ] && mv "depleted/${meta.id}_depleted_R_1.fastq.gz" "depleted/${meta.id}_depleted_R1.fastq.gz" || true
        [ -f "depleted/${meta.id}_depleted_R_2.fastq.gz" ] && mv "depleted/${meta.id}_depleted_R_2.fastq.gz" "depleted/${meta.id}_depleted_R2.fastq.gz" || true

        # Original pipeline behavior: if depletion yields no usable reads, fall back to trimmed/raw reads.
        [ -s "depleted/${meta.id}_depleted_R1.fastq.gz" ] || {
            cp "${reads[0]}" "depleted/${meta.id}_depleted_R1.fastq.gz"
            cp "${reads[1]}" "depleted/${meta.id}_depleted_R2.fastq.gz"
        }

        python3 "${projectDir}/bin/read_set_stats.py" \
            --sample-id "${meta.id}" \
            --layout "${meta.layout}" \
            --seq-type "${meta.seq_type}" \
            --input-fastq "${reads[0]}" "${reads[1]}" \
            --output-fastq "depleted/${meta.id}_depleted_R1.fastq.gz" "depleted/${meta.id}_depleted_R2.fastq.gz" \
            --output-tsv "reports/${meta.id}.host_depletion.stats.tsv"

        {
            printf 'sample_id\tlayout\tseq_type\tread1\tread2\tlog_file\tstats_tsv\n'
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${meta.id}" \
                "${meta.layout}" \
                "${meta.seq_type}" \
                "depleted/${meta.id}_depleted_R1.fastq.gz" \
                "depleted/${meta.id}_depleted_R2.fastq.gz" \
                "reports/${meta.id}.bowtie2_host_depletion.log" \
                "reports/${meta.id}.host_depletion.stats.tsv"
        } > "${meta.id}_host_depletion_manifest.tsv"
        """
    }
    else {
        """
        mkdir -p depleted reports human
        cp -L "${human_fasta}" human/
        cp -L ${human_index_files.collect { "\"${it}\"" }.join(' ')} human/

        bowtie2 \
            -x "${indexPrefix}" \
            -U "${reads[0]}" \
            --un-gz "depleted/${meta.id}_depleted.fastq.gz" \
            -p ${task.cpus} \
            --very-sensitive \
            --no-unal \
            -S /dev/null \
            2> "reports/${meta.id}.bowtie2_host_depletion.log" || true

        # Original pipeline behavior: if depletion yields no usable reads, fall back to trimmed/raw reads.
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
                "reports/${meta.id}.bowtie2_host_depletion.log" \
                "reports/${meta.id}.host_depletion.stats.tsv"
        } > "${meta.id}_host_depletion_manifest.tsv"
        """
    }
}
