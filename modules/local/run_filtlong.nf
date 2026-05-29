process RUN_FILTLONG {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/preprocessed_reads/filtlong", pattern: 'filtered/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/filtlong", pattern: 'reports/*.filtlong.stats.tsv', mode: 'copy', overwrite: true

    input:
    val ready
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('filtered/*'), emit: cleaned_reads
    tuple val(meta), path('reports/*'), emit: reports
    tuple val(meta), path("${meta.id}_filtlong_manifest.tsv"), emit: manifests

    script:
    def minLen = params.min_len_long as String
    def minQual = (params.filtlong_min_mean_q != null ? params.filtlong_min_mean_q : params.min_qual) as String
    def maxLen = params.max_len_long as Integer
    def maxLenArg = maxLen > 0 ? "--max_length ${maxLen}" : ""

    """
    mkdir -p filtered reports

    filtlong \
        --min_length ${minLen} \
        --min_mean_q ${minQual} \
        ${maxLenArg} \
        "${reads[0]}" \
        > >(gzip -c > "filtered/${meta.id}.filtered.fastq.gz") \
        2> "reports/${meta.id}.filtlong.log"

    python3 "${projectDir}/bin/fastq_stats.py" \
        --sample-id "${meta.id}" \
        --layout "${meta.layout}" \
        --seq-type "${meta.seq_type}" \
        --input-fastq "${reads[0]}" \
        --output-fastq "filtered/${meta.id}.filtered.fastq.gz" \
        --output-tsv "reports/${meta.id}.filtlong.stats.tsv" \
        --min-length "${minLen}" \
        --min-mean-q "${minQual}" \
        --max-length "${maxLen}"

    {
        printf 'sample_id\\tlayout\\tseq_type\\tread_count\\tread1\\tmin_length\\tmin_mean_q\\tmax_length\\tlog_file\\tstats_tsv\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
            "${meta.id}" \
            "${meta.layout}" \
            "${meta.seq_type}" \
            "1" \
            "filtered/${meta.id}.filtered.fastq.gz" \
            "${minLen}" \
            "${minQual}" \
            "${maxLen}" \
            "reports/${meta.id}.filtlong.log" \
            "reports/${meta.id}.filtlong.stats.tsv"
    } > "${meta.id}_filtlong_manifest.tsv"
    """

    stub:
    """
    mkdir -p filtered reports
    printf '@${meta.id}\\nACGTACGTACGT\\n+\\nIIIIIIIIIIII\\n' | gzip -c > "filtered/${meta.id}.filtered.fastq.gz"
    echo 'stub Filtlong log' > "reports/${meta.id}.filtlong.log"
    {
        printf 'sample_id\\tlayout\\tseq_type\\tinput_reads\\toutput_reads\\tread_retention_pct\\tinput_mean_len\\toutput_mean_len\\tmin_length\\tmin_mean_q\\tmax_length\\n'
        printf '%s\\t%s\\t%s\\t1\\t1\\t100.00\\t12\\t12\\t%s\\t%s\\t%s\\n' "${meta.id}" "${meta.layout}" "${meta.seq_type}" "${params.min_len_long}" "${params.filtlong_min_mean_q != null ? params.filtlong_min_mean_q : params.min_qual}" "${params.max_len_long}"
    } > "reports/${meta.id}.filtlong.stats.tsv"
    {
        printf 'sample_id\\tlayout\\tseq_type\\tread_count\\tread1\\tmin_length\\tmin_mean_q\\tmax_length\\tlog_file\\tstats_tsv\\n'
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "${meta.layout}" "${meta.seq_type}" "1" "filtered/${meta.id}.filtered.fastq.gz" "${params.min_len_long}" "${params.filtlong_min_mean_q != null ? params.filtlong_min_mean_q : params.min_qual}" "${params.max_len_long}" "reports/${meta.id}.filtlong.log" "reports/${meta.id}.filtlong.stats.tsv"
    } > "${meta.id}_filtlong_manifest.tsv"
    """
}

