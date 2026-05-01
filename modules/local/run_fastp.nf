process RUN_FASTP {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/preprocessed_reads/fastp", pattern: 'trimmed/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/fastp", pattern: 'reports/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/fastp", pattern: '*_fastp_manifest.tsv', mode: 'copy', overwrite: true

    input:
    val ready
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('trimmed/*'), emit: cleaned_reads
    tuple val(meta), path('reports/*'), emit: reports
    tuple val(meta), path("${meta.id}_fastp_manifest.tsv"), emit: manifests

    script:
    def minLen = params.min_len_short as String
    def minQual = params.min_qual as String
    def adapterArg = params.adapter_fasta ? "--adapter_fasta \"${params.adapter_fasta}\"" : ""

    if( meta.layout == 'paired' ) {
        """
        mkdir -p trimmed reports

        fastp \
            --thread ${task.cpus} \
            --in1 "${reads[0]}" \
            --in2 "${reads[1]}" \
            --out1 "trimmed/${meta.id}_R1.trimmed.fastq.gz" \
            --out2 "trimmed/${meta.id}_R2.trimmed.fastq.gz" \
            --html "reports/${meta.id}.fastp.html" \
            --json "reports/${meta.id}.fastp.json" \
            --qualified_quality_phred ${minQual} \
            --length_required ${minLen} ${adapterArg}

        {
            printf 'sample_id\\tlayout\\tseq_type\\tread_count\\tread1\\tread2\\thtml\\tjson\\n'
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
                "${meta.id}" \
                "${meta.layout}" \
                "${meta.seq_type}" \
                "2" \
                "trimmed/${meta.id}_R1.trimmed.fastq.gz" \
                "trimmed/${meta.id}_R2.trimmed.fastq.gz" \
                "reports/${meta.id}.fastp.html" \
                "reports/${meta.id}.fastp.json"
        } > "${meta.id}_fastp_manifest.tsv"
        """
    }
    else {
        """
        mkdir -p trimmed reports

        fastp \
            --thread ${task.cpus} \
            --in1 "${reads[0]}" \
            --out1 "trimmed/${meta.id}.trimmed.fastq.gz" \
            --html "reports/${meta.id}.fastp.html" \
            --json "reports/${meta.id}.fastp.json" \
            --qualified_quality_phred ${minQual} \
            --length_required ${minLen} ${adapterArg}

        {
            printf 'sample_id\\tlayout\\tseq_type\\tread_count\\tread1\\tread2\\thtml\\tjson\\n'
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
                "${meta.id}" \
                "${meta.layout}" \
                "${meta.seq_type}" \
                "1" \
                "trimmed/${meta.id}.trimmed.fastq.gz" \
                "-" \
                "reports/${meta.id}.fastp.html" \
                "reports/${meta.id}.fastp.json"
        } > "${meta.id}_fastp_manifest.tsv"
        """
    }
}
