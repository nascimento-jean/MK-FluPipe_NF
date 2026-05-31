process RUN_SNPEFF_ANNOTATION {
    tag 'snpeff-annotation'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}", pattern: 'snpeff_annotation/*', mode: 'copy', overwrite: true

    input:
    path refseq_segments_dir
    path blast_summary
    path sample_dirs, stageAs: 'input??/*'

    output:
    path('snpeff_annotation/snpeff_annotation.tsv'), emit: report
    path('reports/snpeff_annotation.log'), emit: log

    script:
    """
    mkdir -p snpeff_annotation reports
    python3 "${projectDir}/bin/run_snpeff_annotation.py" \
        --sample-dirs input*/* \
        --refseq-dir "${refseq_segments_dir}" \
        --blast-summary "${blast_summary}" \
        --output-dir snpeff_annotation \
        --log-file reports/snpeff_annotation.log
    """

    stub:
    """
    mkdir -p snpeff_annotation reports
    printf 'sample\ttype\tsubtype\tsegment\tsegment_name\taccession\tpos\tref\talt\teffect\timpact\tgene\tfeature\thgvs_c\thgvs_p\tfrequency\tdepth\tstatus\tmessage\n' > snpeff_annotation/snpeff_annotation.tsv
    for sample_dir in input*/*; do
        sample=\$(basename "\$sample_dir")
        printf '%s\tA\tH3N2\t4\tHA\tNC_STUB\t1\tA\tG\tmissense_variant\tMODERATE\tHA\ttranscript\tc.1A>G\tp.Met1Ile\t0.5000\t100\tannotated\t-\n' "\$sample" >> snpeff_annotation/snpeff_annotation.tsv
    done
    echo 'stub SnpEff annotation' > reports/snpeff_annotation.log
    """
}
