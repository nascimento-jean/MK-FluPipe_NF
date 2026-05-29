process PREPARE_ANTIVIRAL_DB {
    tag 'antiviral-db'
    label 'mk_flu_tools'
    publishDir "${params.blast_db_dir}/antiviral_resistance", pattern: 'antiviral_resistance/*', mode: 'copy', overwrite: true

    output:
    path('antiviral_resistance/flu_antiviral_markers.tsv'), emit: db
    path('antiviral_resistance/prepare_antiviral_db.log'), emit: log

    script:
    """
    mkdir -p antiviral_resistance

    python "${projectDir}/bin/prepare_antiviral_db.py" \
        --output-dir antiviral_resistance \
        --log-file antiviral_resistance/prepare_antiviral_db.log
    """

    stub:
    """
    mkdir -p antiviral_resistance
    printf 'gene\\taa_position\\twt_who\\tmut_who\\tdrug\\tsignificance\\tnomenclature\\n' > antiviral_resistance/flu_antiviral_markers.tsv
    echo 'stub antiviral db' > antiviral_resistance/prepare_antiviral_db.log
    """
}
