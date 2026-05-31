process RUN_FUNCTIONAL_ANNOTATION {
    tag 'functional-annotation'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}", pattern: 'functional_annotation/*', mode: 'copy', overwrite: true

    input:
    path fullvar_summary
    path typing_summary

    output:
    path('functional_annotation/functional_annotation.tsv'), emit: report

    script:
    """
    mkdir -p functional_annotation
    python3 "${projectDir}/bin/annotate_functional_variants.py" \
        --fullvar-tsv "${fullvar_summary}" \
        --typing-tsv "${typing_summary}" \
        --output-tsv functional_annotation/functional_annotation.tsv
    """
}
