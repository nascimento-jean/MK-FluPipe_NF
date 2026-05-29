process PRECHECK_SAMTOOLS_DEPTH {
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    """
    samtools --version >/dev/null
    """

    stub:
    """
    true
    """
}
