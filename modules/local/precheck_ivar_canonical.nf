process PRECHECK_IVAR_CANONICAL {
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    """
    ivar version >/dev/null 2>&1 || ivar -h >/dev/null
    minimap2 --version >/dev/null
    samtools --version >/dev/null
    """
}
