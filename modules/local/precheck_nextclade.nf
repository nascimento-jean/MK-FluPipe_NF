process PRECHECK_NEXTCLADE {
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    """
    nextclade --version >/dev/null
    """
}
