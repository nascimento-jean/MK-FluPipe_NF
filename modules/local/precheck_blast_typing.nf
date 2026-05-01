process PRECHECK_BLAST_TYPING {
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    """
    blastn -version >/dev/null
    makeblastdb -version >/dev/null
    """
}
