process PRECHECK_MEDAKA_VARIANTS {
    label 'medaka_tools'

    output:
    val(true), emit: ok

    script:
    """
    medaka_haploid_variant -h >/dev/null 2>&1 || medaka_variant -h >/dev/null 2>&1
    samtools --version >/dev/null
    """
}
