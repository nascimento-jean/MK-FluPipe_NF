process PRECHECK_MULTIQC {
    tag 'multiqc-check'
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    """
    if ! command -v multiqc >/dev/null 2>&1; then
        echo "MultiQC is not available in PATH." >&2
        exit 1
    fi
    """
}
