process PRECHECK_IRMA_SHORT {
    tag 'irma-short-check'
    label 'irma_tools'

    output:
    val(true), emit: ok

    script:
    """
    if ! command -v IRMA >/dev/null 2>&1; then
        echo "IRMA is not available in PATH." >&2
        echo "Expose the IRMA binary in the execution environment before running the assembly stage." >&2
        exit 1
    fi
    """
}

