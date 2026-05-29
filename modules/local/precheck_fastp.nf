process PRECHECK_FASTP {
    tag 'fastp-check'
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    """
    if ! command -v fastp >/dev/null 2>&1; then
        echo "fastp is not available in PATH." >&2
        echo "Check the configured Conda environment in --mk_flu_conda_env or install fastp manually." >&2
        exit 1
    fi
    """

    stub:
    """
    true
    """
}
