process PRECHECK_FASTQC {
    tag 'fastqc-check'
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    """
    if ! command -v fastqc >/dev/null 2>&1; then
        echo "FastQC is not available in PATH." >&2
        echo "Check the configured Conda environment in --mk_flu_conda_env or install FastQC manually." >&2
        exit 1
    fi
    """

    stub:
    """
    true
    """
}
