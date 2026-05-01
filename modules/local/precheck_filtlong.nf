process PRECHECK_FILTLONG {
    tag 'filtlong-check'
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    """
    if ! command -v filtlong >/dev/null 2>&1; then
        echo "filtlong is not available in PATH." >&2
        echo "Check the configured Conda environment in --mk_flu_conda_env or install filtlong manually." >&2
        exit 1
    fi
    """
}

