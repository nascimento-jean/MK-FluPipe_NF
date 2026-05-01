process PRECHECK_BOWTIE2_HOST_DEPLETION {
    tag 'bowtie2-host-check'
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    """
    for cmd in bowtie2 bowtie2-build; do
        if ! command -v "\${cmd}" >/dev/null 2>&1; then
            echo "\${cmd} is not available in PATH." >&2
            echo "Check the configured Conda environment in --mk_flu_conda_env." >&2
            exit 1
        fi
    done
    """
}

