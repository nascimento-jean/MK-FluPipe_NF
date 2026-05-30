process PRECHECK_FULLVARCALL {
    tag 'fullvarcall-check'
    label 'mk_flu_tools'

    output:
    val(true), emit: ok

    script:
    def seqMode = params.seq_type as String

    """
    required_tools=(ivar samtools)

    if [[ "${seqMode}" == "long" ]]; then
        required_tools+=(minimap2)
    else
        required_tools+=(bowtie2 bowtie2-build)
    fi

    for tool in "\${required_tools[@]}"; do
        command -v "\$tool" >/dev/null 2>&1 || {
            echo "Required tool not found in PATH: \$tool" >&2
            exit 1
        }
    done
    """

    stub:
    """
    true
    """
}
