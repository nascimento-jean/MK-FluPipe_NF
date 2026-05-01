process PREPARE_HUMAN_BOWTIE2_INDEX {
    tag 'human-reference'
    label 'mk_flu_tools'

    output:
    path 'human/GRCh38_no_alt.fna', emit: human_fasta
    path 'human/GRCh38*.bt2*', emit: index_files

    script:
    def humanDbDir = params.human_db_dir as String
    def humanFastaName = params.human_fasta_name as String
    def humanIndexPrefix = params.human_index_prefix as String
    def humanGenomeUrl = params.human_genome_url as String

    """
    mkdir -p human
    mkdir -p "${humanDbDir}"

    DB_DIR="${humanDbDir}"
    HUMAN_FA_NAME="${humanFastaName}"
    HUMAN_INDEX_PREFIX="${humanIndexPrefix}"
    HUMAN_FA="\${DB_DIR}/\${HUMAN_FA_NAME}"
    HUMAN_FA_GZ="\${HUMAN_FA}.gz"
    HUMAN_INDEX="\${DB_DIR}/\${HUMAN_INDEX_PREFIX}"
    HUMAN_URL="${humanGenomeUrl}"

    if [[ ! -f "\${HUMAN_FA}" ]]; then
        if command -v wget >/dev/null 2>&1; then
            wget -q --show-progress -O "\${HUMAN_FA_GZ}" "\${HUMAN_URL}"
        elif command -v curl >/dev/null 2>&1; then
            curl -L "\${HUMAN_URL}" -o "\${HUMAN_FA_GZ}"
        else
            echo "Neither wget nor curl is available to download the human genome." >&2
            exit 1
        fi
        gunzip -f "\${HUMAN_FA_GZ}"
    fi

    if [[ ! -f "\${HUMAN_INDEX}.1.bt2" && ! -f "\${HUMAN_INDEX}.1.bt2l" ]]; then
        bowtie2-build \
            --threads ${task.cpus} \
            "\${HUMAN_FA}" \
            "\${HUMAN_INDEX}"
    fi

    ln -sf "\${HUMAN_FA}" "human/\${HUMAN_FA_NAME}"
    for idx in "\${HUMAN_INDEX}"*.bt2*; do
        ln -sf "\${idx}" "human/\$(basename "\${idx}")"
    done
    """
}

