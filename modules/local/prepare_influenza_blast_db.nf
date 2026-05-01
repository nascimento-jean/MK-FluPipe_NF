process PREPARE_INFLUENZA_BLAST_DB {
    tag 'influenza-blast-db'
    label 'mk_flu_tools'

    output:
    path 'blast_db', emit: blast_db_dir

    script:
    def blastDbDir = params.blast_db_dir as String
    def blastDbFasta = params.blast_db_fasta as String
    def blastDbPrefix = params.blast_db_prefix as String
    def blastDbTimestamp = params.blast_db_timestamp as String
    def blastDbUrl = params.blast_db_url as String
    def blastDbMaxDays = params.blast_db_max_days as Integer

    """
    mkdir -p blast_db
    mkdir -p "${blastDbDir}"

    DB_DIR="${blastDbDir}"
    DB_FILE="${blastDbFasta}"
    DB_PREFIX="${blastDbPrefix}"
    DB_TIMESTAMP="${blastDbTimestamp}"
    DB_URL="${blastDbUrl}"
    DB_MAX_DAYS=${blastDbMaxDays}

    days_since_modified() {
        local file="\$1"
        local now mod
        now=\$(date +%s)
        mod=\$(date -r "\$file" +%s 2>/dev/null || echo 0)
        echo \$(( (now - mod) / 86400 ))
    }

    download_db() {
        local gz_file="\${DB_DIR}/influenza.fna.gz"
        if command -v wget >/dev/null 2>&1; then
            wget -q --show-progress -O "\$gz_file" "\${DB_URL}" || return 1
        elif command -v curl >/dev/null 2>&1; then
            curl -L "\${DB_URL}" -o "\$gz_file" || return 1
        else
            echo "Neither wget nor curl is available to download the influenza database." >&2
            return 1
        fi
        gunzip -f "\$gz_file"
        makeblastdb \
            -in "\${DB_FILE}" \
            -dbtype nucl \
            -out "\${DB_PREFIX}" \
            -title "NCBI_Influenza" \
            -parse_seqids
        date '+%Y-%m-%d' > "\${DB_TIMESTAMP}"
    }

    if [[ ! -f "\${DB_FILE}" || ! -f "\${DB_PREFIX}.nhr" ]]; then
        download_db
    else
        days_old=0
        if [[ -f "\${DB_TIMESTAMP}" ]]; then
            days_old=\$(days_since_modified "\${DB_TIMESTAMP}")
        else
            days_old=\$((DB_MAX_DAYS + 1))
        fi

        if [[ "\$days_old" -gt "\${DB_MAX_DAYS}" ]]; then
            if ! download_db; then
                echo "Database refresh failed; continuing with existing database." >&2
            fi
        fi
    fi

    ln -sf "\${DB_FILE}" "blast_db/influenza.fna"
    for idx in "\${DB_PREFIX}"*; do
        [ -f "\$idx" ] || continue
        ln -sf "\$idx" "blast_db/\$(basename "\$idx")"
    done
    """
}
