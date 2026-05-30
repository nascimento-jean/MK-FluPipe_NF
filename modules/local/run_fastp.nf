process RUN_FASTP {
    tag { meta.id }
    label 'mk_flu_tools'
    publishDir "${params.output_dir}/preprocessed_reads/fastp", pattern: 'trimmed/*', mode: 'copy', overwrite: true
    publishDir "${params.output_dir}/qc_reports/fastp", pattern: 'reports/*.fastp.*', mode: 'copy', overwrite: true

    input:
    val ready
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('trimmed/*'), emit: cleaned_reads
    tuple val(meta), path('reports/*'), emit: reports
    tuple val(meta), path("${meta.id}_fastp_manifest.tsv"), emit: manifests

    script:
    def minLen = params.min_len_short as String
    def minQual = params.min_qual as String
    def adapterArg = params.adapter_fasta ? "--adapter_fasta \"${params.adapter_fasta}\"" : ""
    def fastpThreads = task.cpus as int
    def pigzThreads  = Math.max(1, fastpThreads)
    // Hard timeout for fastp itself (in seconds). If fastp deadlocks on its
    // internal pthread futex (a known issue: see fastp v1.3.2/v1.3.3 release
    // notes "Fix possible hang"), this forces it to die so Nextflow can retry
    // the task instead of hanging forever. Configurable via params.fastp_timeout
    // in nextflow.config (default 1800s = 30 min).
    def fastpTimeoutSec = (params.fastp_timeout ?: 1800) as int
    // Optional startup watchdog. It is disabled by default because fastp can
    // spend several minutes detecting adapters/filtering before creating final
    // output files, especially for larger samples.
    def fastpStartupTimeoutSec = (params.fastp_startup_timeout ?: 0) as int

    if( meta.layout == 'paired' ) {
        """
        set -euo pipefail
        mkdir -p trimmed reports

        run_fastp_with_startup_watch() {
            "\$@" &
            fastp_pid=\$!
            waited=0

            while kill -0 "\$fastp_pid" 2>/dev/null; do
                for expected in "\${expected_outputs[@]}"; do
                    if [ -s "\$expected" ]; then
                        wait "\$fastp_pid"
                        return \$?
                    fi
                done

                if [ ${fastpStartupTimeoutSec} -gt 0 ] && [ "\$waited" -ge ${fastpStartupTimeoutSec} ]; then
                    echo "fastp produced no output for ${fastpStartupTimeoutSec}s; killing stalled process." >&2
                    pkill -TERM -P "\$fastp_pid" 2>/dev/null || true
                    kill -TERM "\$fastp_pid" 2>/dev/null || true
                    sleep 5
                    pkill -KILL -P "\$fastp_pid" 2>/dev/null || true
                    kill -KILL "\$fastp_pid" 2>/dev/null || true
                    wait "\$fastp_pid" 2>/dev/null || true
                    exit 124
                fi

                sleep 5
                waited=\$((waited + 5))
            done

            wait "\$fastp_pid"
        }

        # ---------------------------------------------------------------------
        # WHY UNCOMPRESSED OUTPUT + EXTERNAL pigz:
        # fastp's internal multi-threaded gzip writer (libdeflate) can deadlock
        # on a pthread condition variable under heavy parallel I/O inside
        # Singularity containers. The deadlock is silent (0% CPU, futex_wait
        # forever). Writing plain .fastq and compressing with a separate pigz
        # process eliminates the cross-thread synchronisation that triggers it.
        #
        # WHY `timeout`:
        # Even with the above mitigation, if fastp ever does deadlock again,
        # this kills it after fastpTimeoutSec seconds with exit 124, letting
        # Nextflow's retry strategy (in nextflow.config) handle it cleanly
        # instead of hanging the whole pipeline forever.
        # ---------------------------------------------------------------------
        expected_outputs=(
            "trimmed/${meta.id}_R1.trimmed.fastq"
            "trimmed/${meta.id}_R2.trimmed.fastq"
            "reports/${meta.id}.fastp.html"
            "reports/${meta.id}.fastp.json"
        )

        run_fastp_with_startup_watch timeout --kill-after=30s ${fastpTimeoutSec}s fastp \\
            --thread ${fastpThreads} \\
            --in1 "${reads[0]}" \\
            --in2 "${reads[1]}" \\
            --out1 "trimmed/${meta.id}_R1.trimmed.fastq" \\
            --out2 "trimmed/${meta.id}_R2.trimmed.fastq" \\
            --html "reports/${meta.id}.fastp.html" \\
            --json "reports/${meta.id}.fastp.json" \\
            --qualified_quality_phred ${minQual} \\
            --length_required ${minLen} \\
            --detect_adapter_for_pe ${adapterArg}

        # Ensure fastp's writes are flushed to the filesystem before pigz starts.
        sync

        # Compress in a separate, single-purpose process. pigz cleanly returns
        # an exit code; no deadlock potential.
        pigz -p ${pigzThreads} -f "trimmed/${meta.id}_R1.trimmed.fastq"
        pigz -p ${pigzThreads} -f "trimmed/${meta.id}_R2.trimmed.fastq"

        {
            printf 'sample_id\\tlayout\\tseq_type\\tread_count\\tread1\\tread2\\thtml\\tjson\\n'
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
                "${meta.id}" \\
                "${meta.layout}" \\
                "${meta.seq_type}" \\
                "2" \\
                "trimmed/${meta.id}_R1.trimmed.fastq.gz" \\
                "trimmed/${meta.id}_R2.trimmed.fastq.gz" \\
                "reports/${meta.id}.fastp.html" \\
                "reports/${meta.id}.fastp.json"
        } > "${meta.id}_fastp_manifest.tsv"
        """
    }
    else {
        """
        set -euo pipefail
        mkdir -p trimmed reports

        run_fastp_with_startup_watch() {
            "\$@" &
            fastp_pid=\$!
            waited=0

            while kill -0 "\$fastp_pid" 2>/dev/null; do
                for expected in "\${expected_outputs[@]}"; do
                    if [ -s "\$expected" ]; then
                        wait "\$fastp_pid"
                        return \$?
                    fi
                done

                if [ ${fastpStartupTimeoutSec} -gt 0 ] && [ "\$waited" -ge ${fastpStartupTimeoutSec} ]; then
                    echo "fastp produced no output for ${fastpStartupTimeoutSec}s; killing stalled process." >&2
                    pkill -TERM -P "\$fastp_pid" 2>/dev/null || true
                    kill -TERM "\$fastp_pid" 2>/dev/null || true
                    sleep 5
                    pkill -KILL -P "\$fastp_pid" 2>/dev/null || true
                    kill -KILL "\$fastp_pid" 2>/dev/null || true
                    wait "\$fastp_pid" 2>/dev/null || true
                    exit 124
                fi

                sleep 5
                waited=\$((waited + 5))
            done

            wait "\$fastp_pid"
        }

        # See paired-end branch above for full rationale on uncompressed
        # output, external pigz, and the timeout wrapper.
        expected_outputs=(
            "trimmed/${meta.id}.trimmed.fastq"
            "reports/${meta.id}.fastp.html"
            "reports/${meta.id}.fastp.json"
        )

        run_fastp_with_startup_watch timeout --kill-after=30s ${fastpTimeoutSec}s fastp \\
            --thread ${fastpThreads} \\
            --in1 "${reads[0]}" \\
            --out1 "trimmed/${meta.id}.trimmed.fastq" \\
            --html "reports/${meta.id}.fastp.html" \\
            --json "reports/${meta.id}.fastp.json" \\
            --qualified_quality_phred ${minQual} \\
            --length_required ${minLen} ${adapterArg}

        sync

        pigz -p ${pigzThreads} -f "trimmed/${meta.id}.trimmed.fastq"

        {
            printf 'sample_id\\tlayout\\tseq_type\\tread_count\\tread1\\tread2\\thtml\\tjson\\n'
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
                "${meta.id}" \\
                "${meta.layout}" \\
                "${meta.seq_type}" \\
                "1" \\
                "trimmed/${meta.id}.trimmed.fastq.gz" \\
                "-" \\
                "reports/${meta.id}.fastp.html" \\
                "reports/${meta.id}.fastp.json"
        } > "${meta.id}_fastp_manifest.tsv"
        """
    }

    stub:
    if( meta.layout == 'paired' ) {
        """
        mkdir -p trimmed reports
        printf '@${meta.id}/1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > "trimmed/${meta.id}_R1.trimmed.fastq.gz"
        printf '@${meta.id}/2\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > "trimmed/${meta.id}_R2.trimmed.fastq.gz"
        echo '<html><body>stub fastp</body></html>' > "reports/${meta.id}.fastp.html"
        printf '{"summary":{"before_filtering":{"total_reads":2,"total_bases":16},"after_filtering":{"total_reads":2,"q30_rate":1.0,"read1_mean_length":8}},"filtering_result":{"passed_filter_reads":2}}\\n' > "reports/${meta.id}.fastp.json"
        {
            printf 'sample_id\\tlayout\\tseq_type\\tread_count\\tread1\\tread2\\thtml\\tjson\\n'
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "${meta.layout}" "${meta.seq_type}" "2" "trimmed/${meta.id}_R1.trimmed.fastq.gz" "trimmed/${meta.id}_R2.trimmed.fastq.gz" "reports/${meta.id}.fastp.html" "reports/${meta.id}.fastp.json"
        } > "${meta.id}_fastp_manifest.tsv"
        """
    }
    else {
        """
        mkdir -p trimmed reports
        printf '@${meta.id}\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > "trimmed/${meta.id}.trimmed.fastq.gz"
        echo '<html><body>stub fastp</body></html>' > "reports/${meta.id}.fastp.html"
        printf '{"summary":{"before_filtering":{"total_reads":1,"total_bases":8},"after_filtering":{"total_reads":1,"q30_rate":1.0,"read1_mean_length":8}},"filtering_result":{"passed_filter_reads":1}}\\n' > "reports/${meta.id}.fastp.json"
        {
            printf 'sample_id\\tlayout\\tseq_type\\tread_count\\tread1\\tread2\\thtml\\tjson\\n'
            printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "${meta.layout}" "${meta.seq_type}" "1" "trimmed/${meta.id}.trimmed.fastq.gz" "-" "reports/${meta.id}.fastp.html" "reports/${meta.id}.fastp.json"
        } > "${meta.id}_fastp_manifest.tsv"
        """
    }
}
