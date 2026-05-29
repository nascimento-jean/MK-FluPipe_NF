process RUN_PHYLOGENY {
    tag 'augur-ha-na'
    label 'mk_flu_tools'

    input:
    path blast_summary
    path segment_files, stageAs: 'segments/*'
    path metadata_csv
    path context_fasta
    path context_metadata
    val min_sequences

    output:
    path 'phylogeny', emit: results
    path 'reports/phylogeny.log', emit: log

    script:
    """
    mkdir -p phylogeny reports

    # Preserve dated trees, separated state/source colors, and offline HTML tree views.
    python3 "${projectDir}/bin/run_phylogeny.py" \
        --blast-summary "${blast_summary}" \
        --metadata "${metadata_csv}" \
        --context-fasta "${context_fasta}" \
        --context-metadata "${context_metadata}" \
        --output-dir phylogeny \
        --log-file reports/phylogeny.log \
        --threads ${task.cpus} \
        --min-sequences "${min_sequences}" \
        segments/*
    """
}
