process DISCOVER_SAMPLES {
    tag "${input_dir}"
    publishDir "${params.output_dir}/bootstrap", mode: 'copy', overwrite: true

    input:
    val input_dir
    val seq_type

    output:
    path 'samplesheet.csv', emit: samplesheet
    path 'sample_summary.json', emit: summary

    script:
    """
    python3 "${projectDir}/bin/discover_samples.py" \
        --input-dir "${input_dir}" \
        --seq-type "${seq_type}" \
        --samplesheet samplesheet.csv \
        --summary sample_summary.json
    """
}

