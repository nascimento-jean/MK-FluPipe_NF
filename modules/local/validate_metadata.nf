process VALIDATE_METADATA {
    tag 'metadata'
    publishDir "${params.output_dir}/bootstrap", mode: 'copy', overwrite: true

    input:
    path samplesheet
    path metadata_csv

    output:
    path 'validated_samplesheet.csv', emit: samplesheet
    path 'validated_metadata.csv', emit: metadata
    path 'metadata_validation.json', emit: summary

    script:
    """
    python3 "${projectDir}/bin/validate_metadata.py" \
        --metadata "${metadata_csv}" \
        --samplesheet "${samplesheet}" \
        --output validated_metadata.csv \
        --summary metadata_validation.json \
        --validated-samplesheet validated_samplesheet.csv
    """
}
