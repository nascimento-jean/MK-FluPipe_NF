process RUN_LEGACY_PIPELINE {
    tag 'legacy-bridge'
    publishDir "${params.output_dir}/legacy_bridge", mode: 'copy', overwrite: true

    input:
    val input_dir
    val output_dir
    val irma_module
    val seq_mode
    val seq_type
    val legacy_script
    val run_fastqc
    val host_depletion
    val run_ivar
    val run_medaka
    val run_antiviral
    val run_h5_virulence
    val run_fullvarcall
    val adapter_fasta
    val min_len_short
    val min_len_long
    val max_len_long
    val min_qual
    val min_coverage
    val max_n_pct
    val min_segments
    val ivar_freq
    val ivar_depth
    val minority_freq
    val coinfection_pct
    val medaka_env
    val gisaid_location
    val gisaid_year

    output:
    path 'legacy_bridge.done', emit: done_flag

    script:
    def seqModeArg = seq_mode ? "\"${seq_mode}\"" : "\"\""
    def adapterArg = adapter_fasta ? "--adapter_fasta \"${adapter_fasta}\"" : ""
    def gisaidLocationArg = gisaid_location ? "--gisaid_location \"${gisaid_location}\"" : ""
    def gisaidYearArg = gisaid_year ? "--gisaid_year \"${gisaid_year}\"" : ""

    """
    if [[ ! -f "${legacy_script}" ]]; then
        echo "Legacy script not found: ${legacy_script}" >&2
        exit 1
    fi

    bash "${legacy_script}" \
        "${input_dir}" \
        "${output_dir}" \
        "${irma_module}" \
        ${seqModeArg} \
        --seq_type "${seq_type}" \
        --run_fastqc "${run_fastqc}" \
        --min_len_short "${min_len_short}" \
        --min_len_long "${min_len_long}" \
        --max_len_long "${max_len_long}" \
        --min_qual "${min_qual}" \
        --host_depletion "${host_depletion}" \
        --min_coverage "${min_coverage}" \
        --max_n_pct "${max_n_pct}" \
        --min_segments "${min_segments}" \
        --ivar "${run_ivar}" \
        --medaka "${run_medaka}" \
        --medaka_env "${medaka_env}" \
        --ivar_freq "${ivar_freq}" \
        --ivar_depth "${ivar_depth}" \
        --minority_freq "${minority_freq}" \
        --coinfection_pct "${coinfection_pct}" \
        --antiviral "${run_antiviral}" \
        --h5_virulence "${run_h5_virulence}" \
        --fullvarcall "${run_fullvarcall}" \
        ${adapterArg} \
        ${gisaidLocationArg} \
        ${gisaidYearArg}

    touch legacy_bridge.done
    """
}

