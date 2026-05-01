process PLAN_RUN {
    tag 'plan-run'
    publishDir "${params.output_dir}/bootstrap", mode: 'copy', overwrite: true

    input:
    path samplesheet
    path summary_json
    val seq_type
    val irma_module
    val run_fastqc
    val host_depletion
    val run_ivar
    val run_medaka
    val run_antiviral
    val run_h5_virulence
    val run_fullvarcall
    val run_legacy_bridge

    output:
    path 'execution_plan.txt', emit: execution_plan

    script:
    """
    python3 "${projectDir}/bin/render_plan.py" \
        --samplesheet "${samplesheet}" \
        --summary "${summary_json}" \
        --output execution_plan.txt \
        --seq-type "${seq_type}" \
        --irma-module "${irma_module}" \
        --run-fastqc "${run_fastqc}" \
        --host-depletion "${host_depletion}" \
        --run-ivar "${run_ivar}" \
        --run-medaka "${run_medaka}" \
        --run-antiviral "${run_antiviral}" \
        --run-h5-virulence "${run_h5_virulence}" \
        --run-fullvarcall "${run_fullvarcall}" \
        --run-legacy-bridge "${run_legacy_bridge}"
    """
}

