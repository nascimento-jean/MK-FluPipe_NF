process RUN_MERGE_FULLVARCALL {
    tag 'merge-fullvarcall'
    label 'mk_flu_tools'
    publishDir "${params.output_dir}", pattern: 'full_variant_calls/*', mode: 'copy', overwrite: true

    input:
    path sample_dirs, stageAs: 'input??/*'

    output:
    path('full_variant_calls/all_samples_protein_mutations.tsv'), emit: summary

    script:
    """
    mkdir -p full_variant_calls
    printf 'sample	type	segment	gene	aa_position	ref_aa	alt_aa	frequency	depth	mutation
' > full_variant_calls/all_samples_protein_mutations.tsv

    count=0
    for sample_dir in input*/*; do
        for tsv in "\$sample_dir"/*_protein_mutations.tsv; do
            [ -f "\$tsv" ] || continue
            tail -n +2 "\$tsv" >> full_variant_calls/all_samples_protein_mutations.tsv
            count=\$((count + 1))
        done
    done

    printf 'Aggregated protein mutation files: %s
' "\$count"
    """

    stub:
    """
    mkdir -p full_variant_calls
    printf 'sample\ttype\tsegment\tgene\taa_position\tref_aa\talt_aa\tfrequency\tdepth\tmutation\n' > full_variant_calls/all_samples_protein_mutations.tsv
    for tsv in input*/*/*_protein_mutations.tsv; do
        [ -f "\$tsv" ] || continue
        tail -n +2 "\$tsv" >> full_variant_calls/all_samples_protein_mutations.tsv
    done
    """
}
