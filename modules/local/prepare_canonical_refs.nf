process PREPARE_CANONICAL_REFS {
    tag 'canonical-refs'
    label 'mk_flu_tools'

    output:
    path 'canonical_refs', emit: refs_dir
    path 'canonical_refs.log', emit: log

    script:
    def canonicalRefsDir = params.canonical_refs_dir as String

    """
    mkdir -p "${canonicalRefsDir}"
    python "${projectDir}/bin/prepare_canonical_refs.py" \
        --output-dir "${canonicalRefsDir}" \
        --log-file canonical_refs.log

    mkdir -p canonical_refs
    for ref in "${canonicalRefsDir}"/*.fa; do
        [ -f "\$ref" ] || continue
        ln -sf "\$ref" "canonical_refs/\$(basename "\$ref")"
    done
    """

    stub:
    """
    mkdir -p canonical_refs
    for ref in N1_NA N2_NA N1_PA N2_PA N1_MP N2_MP; do
        printf '>%s\\nACGTACGTACGT\\n' "\$ref" > "canonical_refs/\${ref}.fa"
    done
    echo 'stub canonical refs' > canonical_refs.log
    """
}
