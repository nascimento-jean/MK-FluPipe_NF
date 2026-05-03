nextflow.enable.dsl = 2

include { DISCOVER_SAMPLES } from './modules/local/discover_samples'
include { PLAN_RUN } from './modules/local/plan_run'
include { PRECHECK_FASTQC } from './modules/local/precheck_fastqc'
include { RUN_FASTQC } from './modules/local/run_fastqc'
include { PRECHECK_FASTP } from './modules/local/precheck_fastp'
include { RUN_FASTP } from './modules/local/run_fastp'
include { PRECHECK_FILTLONG } from './modules/local/precheck_filtlong'
include { RUN_FILTLONG } from './modules/local/run_filtlong'
include { PREPARE_HUMAN_BOWTIE2_INDEX } from './modules/local/prepare_human_bowtie2_index'
include { PRECHECK_BOWTIE2_HOST_DEPLETION } from './modules/local/precheck_bowtie2_host_depletion'
include { RUN_HOST_DEPLETION_BOWTIE2 } from './modules/local/run_host_depletion_bowtie2'
include { PRECHECK_MINIMAP2_HOST_DEPLETION } from './modules/local/precheck_minimap2_host_depletion'
include { RUN_HOST_DEPLETION_MINIMAP2 } from './modules/local/run_host_depletion_minimap2'
include { PRECHECK_IRMA_SHORT } from './modules/local/precheck_irma_short'
include { RUN_IRMA_SHORT } from './modules/local/run_irma_short'
include { PRECHECK_IRMA_LONG } from './modules/local/precheck_irma_long'
include { RUN_IRMA_LONG } from './modules/local/run_irma_long'
include { RUN_EXTRACT_SEGMENTS } from './modules/local/run_extract_segments'
include { RUN_MERGE_SEGMENTS } from './modules/local/run_merge_segments'
include { RUN_ASSEMBLY_QC } from './modules/local/run_assembly_qc'
include { RUN_MERGE_ASSEMBLY_QC } from './modules/local/run_merge_assembly_qc'
include { PRECHECK_SAMTOOLS_DEPTH } from './modules/local/precheck_samtools_depth'
include { RUN_SAMTOOLS_DEPTH } from './modules/local/run_samtools_depth'
include { RUN_MERGE_DEPTH_SUMMARY } from './modules/local/run_merge_depth_summary'
include { PRECHECK_BLAST_TYPING } from './modules/local/precheck_blast_typing'
include { PREPARE_INFLUENZA_BLAST_DB } from './modules/local/prepare_influenza_blast_db'
include { RUN_BLAST_TYPING } from './modules/local/run_blast_typing'
include { PRECHECK_NEXTCLADE } from './modules/local/precheck_nextclade'
include { RUN_NEXTCLADE } from './modules/local/run_nextclade'
include { PRECHECK_MULTIQC } from './modules/local/precheck_multiqc'
include { RUN_MULTIQC } from './modules/local/run_multiqc'
include { PRECHECK_IVAR_CANONICAL } from './modules/local/precheck_ivar_canonical'
include { PREPARE_CANONICAL_REFS } from './modules/local/prepare_canonical_refs'
include { PRECHECK_FULLVARCALL } from './modules/local/precheck_fullvarcall'
include { PREPARE_REFSEQ_SEGMENTS } from './modules/local/prepare_refseq_segments'
include { RUN_FULLVARCALL_SAMPLE } from './modules/local/run_fullvarcall_sample'
include { RUN_MERGE_FULLVARCALL } from './modules/local/run_merge_fullvarcall'
include { RUN_IVAR_CANONICAL } from './modules/local/run_ivar_canonical'
include { PREPARE_ANTIVIRAL_DB } from './modules/local/prepare_antiviral_db'
include { RUN_ANTIVIRAL_RESISTANCE } from './modules/local/run_antiviral_resistance'
include { RUN_MEDAKA_CANONICAL } from './modules/local/run_medaka_canonical'
include { RUN_ANTIVIRAL_RESISTANCE_LONG } from './modules/local/run_antiviral_resistance_long'
include { RUN_H5_VIRULENCE } from './modules/local/run_h5_virulence'
include { PRECHECK_MEDAKA_VARIANTS } from './modules/local/precheck_medaka_variants'
include { RUN_MEDAKA_VARIANTS } from './modules/local/run_medaka_variants'
include { RUN_COINFECTION } from './modules/local/run_coinfection'
include { RUN_SURVEILLANCE_OUTPUTS } from './modules/local/run_surveillance_outputs'
include { RUN_LEGACY_PIPELINE } from './modules/local/run_legacy_pipeline'

def validateParams() {
    def errors = []
    def allowedSeqTypes = ['auto', 'short_paired', 'short_single', 'long']

    if( !params.input_dir ) {
        errors << "Missing required parameter: --input_dir"
    }

    if( !params.irma_module ) {
        errors << "Missing required parameter: --irma_module"
    }

    if( !allowedSeqTypes.contains(params.seq_type as String) ) {
        errors << "Invalid --seq_type '${params.seq_type}'. Allowed: ${allowedSeqTypes.join(', ')}"
    }

    if( params.run_legacy_bridge && !params.legacy_script ) {
        errors << "When --run_legacy_bridge is true, --legacy_script must be provided"
    }

    if( errors ) {
        errors.each { log.error it }
        System.exit(1)
    }
}

def yesNo(boolean value) {
    return value ? 'yes' : 'no'
}

def boolString(boolean value) {
    return value ? 'true' : 'false'
}

def sampleTuple(row) {
    def meta = [
        id      : row.sample_id as String,
        layout  : row.layout as String,
        seq_type: row.seq_type as String,
    ]

    def reads = [ file(row.r1 as String) ]
    if( row.r2 ) {
        reads << file(row.r2 as String)
    }

    return tuple(meta, reads)
}

def normalizeReadSet(reads) {
    def items = reads instanceof List ? reads.toList() : [reads]
    return items.sort { a, b -> a.getFileName().toString() <=> b.getFileName().toString() }
}

def hasUsableReads(meta, reads) {
    def items = normalizeReadSet(reads)
    if( meta.layout == 'paired' ) {
        return items.size() == 2 && items.every { it.exists() && it.size() > 0 }
    }
    return items.size() >= 1 && items[0].exists() && items[0].size() > 0
}

def chooseBestReads(meta, rawReads, trimmedReads, depletedReads = null) {
    def rawSet = normalizeReadSet(rawReads)
    def trimmedSet = trimmedReads != null ? normalizeReadSet(trimmedReads) : null
    def depletedSet = depletedReads != null ? normalizeReadSet(depletedReads) : null

    if( depletedSet != null && hasUsableReads(meta, depletedSet) ) {
        return depletedSet
    }
    if( trimmedSet != null && hasUsableReads(meta, trimmedSet) ) {
        return trimmedSet
    }
    return rawSet
}

workflow {
    validateParams()

    log.info "Starting MK Flu-Pipe Nextflow MVP"
    log.info "Input directory : ${params.input_dir}"
    log.info "Output directory: ${params.output_dir}"
    log.info "IRMA module     : ${params.irma_module}"
    log.info "Seq type        : ${params.seq_type}"
    log.info "Legacy bridge   : ${params.run_legacy_bridge}"

    DISCOVER_SAMPLES(
        params.input_dir,
        params.seq_type
    )

    PLAN_RUN(
        DISCOVER_SAMPLES.out.samplesheet,
        DISCOVER_SAMPLES.out.summary,
        params.seq_type as String,
        params.irma_module as String,
        boolString(params.run_fastqc as boolean),
        boolString(params.host_depletion as boolean),
        boolString(params.run_ivar as boolean),
        boolString(params.run_medaka as boolean),
        boolString(params.run_antiviral as boolean),
        boolString(params.run_h5_virulence as boolean),
        boolString(params.run_fullvarcall as boolean),
        boolString(params.run_legacy_bridge as boolean)
    )

    def routedSamples = DISCOVER_SAMPLES.out.samplesheet
        .splitCsv(header: true)
        .map { row -> sampleTuple(row) }
        .branch { meta, reads ->
            long_reads: meta.seq_type == 'long'
            short_reads: meta.seq_type != 'long'
        }

    def irmaConsensus = Channel.empty()
    def irmaDirs = Channel.empty()
    def shortReadsForVariants = Channel.empty()
    def longReadsForVariants = Channel.empty()
    def multiqcArtifacts = Channel.empty()
    def dashboardQcArtifacts = Channel.empty()

    if( params.run_fastqc ) {
        PRECHECK_FASTQC()

        RUN_FASTQC(
            PRECHECK_FASTQC.out.ok,
            routedSamples.short_reads.mix(routedSamples.long_reads)
        )

        multiqcArtifacts = multiqcArtifacts.mix(
            RUN_FASTQC.out.report_dirs.map { meta, reportDir -> reportDir }
        )
    }

    if( params.seq_type != 'long' ) {
        PRECHECK_FASTP()

        RUN_FASTP(
            PRECHECK_FASTP.out.ok,
            routedSamples.short_reads
        )

        multiqcArtifacts = multiqcArtifacts.mix(
            RUN_FASTP.out.reports.map { meta, report -> report }
        )
        dashboardQcArtifacts = dashboardQcArtifacts.mix(
            RUN_FASTP.out.reports.map { meta, report -> report }
        )
    }

    if( params.seq_type == 'long' ) {
        PRECHECK_FILTLONG()

        RUN_FILTLONG(
            PRECHECK_FILTLONG.out.ok,
            routedSamples.long_reads
        )

        multiqcArtifacts = multiqcArtifacts.mix(
            RUN_FILTLONG.out.reports.map { meta, report -> report }
        )
        dashboardQcArtifacts = dashboardQcArtifacts.mix(
            RUN_FILTLONG.out.reports.map { meta, report -> report }
        )
    }

    if( params.host_depletion ) {
        PREPARE_HUMAN_BOWTIE2_INDEX()
    }

    if( params.host_depletion && params.seq_type != 'long' ) {
        PRECHECK_BOWTIE2_HOST_DEPLETION()

        RUN_HOST_DEPLETION_BOWTIE2(
            PRECHECK_BOWTIE2_HOST_DEPLETION.out.ok,
            PREPARE_HUMAN_BOWTIE2_INDEX.out.human_fasta,
            PREPARE_HUMAN_BOWTIE2_INDEX.out.index_files,
            RUN_FASTP.out.cleaned_reads
        )

        dashboardQcArtifacts = dashboardQcArtifacts.mix(
            RUN_HOST_DEPLETION_BOWTIE2.out.reports.map { meta, report -> report }
        )
    }

    if( params.host_depletion && params.seq_type == 'long' ) {
        PRECHECK_MINIMAP2_HOST_DEPLETION()

        RUN_HOST_DEPLETION_MINIMAP2(
            PRECHECK_MINIMAP2_HOST_DEPLETION.out.ok,
            PREPARE_HUMAN_BOWTIE2_INDEX.out.human_fasta,
            RUN_FILTLONG.out.cleaned_reads
        )

        dashboardQcArtifacts = dashboardQcArtifacts.mix(
            RUN_HOST_DEPLETION_MINIMAP2.out.reports.map { meta, report -> report }
        )
    }


    if( params.seq_type != 'long' ) {
        PRECHECK_IRMA_SHORT()

        def shortRawById = routedSamples.short_reads
            .map { meta, reads -> tuple(meta.id, meta, reads) }
        def shortTrimmedById = RUN_FASTP.out.cleaned_reads
            .map { meta, reads -> tuple(meta.id, reads) }

        def shortReadsForIrma
        if( params.host_depletion ) {
            def shortDepletedById = RUN_HOST_DEPLETION_BOWTIE2.out.depleted_reads
                .map { meta, reads -> tuple(meta.id, reads) }

            shortReadsForIrma = shortRawById
                .join(shortTrimmedById)
                .join(shortDepletedById)
                .map { sampleId, meta, rawReads, trimmedReads, depletedReads ->
                    tuple(meta, chooseBestReads(meta, rawReads, trimmedReads, depletedReads))
                }
        }
        else {
            shortReadsForIrma = shortRawById
                .join(shortTrimmedById)
                .map { sampleId, meta, rawReads, trimmedReads ->
                    tuple(meta, chooseBestReads(meta, rawReads, trimmedReads, null))
                }
        }

        shortReadsForVariants = shortReadsForIrma

        RUN_IRMA_SHORT(
            PRECHECK_IRMA_SHORT.out.ok,
            shortReadsForIrma
        )

        irmaConsensus = irmaConsensus.mix(RUN_IRMA_SHORT.out.consensus_fastas)
        irmaDirs = irmaDirs.mix(RUN_IRMA_SHORT.out.irma_dirs)
    }

    if( params.seq_type == 'long' ) {
        PRECHECK_IRMA_LONG()

        def longRawById = routedSamples.long_reads
            .map { meta, reads -> tuple(meta.id, meta, reads) }
        def longFilteredById = RUN_FILTLONG.out.cleaned_reads
            .map { meta, reads -> tuple(meta.id, reads) }

        def longReadsForIrma
        if( params.host_depletion ) {
            def longDepletedById = RUN_HOST_DEPLETION_MINIMAP2.out.depleted_reads
                .map { meta, reads -> tuple(meta.id, reads) }

            longReadsForIrma = longRawById
                .join(longFilteredById)
                .join(longDepletedById)
                .map { sampleId, meta, rawReads, filteredReads, depletedReads ->
                    tuple(meta, chooseBestReads(meta, rawReads, filteredReads, depletedReads))
                }
        }
        else {
            longReadsForIrma = longRawById
                .join(longFilteredById)
                .map { sampleId, meta, rawReads, filteredReads ->
                    tuple(meta, chooseBestReads(meta, rawReads, filteredReads, null))
                }
        }

        longReadsForVariants = longReadsForIrma

        RUN_IRMA_LONG(
            PRECHECK_IRMA_LONG.out.ok,
            longReadsForIrma
        )

        irmaConsensus = irmaConsensus.mix(RUN_IRMA_LONG.out.consensus_fastas)
        irmaDirs = irmaDirs.mix(RUN_IRMA_LONG.out.irma_dirs)
    }

    RUN_EXTRACT_SEGMENTS(irmaConsensus)
    RUN_MERGE_SEGMENTS(RUN_EXTRACT_SEGMENTS.out.segment_files.flatten().collect())

    def assemblyQcInput = irmaDirs
        .map { meta, irmaDir -> tuple(meta.id, meta, irmaDir) }
        .join(
            irmaConsensus.map { meta, consensusFasta -> tuple(meta.id, meta, consensusFasta) }
        )
        .map { sampleId, metaLeft, irmaDir, metaRight, consensusFasta ->
            tuple(metaLeft, irmaDir, consensusFasta)
        }

    RUN_ASSEMBLY_QC(assemblyQcInput)
    RUN_MERGE_ASSEMBLY_QC(RUN_ASSEMBLY_QC.out.qc_rows.collect())

    PRECHECK_SAMTOOLS_DEPTH()
    RUN_SAMTOOLS_DEPTH(
        PRECHECK_SAMTOOLS_DEPTH.out.ok,
        irmaDirs
    )
    RUN_MERGE_DEPTH_SUMMARY(RUN_SAMTOOLS_DEPTH.out.depth_stats.collect())

    PRECHECK_MULTIQC()
    RUN_MULTIQC(
        PRECHECK_MULTIQC.out.ok,
        multiqcArtifacts.collect()
    )

    PRECHECK_BLAST_TYPING()
    PREPARE_INFLUENZA_BLAST_DB()
    RUN_BLAST_TYPING(
        PRECHECK_BLAST_TYPING.out.ok,
        PREPARE_INFLUENZA_BLAST_DB.out.blast_db_dir,
        RUN_EXTRACT_SEGMENTS.out.segment_files.flatten().collect()
    )

    PRECHECK_NEXTCLADE()
    RUN_NEXTCLADE(
        PRECHECK_NEXTCLADE.out.ok,
        RUN_BLAST_TYPING.out.summary,
        RUN_EXTRACT_SEGMENTS.out.segment_files
            .flatten()
            .filter { it.getFileName().toString().endsWith('_segment_4.fasta') }
            .collect()
    )

    def blastTypingRows = RUN_BLAST_TYPING.out.summary
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample as String, row.type_blast as String, row.subtype_HA as String, row.subtype_NA as String) }

    if( (params.seq_type != 'long' && params.run_ivar) || (params.seq_type == 'long' && params.run_antiviral && params.run_medaka) ) {
        PREPARE_CANONICAL_REFS()
    }

    if( params.run_fullvarcall ) {
        PRECHECK_FULLVARCALL()
        PREPARE_REFSEQ_SEGMENTS(
            RUN_BLAST_TYPING.out.summary,
            params.seq_type as String
        )
    }

    if( params.seq_type != 'long' && params.run_ivar ) {
        PRECHECK_IVAR_CANONICAL()

        def ivarInput = shortReadsForVariants
            .map { meta, reads -> tuple(meta.id, meta, reads) }
            .join(blastTypingRows)
            .map { sampleId, meta, reads, fluType, subtypeHa, subtypeNa ->
                tuple(meta, reads, fluType, subtypeHa, subtypeNa)
            }

        RUN_IVAR_CANONICAL(
            PRECHECK_IVAR_CANONICAL.out.ok,
            PREPARE_CANONICAL_REFS.out.refs_dir,
            ivarInput
        )

        if( params.run_antiviral ) {
            PREPARE_ANTIVIRAL_DB()
            RUN_ANTIVIRAL_RESISTANCE(
                PREPARE_ANTIVIRAL_DB.out.db,
                PREPARE_CANONICAL_REFS.out.refs_dir,
                RUN_BLAST_TYPING.out.summary,
                RUN_IVAR_CANONICAL.out.sample_dirs.collect()
            )
        }
    }

    if( params.seq_type == 'long' && params.run_medaka ) {
        PRECHECK_MEDAKA_VARIANTS()
        RUN_MEDAKA_VARIANTS(
            PRECHECK_MEDAKA_VARIANTS.out.ok,
            irmaDirs
        )

        if( params.run_antiviral ) {
            PREPARE_ANTIVIRAL_DB()

            def medakaCanonicalInput = longReadsForVariants
                .map { meta, reads -> tuple(meta.id, meta, reads) }
                .join(blastTypingRows)
                .map { sampleId, meta, reads, fluType, subtypeHa, subtypeNa ->
                    tuple(meta, reads, fluType, subtypeHa, subtypeNa)
                }

            RUN_MEDAKA_CANONICAL(
                PRECHECK_MEDAKA_VARIANTS.out.ok,
                PREPARE_CANONICAL_REFS.out.refs_dir,
                medakaCanonicalInput
            )

            RUN_ANTIVIRAL_RESISTANCE_LONG(
                PREPARE_ANTIVIRAL_DB.out.db,
                PREPARE_CANONICAL_REFS.out.refs_dir,
                RUN_BLAST_TYPING.out.summary,
                RUN_MEDAKA_CANONICAL.out.sample_dirs.collect()
            )
        }
    }

    if( params.run_fullvarcall ) {
        def fullVarcallReads = params.seq_type == 'long' ? longReadsForVariants : shortReadsForVariants
        def fullVarcallInput = fullVarcallReads
            .map { meta, reads -> tuple(meta.id, meta, reads) }
            .join(blastTypingRows)
            .map { sampleId, meta, reads, fluType, subtypeHa, subtypeNa ->
                tuple(meta, reads, fluType, subtypeHa, subtypeNa)
            }

        RUN_FULLVARCALL_SAMPLE(
            PRECHECK_FULLVARCALL.out.ok,
            PREPARE_REFSEQ_SEGMENTS.out.refs_dir,
            params.seq_type as String,
            fullVarcallInput
        )

        RUN_MERGE_FULLVARCALL(
            RUN_FULLVARCALL_SAMPLE.out.sample_dirs.collect()
        )
    }

    if( params.run_h5_virulence ) {
        RUN_H5_VIRULENCE(
            RUN_BLAST_TYPING.out.summary,
            irmaDirs
                .map { meta, irmaDir -> irmaDir }
                .collect()
        )
    }

    RUN_COINFECTION(
        irmaDirs
            .map { meta, irmaDir -> irmaDir }
            .collect()
    )

    def surveillanceDependencies = RUN_BLAST_TYPING.out.summary
        .mix(RUN_NEXTCLADE.out.summary)
        .mix(RUN_MERGE_ASSEMBLY_QC.out.report)
        .mix(RUN_MERGE_DEPTH_SUMMARY.out.summary)
        .mix(RUN_COINFECTION.out.report)
        .mix(RUN_MULTIQC.out.report)
        .mix(RUN_MULTIQC.out.data_dir)
        .mix(dashboardQcArtifacts)

    if( params.run_antiviral ) {
        surveillanceDependencies = surveillanceDependencies.mix(
            params.seq_type == 'long'
                ? RUN_ANTIVIRAL_RESISTANCE_LONG.out.report
                : RUN_ANTIVIRAL_RESISTANCE.out.report
        )
    }

    if( params.run_h5_virulence ) {
        surveillanceDependencies = surveillanceDependencies.mix(RUN_H5_VIRULENCE.out.report)
    }

    if( params.run_fullvarcall ) {
        surveillanceDependencies = surveillanceDependencies.mix(RUN_MERGE_FULLVARCALL.out.summary)
    }

    RUN_SURVEILLANCE_OUTPUTS(
        surveillanceDependencies.collect(),
        irmaConsensus
            .map { meta, consensusFasta -> consensusFasta }
            .collect(),
        irmaDirs
            .map { meta, irmaDir -> irmaDir }
            .collect(),
        params.irma_module as String,
        '0.1.0',
        (params.gisaid_location ?: '') as String,
        (params.gisaid_year ?: '') as String
    )

    if( params.run_legacy_bridge ) {
        RUN_LEGACY_PIPELINE(
            params.input_dir as String,
            params.output_dir as String,
            params.irma_module as String,
            (params.seq_mode ?: '') as String,
            params.seq_type as String,
            params.legacy_script as String,
            yesNo(params.run_fastqc as boolean),
            yesNo(params.host_depletion as boolean),
            yesNo(params.run_ivar as boolean),
            yesNo(params.run_medaka as boolean),
            yesNo(params.run_antiviral as boolean),
            yesNo(params.run_h5_virulence as boolean),
            yesNo(params.run_fullvarcall as boolean),
            params.adapter_fasta as String,
            params.min_len_short as String,
            params.min_len_long as String,
            params.max_len_long as String,
            params.min_qual as String,
            params.min_coverage as String,
            params.max_n_pct as String,
            params.min_segments as String,
            params.ivar_freq as String,
            params.ivar_depth as String,
            params.minority_freq as String,
            params.coinfection_pct as String,
            params.medaka_env as String,
            (params.gisaid_location ?: '') as String,
            (params.gisaid_year ?: '') as String
        )
    }
}
