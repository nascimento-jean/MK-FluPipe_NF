nextflow.enable.dsl = 2

include { DISCOVER_SAMPLES } from './modules/local/discover_samples'
include { VALIDATE_METADATA } from './modules/local/validate_metadata'
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
include { RUN_MERGE_IRMA_STATUS } from './modules/local/run_merge_irma_status'
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
include { RUN_FUNCTIONAL_ANNOTATION } from './modules/local/run_functional_annotation'
include { RUN_SNPEFF_ANNOTATION } from './modules/local/run_snpeff_annotation'
include { RUN_IVAR_CANONICAL } from './modules/local/run_ivar_canonical'
include { PREPARE_ANTIVIRAL_DB } from './modules/local/prepare_antiviral_db'
include { RUN_ANTIVIRAL_RESISTANCE } from './modules/local/run_antiviral_resistance'
include { RUN_MEDAKA_CANONICAL } from './modules/local/run_medaka_canonical'
include { RUN_ANTIVIRAL_RESISTANCE_LONG } from './modules/local/run_antiviral_resistance_long'
include { RUN_H5_VIRULENCE } from './modules/local/run_h5_virulence'
include { PRECHECK_MEDAKA_VARIANTS } from './modules/local/precheck_medaka_variants'
include { RUN_MEDAKA_VARIANTS } from './modules/local/run_medaka_variants'
include { RUN_COINFECTION } from './modules/local/run_coinfection'
include { RUN_PHYLOGENY } from './modules/local/run_phylogeny'
include { RUN_SURVEILLANCE_OUTPUTS } from './modules/local/run_surveillance_outputs'
include { RUN_LEGACY_PIPELINE } from './modules/local/run_legacy_pipeline'

def isIntegerLike(value) {
    return value != null && value.toString() ==~ /-?\d+/
}

def isDecimalLike(value) {
    return value != null && value.toString() ==~ /-?(\d+(\.\d*)?|\.\d+)/
}

def paramBool(value) {
    if( value instanceof Boolean ) {
        return value
    }
    if( value == null ) {
        return false
    }
    def normalized = value.toString().trim().toLowerCase()
    return ['true', '1', 'yes', 'y', 'on'].contains(normalized)
}

def validateParams() {
    def errors = []
    def warnings = []
    def allowedSeqTypes = ['auto', 'short', 'short_paired', 'short_single', 'long']
    def allowedSeqModes = ['', 'illumina_paired', 'sra_paired', 'generic_paired', 'single']
    def seqType = params.seq_type as String
    def seqMode = (params.seq_mode ?: '') as String
    def irmaModule = (params.irma_module ?: '') as String
    def runIvar = paramBool(params.run_ivar)
    def runMedaka = paramBool(params.run_medaka)
    def runAntiviral = paramBool(params.run_antiviral)
    def runFullvarcall = paramBool(params.run_fullvarcall)
    def runSnpeff = paramBool(params.run_snpeff)
    def runPhylogeny = paramBool(params.run_phylogeny)
    def runLegacyBridge = paramBool(params.run_legacy_bridge)

    if( !params.input_dir ) {
        errors << "Missing required parameter: --input_dir"
    }
    else {
        def inputDir = new File(params.input_dir as String)
        if( !inputDir.exists() ) {
            errors << "Input directory does not exist: ${params.input_dir}"
        }
        else if( !inputDir.isDirectory() ) {
            errors << "Input path is not a directory: ${params.input_dir}"
        }
    }

    if( !params.irma_module ) {
        errors << "Missing required parameter: --irma_module"
    }

    if( !allowedSeqTypes.contains(seqType) ) {
        errors << "Invalid --seq_type '${params.seq_type}'. Allowed: ${allowedSeqTypes.join(', ')}"
    }

    if( !allowedSeqModes.contains(seqMode) ) {
        errors << "Invalid --seq_mode '${params.seq_mode}'. Allowed: ${allowedSeqModes.findAll { it }.join(', ')} or empty"
    }

    if( seqType == 'long' && runIvar ) {
        errors << "--run_ivar is only supported for short-read runs. Use --run_medaka true for long reads."
    }

    if( seqType != 'long' && runMedaka ) {
        errors << "--run_medaka is only supported for long-read runs. Use --run_ivar true for short reads."
    }

    if( seqType == 'long' && irmaModule != 'FLU-minion' ) {
        errors << "Long-read runs require --irma_module FLU-minion"
    }

    if( seqType != 'long' && irmaModule == 'FLU-minion' ) {
        errors << "--irma_module FLU-minion requires --seq_type long. The pipeline does not infer long reads from --seq_type auto."
    }

    if( seqType == 'long' && seqMode != '' ) {
        errors << "--seq_mode is only used for short-read sample discovery. Leave --seq_mode empty for long-read runs."
    }

    if( seqType == 'long' && runAntiviral && !runMedaka ) {
        errors << "Long-read antiviral resistance analysis requires --run_medaka true because canonical long-read variants are produced with Medaka."
    }

    if( runSnpeff && !runFullvarcall ) {
        errors << "--run_snpeff true requires --run_fullvarcall true because SnpEff annotation uses full variant call outputs"
    }

    if( seqType == 'long' && params.adapter_fasta ) {
        warnings << "--adapter_fasta is ignored for long-read runs because preprocessing uses Filtlong, not fastp."
    }

    if( seqMode == 'single' && seqType == 'short_paired' ) {
        errors << "--seq_mode single cannot be combined with --seq_type short_paired"
    }

    if( params.adapter_fasta ) {
        def adapterPath = new File(params.adapter_fasta as String)
        if( !adapterPath.exists() || !adapterPath.isFile() ) {
            errors << "Adapter FASTA does not exist or is not a file: ${params.adapter_fasta}"
        }
    }

    if( params.metadata_csv ) {
        def metadataPath = new File(params.metadata_csv as String)
        if( !metadataPath.exists() || !metadataPath.isFile() ) {
            errors << "Metadata CSV does not exist or is not a file: ${params.metadata_csv}"
        }
        else if( !metadataPath.name.toLowerCase().endsWith('.csv') ) {
            errors << "--metadata_csv must point to a .csv file"
        }
    }

    if( runPhylogeny && !params.metadata_csv ) {
        errors << "--run_phylogeny true requires --metadata_csv so tree tips have collection dates and locations"
    }

    def phylogenyContextFasta = (params.phylogeny_context_fasta ?: '').toString().trim()
    def phylogenyContextMetadata = (params.phylogeny_context_metadata ?: '').toString().trim()
    if( phylogenyContextFasta || phylogenyContextMetadata ) {
        if( !runPhylogeny ) {
            warnings << "Phylogeny context files are ignored unless --run_phylogeny true is supplied."
        }
        if( !phylogenyContextFasta || !phylogenyContextMetadata ) {
            errors << "--phylogeny_context_fasta and --phylogeny_context_metadata must be supplied together"
        }
        else {
            def contextFastaPath = new File(phylogenyContextFasta)
            def contextMetadataPath = new File(phylogenyContextMetadata)
            if( !contextFastaPath.exists() || !contextFastaPath.isFile() ) {
                errors << "Phylogeny context FASTA does not exist or is not a file: ${phylogenyContextFasta}"
            }
            if( !contextMetadataPath.exists() || !contextMetadataPath.isFile() ) {
                errors << "Phylogeny context metadata does not exist or is not a file: ${phylogenyContextMetadata}"
            }
        }
    }

    if( runLegacyBridge && !params.legacy_script ) {
        errors << "When --run_legacy_bridge is true, --legacy_script must be provided"
    }
    else if( runLegacyBridge ) {
        def legacyPath = new File(params.legacy_script as String)
        if( !legacyPath.exists() || !legacyPath.isFile() ) {
            errors << "Legacy script does not exist or is not a file: ${params.legacy_script}"
        }
    }

    if( params.gisaid_year ) {
        def yearText = params.gisaid_year.toString()
        if( !(yearText ==~ /\d{4}/) ) {
            errors << "--gisaid_year must be a four-digit year, for example 2026"
        }
    }

    def positiveIntParams = [
        'min_len_short',
        'min_len_long',
        'min_coverage',
        'min_segments',
        'ivar_depth',
        'queue_size',
        'fastqc_threads',
        'fastp_threads',
        'host_depletion_threads',
        'irma_threads',
        'fastp_max_forks',
        'fastp_timeout',
        'host_depletion_max_forks',
        'irma_max_forks',
        'phylogeny_min_sequences',
        'phylogeny_threads'
    ]
    positiveIntParams.each { name ->
        def value = params[name]
        if( !isIntegerLike(value) || value.toString().toInteger() < 1 ) {
            errors << "--${name} must be a positive integer"
        }
    }

    if( !isIntegerLike(params.max_len_long) || params.max_len_long.toString().toInteger() < 0 ) {
        errors << "--max_len_long must be 0 or a positive integer"
    }
    if( !isIntegerLike(params.fastp_startup_timeout) || params.fastp_startup_timeout.toString().toInteger() < 0 ) {
        errors << "--fastp_startup_timeout must be 0 or a positive integer"
    }
    else if( isIntegerLike(params.min_len_long) && params.max_len_long.toString().toInteger() > 0 && params.max_len_long.toString().toInteger() < params.min_len_long.toString().toInteger() ) {
        errors << "--max_len_long must be greater than or equal to --min_len_long, or 0 to disable the upper limit"
    }

    def fractionParams = ['ivar_freq', 'minority_freq']
    fractionParams.each { name ->
        if( !isDecimalLike(params[name]) ) {
            errors << "--${name} must be a number between 0 and 1"
            return
        }
        def value = new BigDecimal(params[name].toString())
        if( value < 0 || value > 1 ) {
            errors << "--${name} must be between 0 and 1"
        }
    }

    def percentageParams = ['max_n_pct', 'coinfection_pct']
    percentageParams.each { name ->
        if( !isDecimalLike(params[name]) ) {
            errors << "--${name} must be a number between 0 and 100"
            return
        }
        def value = new BigDecimal(params[name].toString())
        if( value < 0 || value > 100 ) {
            errors << "--${name} must be between 0 and 100"
        }
    }

    if( !isIntegerLike(params.min_qual) || params.min_qual.toString().toInteger() < 0 || params.min_qual.toString().toInteger() > 93 ) {
        errors << "--min_qual must be between 0 and 93"
    }

    if( params.max_cpus != null && params.max_cpus.toString().trim() != '' && (!isIntegerLike(params.max_cpus) || params.max_cpus.toString().toInteger() < 1) ) {
        errors << "--max_cpus must be a positive integer when provided"
    }

    if( params.filtlong_min_mean_q != null && params.filtlong_min_mean_q.toString().trim() != '' && !isDecimalLike(params.filtlong_min_mean_q) ) {
        errors << "--filtlong_min_mean_q must be numeric when provided"
    }
    else if( params.filtlong_min_mean_q != null && params.filtlong_min_mean_q.toString().trim() != '' && new BigDecimal(params.filtlong_min_mean_q.toString()) < 0 ) {
        errors << "--filtlong_min_mean_q must be greater than or equal to 0"
    }

    if( errors ) {
        log.error "Parameter validation failed. See nextflow_schema.json for the documented parameter contract."
        errors.each { log.error it }
        System.exit(1)
    }

    warnings.each { log.warn it }
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
    def seqType = params.seq_type as String
    def runFastqc = paramBool(params.run_fastqc)
    def hostDepletion = paramBool(params.host_depletion)
    def runIvar = paramBool(params.run_ivar)
    def runMedaka = paramBool(params.run_medaka)
    def runAntiviral = paramBool(params.run_antiviral)
    def runH5Virulence = paramBool(params.run_h5_virulence)
    def runFullvarcall = paramBool(params.run_fullvarcall)
    def runSnpeff = paramBool(params.run_snpeff)
    def runPhylogeny = paramBool(params.run_phylogeny)
    def runLegacyBridge = paramBool(params.run_legacy_bridge)
    def metadataCsv = (params.metadata_csv ?: '').toString().trim()
    def phylogenyContextFasta = (params.phylogeny_context_fasta ?: '').toString().trim()
    def phylogenyContextMetadata = (params.phylogeny_context_metadata ?: '').toString().trim()

    log.info "Starting MK Flu-Pipe Nextflow MVP"
    log.info "Input directory : ${params.input_dir}"
    log.info "Output directory: ${params.output_dir}"
    log.info "IRMA module     : ${params.irma_module}"
    log.info "Seq type        : ${params.seq_type}"
    log.info "Seq mode        : ${params.seq_mode ?: 'auto'}"
    log.info "SnpEff          : ${runSnpeff ? 'enabled' : 'disabled'}"
    log.info "Phylogeny       : ${runPhylogeny ? 'HA/NA with Augur' : 'disabled'}"
    log.info "Legacy bridge   : ${runLegacyBridge}"

    DISCOVER_SAMPLES(
        params.input_dir,
        params.seq_type,
        (params.seq_mode ?: '') as String
    )

    def analysisSamplesheet = DISCOVER_SAMPLES.out.samplesheet
    def validatedMetadata = Channel.empty()
    if( metadataCsv ) {
        VALIDATE_METADATA(
            DISCOVER_SAMPLES.out.samplesheet,
            file(metadataCsv)
        )
        analysisSamplesheet = VALIDATE_METADATA.out.samplesheet
        validatedMetadata = VALIDATE_METADATA.out.metadata
    }

    PLAN_RUN(
        analysisSamplesheet,
        DISCOVER_SAMPLES.out.summary,
        params.seq_type as String,
        params.irma_module as String,
        boolString(runFastqc),
        boolString(hostDepletion),
        boolString(runIvar),
        boolString(runMedaka),
        boolString(runAntiviral),
        boolString(runH5Virulence),
        boolString(runFullvarcall),
        boolString(runLegacyBridge)
    )

    def routedSamples = analysisSamplesheet
        .splitCsv(header: true)
        .map { row -> sampleTuple(row) }
        .branch { meta, reads ->
            long_reads: meta.seq_type == 'long'
            short_reads: meta.seq_type != 'long'
        }

    def irmaConsensus = Channel.empty()
    def irmaDirs = Channel.empty()
    def irmaStatuses = Channel.empty()
    def shortReadsForVariants = Channel.empty()
    def longReadsForVariants = Channel.empty()
    def multiqcArtifacts = Channel.empty()
    def dashboardQcArtifacts = Channel.empty()

    if( runFastqc ) {
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

    if( hostDepletion ) {
        PREPARE_HUMAN_BOWTIE2_INDEX()
    }

    if( hostDepletion && seqType != 'long' ) {
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

    if( hostDepletion && seqType == 'long' ) {
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


    if( seqType != 'long' ) {
        PRECHECK_IRMA_SHORT()

        def shortRawById = routedSamples.short_reads
            .map { meta, reads -> tuple(meta.id, meta, reads) }
        def shortTrimmedById = RUN_FASTP.out.cleaned_reads
            .map { meta, reads -> tuple(meta.id, reads) }

        def shortReadsForIrma
        if( hostDepletion ) {
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
        irmaStatuses = irmaStatuses.mix(RUN_IRMA_SHORT.out.statuses)
    }

    if( seqType == 'long' ) {
        PRECHECK_IRMA_LONG()

        def longRawById = routedSamples.long_reads
            .map { meta, reads -> tuple(meta.id, meta, reads) }
        def longFilteredById = RUN_FILTLONG.out.cleaned_reads
            .map { meta, reads -> tuple(meta.id, reads) }

        def longReadsForIrma
        if( hostDepletion ) {
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
        irmaStatuses = irmaStatuses.mix(RUN_IRMA_LONG.out.statuses)
    }

    RUN_MERGE_IRMA_STATUS(
        irmaStatuses
            .map { meta, statusFile -> statusFile }
            .collect()
    )

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

    def phylogenyResults = Channel.empty()
    if( runPhylogeny ) {
        def contextFasta = phylogenyContextFasta
            ? file(phylogenyContextFasta)
            : file("${projectDir}/assets/phylogeny/empty_context.fasta")
        def contextMetadata = phylogenyContextMetadata
            ? file(phylogenyContextMetadata)
            : file("${projectDir}/assets/phylogeny/empty_context_metadata.csv")

        RUN_PHYLOGENY(
            RUN_BLAST_TYPING.out.summary,
            RUN_EXTRACT_SEGMENTS.out.segment_files
                .flatten()
                .filter { path ->
                    def name = path.getFileName().toString()
                    name.endsWith('_segment_4.fasta') || name.endsWith('_segment_6.fasta')
                }
                .collect(),
            validatedMetadata,
            contextFasta,
            contextMetadata,
            params.phylogeny_min_sequences as int
        )
        phylogenyResults = RUN_PHYLOGENY.out.results
    }

    def blastTypingRows = RUN_BLAST_TYPING.out.summary
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample as String, row.type_blast as String, row.subtype_HA as String, row.subtype_NA as String) }

    if( (seqType != 'long' && runIvar) || (seqType == 'long' && runAntiviral && runMedaka) ) {
        PREPARE_CANONICAL_REFS()
    }

    if( runFullvarcall ) {
        PRECHECK_FULLVARCALL()
        PREPARE_REFSEQ_SEGMENTS(
            RUN_BLAST_TYPING.out.summary,
            params.seq_type as String
        )
    }

    if( seqType != 'long' && runIvar ) {
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

        if( runAntiviral ) {
            PREPARE_ANTIVIRAL_DB()
            RUN_ANTIVIRAL_RESISTANCE(
                PREPARE_ANTIVIRAL_DB.out.db,
                PREPARE_CANONICAL_REFS.out.refs_dir,
                RUN_BLAST_TYPING.out.summary,
                RUN_IVAR_CANONICAL.out.sample_dirs.collect()
            )
        }
    }

    if( seqType == 'long' && runMedaka ) {
        PRECHECK_MEDAKA_VARIANTS()
        RUN_MEDAKA_VARIANTS(
            PRECHECK_MEDAKA_VARIANTS.out.ok,
            irmaDirs
        )

        if( runAntiviral ) {
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

    if( runFullvarcall ) {
        def fullVarcallReads = seqType == 'long' ? longReadsForVariants : shortReadsForVariants
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

        RUN_FUNCTIONAL_ANNOTATION(
            RUN_MERGE_FULLVARCALL.out.summary,
            RUN_BLAST_TYPING.out.summary
        )

        if( runSnpeff ) {
            RUN_SNPEFF_ANNOTATION(
                PREPARE_REFSEQ_SEGMENTS.out.refs_dir,
                RUN_BLAST_TYPING.out.summary,
                RUN_FULLVARCALL_SAMPLE.out.sample_dirs.collect()
            )
        }
    }

    if( runH5Virulence ) {
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

    if( runAntiviral ) {
        surveillanceDependencies = surveillanceDependencies.mix(
            seqType == 'long'
                ? RUN_ANTIVIRAL_RESISTANCE_LONG.out.report
                : RUN_ANTIVIRAL_RESISTANCE.out.report
        )
    }

    if( runH5Virulence ) {
        surveillanceDependencies = surveillanceDependencies.mix(RUN_H5_VIRULENCE.out.report)
    }

    if( runFullvarcall ) {
        surveillanceDependencies = surveillanceDependencies.mix(RUN_MERGE_FULLVARCALL.out.summary)
        surveillanceDependencies = surveillanceDependencies.mix(RUN_FUNCTIONAL_ANNOTATION.out.report)
        if( runSnpeff ) {
            surveillanceDependencies = surveillanceDependencies.mix(RUN_SNPEFF_ANNOTATION.out.report)
        }
    }

    if( metadataCsv ) {
        surveillanceDependencies = surveillanceDependencies.mix(validatedMetadata)
    }

    if( runPhylogeny ) {
        surveillanceDependencies = surveillanceDependencies.mix(phylogenyResults)
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

    if( runLegacyBridge ) {
        RUN_LEGACY_PIPELINE(
            params.input_dir as String,
            params.output_dir as String,
            params.irma_module as String,
            (params.seq_mode ?: '') as String,
            params.seq_type as String,
            params.legacy_script as String,
            yesNo(runFastqc),
            yesNo(hostDepletion),
            yesNo(runIvar),
            yesNo(runMedaka),
            yesNo(runAntiviral),
            yesNo(runH5Virulence),
            yesNo(runFullvarcall),
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
