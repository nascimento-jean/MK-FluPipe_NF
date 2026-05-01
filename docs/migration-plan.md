# Migration Plan

## Goal

Transform the current MK Flu-Pipe Bash pipeline into a maintainable DSL2
Nextflow project while preserving the validated scientific logic.

## Recommended migration phases

1. Bootstrap

- create Nextflow structure
- validate parameters
- discover samples and infer layout
- write samplesheet and execution plan
- optionally bridge to the legacy Bash script

2. Pre-assembly QC

- FastQC
- fastp for short reads
- Filtlong for long reads
- host depletion

3. Assembly

- IRMA assembly as a dedicated process or subworkflow
- assembly checkpoints handled by Nextflow caching instead of manual files

4. Typing and clade assignment

- segment extraction
- BLAST classification
- Nextclade analysis

5. Variant analysis

- samtools depth
- iVar canonical variant calling
- Medaka ONT branch
- co-infection analysis
- antiviral resistance
- H5 virulence markers
- full-genome RefSeq + GFF3 protein mutation branch

6. Reporting

- consolidate outputs
- run summary TSV/JSON
- surveillance HTML
- GISAID-ready FASTA and metadata

## Suggested module map

- `modules/local/discover_samples.nf`
- `modules/local/fastqc.nf`
- `modules/local/trim_short_reads.nf`
- `modules/local/trim_long_reads.nf`
- `modules/local/host_depletion_short.nf`
- `modules/local/host_depletion_long.nf`
- `modules/local/irma_assembly.nf`
- `modules/local/extract_segments.nf`
- `modules/local/blast_typing.nf`
- `modules/local/nextclade.nf`
- `modules/local/depth_summary.nf`
- `modules/local/ivar_variants.nf`
- `modules/local/medaka_variants.nf`
- `modules/local/coinfection_report.nf`
- `modules/local/antiviral_report.nf`
- `modules/local/h5_virulence_report.nf`
- `modules/local/fullvarcall.nf`
- `modules/local/surveillance_outputs.nf`

## Design choices for this MVP

- keep the first delivery small enough to test quickly in WSL
- preserve the original Bash pipeline as a fallback bridge
- replace manual checkpoint files with Nextflow process caching over time
- move complex inline Bash/Python logic into small helper scripts when needed
- prepare Conda-backed process environments early to improve portability to Ubuntu

## Implemented now

- bootstrap discovery of samples from FASTQ input
- execution-plan rendering
- routing foundation for short-read and long-read branches
- FastQC raw-read stage as the first real migrated Nextflow process
- fastp short-read trimming stage with paired-end and single-end support
- Filtlong long-read QC stage with configurable min/max length and mean quality
- Conda labels wired to existing `mk_flu` and `medaka_env` directories for Linux execution
- Long-read QC now includes before/after FASTQ statistics for retention tracking
- Short-read host depletion now prepares GRCh38 automatically and removes host reads with Bowtie2
- Long-read host depletion now reuses the prepared human FASTA with minimap2 + samtools
- Short-read IRMA assembly now consumes depleted reads when host depletion is active, otherwise fastp outputs
- Long-read IRMA assembly now consumes depleted reads when host depletion is active, otherwise Filtlong outputs
