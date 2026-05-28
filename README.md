<div align="center">
  <img src="docs/mk_flupipe_nextflow_workflow.svg" alt="MK Flu-Pipe Nextflow workflow" width="1100" />

# MK Flu-Pipe Nextflow

**A reproducible DSL2 Nextflow workflow for Influenza short-read and long-read genomic surveillance**

[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-23aa62?style=for-the-badge)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Containers-Docker-2496ED?style=for-the-badge)](https://www.docker.com/)
[![Singularity](https://img.shields.io/badge/Containers-Singularity%20%2F%20Apptainer-1f6feb?style=for-the-badge)](https://sylabs.io/docs/)
[![IRMA](https://img.shields.io/badge/IRMA-v1.3.2-8a2be2?style=for-the-badge)](https://hub.docker.com/r/cdcgov/irma)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20100567.svg)](https://doi.org/10.5281/zenodo.20100567)

</div>

---

## Contents
- [1. Overview](#1-overview)
- [2. What the pipeline does](#2-what-the-pipeline-does)
- [3. Current implementation status](#3-current-implementation-status)
- [4. Requirements](#4-requirements)
- [5. Installation and setup](#5-installation-and-setup)
- [6. Container strategy](#6-container-strategy)
- [7. Running the pipeline](#7-running-the-pipeline)
- [8. Parameters](#8-parameters)
- [9. Outputs](#9-outputs)
- [10. Databases and cache behavior](#10-databases-and-cache-behavior)
- [11. Frequently asked questions](#11-frequently-asked-questions)
- [12. Automated checks](#12-automated-checks)
- [13. Citation](#13-citation)

## 1. Overview

`MK Flu-Pipe Nextflow` is the Nextflow DSL2 implementation of the MK Flu-Pipe Influenza workflow. It supports Illumina short-read and Oxford Nanopore long-read data, from raw FASTQ files to final surveillance tables, consensus FASTA files, variant summaries, GISAID-ready files, MultiQC reports, and an interactive HTML dashboard.

The workflow supports:
- short-read Illumina data;
- long-read ONT data;
- automatic sample discovery;
- Docker and Singularity / Apptainer execution;
- Linux, Ubuntu, and WSL environments;
- tunable QC, preprocessing, assembly, variant calling, GISAID, and resource parameters.

## 2. What the pipeline does

### Short-read branch

1. Sample discovery and automatic run planning.
2. Raw read QC with `FastQC`.
3. Read preprocessing with `fastp`.
4. Optional host depletion with `Bowtie2`.
5. Influenza assembly with `IRMA` using `FLU` or `FLU-utr`.
6. Segment extraction and multi-sample segment FASTA generation.
7. Assembly QC and depth summaries.
8. Typing and subtyping with `BLAST`.
9. Clade assignment with `Nextclade`.
10. Canonical short-read variant calling with `iVar`.
11. Antiviral resistance screening.
12. H5 virulence marker screening when applicable.
13. Full protein mutation calling against RefSeq segment references and GFF3 annotations.
14. Coinfection / subtype mixing analysis.
15. Final dashboard, surveillance tables, MultiQC copy, FASTA exports, and optional GISAID-ready files.

### Long-read branch

1. Sample discovery and automatic run planning.
2. Raw read QC with `FastQC`.
3. Read preprocessing with `Filtlong`.
4. Optional host depletion with `minimap2`.
5. Influenza assembly with `IRMA` using `FLU-minion`.
6. Segment extraction and multi-sample segment FASTA generation.
7. Assembly QC and depth summaries.
8. Typing and subtyping with `BLAST`.
9. Clade assignment with `Nextclade`.
10. Canonical long-read variant calling with `Medaka`.
11. Antiviral resistance screening.
12. H5 virulence marker screening when applicable.
13. Full protein mutation calling against RefSeq segment references and GFF3 annotations.
14. Coinfection / subtype mixing analysis.
15. Final dashboard, surveillance tables, MultiQC copy, FASTA exports, and optional GISAID-ready files.

## 3. Current implementation status

The following modules are implemented in the current Nextflow version:
- sample discovery and run planning;
- `FastQC`;
- `fastp`;
- `Filtlong`;
- host depletion with `Bowtie2` and `minimap2`;
- `IRMA` short-read and long-read branches;
- IRMA failure status reporting without stopping the full run when a sample fails QC/assembly;
- segment extraction and segment merging;
- assembly QC and `samtools depth` summaries;
- `BLAST` typing;
- `Nextclade`;
- canonical short-read variant calling with `iVar`;
- canonical long-read variant calling with `Medaka`;
- antiviral resistance analysis;
- H5 virulence analysis;
- full protein mutation calling;
- coinfection and subtype mixing analysis;
- MultiQC aggregation;
- GISAID-ready FASTA and metadata exports;
- interactive HTML dashboard and final surveillance output generation.

## 4. Requirements

Recommended environment:
- Linux, Ubuntu, or WSL;
- `Nextflow` installed and available in `PATH`;
- either `Docker` or `Singularity / Apptainer`;
- internet access during the first run for database and image downloads;
- at least 8 GB RAM for very small tests;
- 16 GB to 32 GB RAM recommended for larger short-read runs;
- sufficient disk space for `work/`, downloaded databases, intermediate FASTQ files, and final results.

## 5. Installation and setup

### 5.1. Clone the repository

```bash
git clone https://github.com/nascimento-jean/MK-FluPipe_NF.git
cd MK-FluPipe_NF
```

### 5.2. Confirm Nextflow

```bash
nextflow -version
```

### 5.3. What is versioned in GitHub

The GitHub repository stores workflow source code, modules, helper scripts, documentation, and container build recipes.

The repository does not store:
- prebuilt Docker images;
- prebuilt Singularity `.sif` images;
- downloaded databases under `mk_flupipe_db/`;
- execution outputs;
- `work/` directories;
- Nextflow cache files.

After cloning, the normal flow is:
1. build or pull the required containers;
2. run the pipeline;
3. let the pipeline download or rebuild required databases automatically.

## 6. Container strategy

The workflow uses three container groups:
- `irma_tools`: resolved to `cdcgov/irma:v1.3.2`;
- `mk_flu_tools`: local image with the main workflow tool stack;
- `medaka_tools`: local image with Medaka-related tools.

### 6.1. Build local Docker images

```bash
bash containers/build_docker_images.sh
```

This creates:
- `mk-flu-pipe/mk_flu_tools:local`
- `mk-flu-pipe/medaka_tools:local`

### 6.2. Build local Singularity / Apptainer images

```bash
bash containers/build_singularity_images.sh
```

This creates:
- `containers/sif/mk_flu_tools_local.sif`
- `containers/sif/medaka_tools_local.sif`

## 7. Running the pipeline

### 7.1. Recommended profiles

Use the base `linux` profile together with one execution backend:

```bash
-profile linux,docker
```

or:

```bash
-profile linux,singularity
```

The `wsl` and `ubuntu` profiles remain as compatibility aliases, but `linux` is the recommended profile.

### 7.2. Short-read example

```bash
nextflow run main.nf \
  -resume \
  -profile linux,docker \
  --input_dir /path/to/FLU/ \
  --output_dir mk-flupipe_short_results \
  --irma_module FLU-utr \
  --seq_type short \
  --host_depletion true \
  --run_ivar true \
  --run_antiviral true \
  --run_h5_virulence true \
  --run_fullvarcall true
```

### 7.3. Long-read example

```bash
nextflow run main.nf \
  -resume \
  -profile linux,singularity \
  --input_dir /path/to/FLU_long/ \
  --output_dir mk-flupipe_long_results \
  --irma_module FLU-minion \
  --seq_type long \
  --host_depletion true \
  --run_medaka true \
  --run_antiviral true \
  --run_h5_virulence true \
  --run_fullvarcall true
```

### 7.4. Example with QC, GISAID, and resource parameters

```bash
nextflow run main.nf \
  -resume \
  -profile linux,docker \
  --input_dir /path/to/FLU/ \
  --output_dir mk-flupipe_results \
  --irma_module FLU-utr \
  --seq_type auto \
  --host_depletion true \
  --adapter_fasta /path/to/adapters.fa \
  --min_len_short 75 \
  --min_qual 20 \
  --min_coverage 50 \
  --max_n_pct 10 \
  --min_segments 4 \
  --ivar_freq 0.03 \
  --ivar_depth 10 \
  --minority_freq 0.20 \
  --coinfection_pct 5.0 \
  --gisaid_location Brazil-AL \
  --gisaid_year 2026 \
  --metadata_csv /path/to/sample_metadata.csv \
  --run_phylogeny true \
  --phylogeny_context_fasta /path/to/context_ha_na.fasta \
  --phylogeny_context_metadata /path/to/context_ha_na.csv \
  --max_cpus 8 \
  --max_memory "24 GB" \
  --queue_size 2
```

### 7.5. Default behavior

All parameters have defaults in `nextflow.config`. If you do not provide a parameter on the command line, the pipeline uses the default value.

The formal parameter contract is available in:

```text
nextflow_schema.json
```

Example parameter files are available in:

```text
params.example.yml
params.long.example.yml
```

## 8. Parameters

### 8.1. Required and core run parameters

| Parameter | Default | Description |
|---|---:|---|
| `--input_dir` | `null` | Folder containing input FASTQ / FASTQ.GZ files. Required for normal runs. |
| `--output_dir` | `${projectDir}/results` | Folder where all outputs will be written. |
| `--irma_module` | `null` | IRMA module to use. Common values are `FLU-utr` for short reads and `FLU-minion` for long reads. |
| `--seq_type` | `auto` | Sequencing type. Use `auto`, `short`, or `long`. |
| `--seq_mode` | empty | Optional discovery mode hint used by sample discovery and legacy compatibility. Usually left empty. |

For ONT long-read runs, use `--seq_type long` and `--irma_module FLU-minion`. The pipeline intentionally does not infer long-read mode from `--seq_type auto`.

### 8.2. Analysis switches

| Parameter | Default | Description |
|---|---:|---|
| `--run_fastqc` | `true` | Run raw-read FastQC. |
| `--host_depletion` | `false` | Enable host depletion. Short reads use `Bowtie2`; long reads use `minimap2`. |
| `--run_ivar` | `false` | Enable canonical short-read variant calling with `iVar`. |
| `--run_medaka` | `false` | Enable canonical long-read variant calling with `Medaka`. |
| `--run_antiviral` | `true` | Enable antiviral resistance analysis. |
| `--run_h5_virulence` | `true` | Enable H5 virulence marker analysis. |
| `--run_fullvarcall` | `false` | Enable full protein mutation calling against RefSeq segment references. |

For long-read antiviral resistance analysis, `--run_medaka true` is required because canonical long-read variants are produced with Medaka.
| `--run_legacy_bridge` | `false` | Run the optional legacy Bash bridge after the Nextflow workflow. Normally disabled. |
| `--legacy_script` | empty | Path to the optional legacy Bash script used only when `--run_legacy_bridge true`. |

### 8.3. Short-read preprocessing parameters

| Parameter | Default | Description |
|---|---:|---|
| `--adapter_fasta` | empty | Optional adapter FASTA passed to `fastp`. If empty, paired-end adapter detection is used. |
| `--min_len_short` | `75` | Minimum read length retained by `fastp` for short reads. |
| `--min_qual` | `20` | Minimum qualified base quality threshold used by `fastp`. |

### 8.4. Long-read preprocessing parameters

| Parameter | Default | Description |
|---|---:|---|
| `--min_len_long` | `200` | Minimum long-read length retained by `Filtlong`. |
| `--max_len_long` | `0` | Maximum long-read length retained by `Filtlong`; `0` disables the upper limit. |
| `--filtlong_min_mean_q` | `null` | Optional minimum mean read quality for `Filtlong`. If null, this filter is not applied. |

### 8.5. Assembly QC and segment filters

| Parameter | Default | Description |
|---|---:|---|
| `--min_coverage` | `50` | Minimum coverage threshold used by assembly QC and related summaries. |
| `--max_n_pct` | `10` | Maximum allowed percentage of `N` bases before a segment/sample is flagged. |
| `--min_segments` | `4` | Minimum number of detected segments expected for downstream reporting. |

### 8.6. Variant, resistance, and coinfection parameters

| Parameter | Default | Description |
|---|---:|---|
| `--ivar_freq` | `0.03` | Minimum allele frequency used by `iVar` variant calling. |
| `--ivar_depth` | `10` | Minimum depth used by `iVar` variant calling. |
| `--minority_freq` | `0.20` | Frequency threshold used for minority variant interpretation. |
| `--coinfection_pct` | `5.0` | Percentage threshold used to flag possible coinfection or subtype mixing. |
| `--medaka_env` | `medaka_env` | Medaka environment name retained for compatibility with legacy logic. Containerized execution does not require activating this environment manually. |

### 8.7. GISAID parameters

| Parameter | Default | Description |
|---|---:|---|
| `--gisaid_location` | empty | Location string used to create GISAID-style isolate names and to enable `GISAID_ready/` outputs. If empty, GISAID-ready files are not generated. |
| `--gisaid_year` | `null` | Year used in GISAID-style isolate names. If empty, the current year is used. |

### 8.8. Sample metadata for phylogenetic extensions

| Parameter | Default | Description |
|---|---:|---|
| `--metadata_csv` | empty | Optional CSV with sample metadata to validate, report, and use for phylogenetic analyses. Required columns are `sample_name` and `collection_date`; optional geographic columns include `country`, `state`, `city`, and `location`. |

When provided, `sample_name` must exactly match every sample ID discovered from the FASTQ filenames and `collection_date` must use ISO format (`YYYY-MM-DD`). Each discovered sample must occur exactly once. For standard Illumina paired files such as `261118000051_S40_L001_R1_001.fastq.gz`, the technical sample-number suffix is omitted and the expected metadata identifier is `261118000051`.

Example:

```csv
sample_name,collection_date,country,state,city
SAMPLE_A,2026-05-24,Brazil,Alagoas,Maceio
```

The validated table is exported as `Surveillance_Outputs/metadata.csv` and is displayed in a `Metadata` dashboard tab.

### 8.9. Optional HA/NA phylogeny with Augur

The optional phylogeny module generates trees only for the biologically informative surface segments `HA` and `NA`. Influenza A trees are separated by subtype (for example, `A_H3_HA` and `A_N2_NA`), while Influenza B produces `B_HA` and `B_NA` trees. The module uses the consensus FASTA generated by the pipeline, so the same implementation is available for validated short-read and long-read runs.

| Parameter | Default | Description |
|---|---:|---|
| `--run_phylogeny` | `false` | Run the optional Augur HA/NA phylogeny module. Requires `--metadata_csv`. |
| `--phylogeny_context_fasta` | empty | Optional contextual HA/NA FASTA selected from public NCBI records or manually downloaded by an authorized GISAID user. |
| `--phylogeny_context_metadata` | empty | CSV/TSV metadata matching context FASTA records. It must be supplied together with `--phylogeny_context_fasta`. |
| `--phylogeny_min_sequences` | `3` | Minimum total sequence count required to generate each independent tree. |
| `--phylogeny_threads` | `4` | Maximum threads requested by Augur/MAFFT/IQ-TREE. |

Context metadata records must include `strain`, `collection_date`, `type`, and `segment`. Influenza A context records must also provide the relevant `subtype_HA` or `subtype_NA`; optional columns are `country`, `state`, `city`, and `source`. The `strain` value must match the corresponding FASTA identifier.

```csv
strain,collection_date,type,segment,subtype_HA,subtype_NA,country,state,source
NCBI_H3_HA_001,2025-10-15,A,HA,H3,N2,Brazil,Alagoas,NCBI
GISAID_B_NA_001,2025-11-03,B,NA,-,-,Brazil,Alagoas,GISAID
```

GISAID sequences are not downloaded automatically by this pipeline. If GISAID context is used, it must be downloaded by an authorized user and analyzed locally in accordance with GISAID terms. NCBI context can be supplied with the same interface; automated NCBI retrieval will be added only after a reproducible context-selection policy is defined.

### 8.10. Resource and concurrency parameters

| Parameter | Default | Description |
|---|---:|---|
| `--max_cpus` | `null` | Global CPU cap per process. If omitted, each process uses its configured default. |
| `--max_memory` | `null` | Global memory cap per process, for example `"24 GB"`. If omitted, each process uses its configured default. |
| `--queue_size` | `8` | Local executor queue size. Controls how many tasks Nextflow may submit concurrently. |
| `--fastqc_threads` | `2` | CPU threads requested by `RUN_FASTQC`. |
| `--fastp_threads` | `2` | CPU threads requested by `RUN_FASTP`. |
| `--host_depletion_threads` | `2` | CPU threads requested by host depletion processes. |
| `--irma_threads` | `4` | CPU threads requested by IRMA processes. |
| `--fastp_max_forks` | `2` | Maximum number of concurrent `RUN_FASTP` tasks. |
| `--fastp_timeout` | `1800` | Hard timeout in seconds for a `fastp` task. Exit code `124` is retried by Nextflow. |
| `--fastp_startup_timeout` | `300` | Startup watchdog in seconds. If `fastp` produces no output during this window, the attempt is killed and retried. |
| `--host_depletion_max_forks` | `2` | Maximum number of concurrent host depletion tasks. |
| `--irma_max_forks` | `2` | Maximum number of concurrent IRMA tasks. |
| `--phylogeny_threads` | `4` | Maximum threads requested by the optional Augur HA/NA tree process. |

### 8.11. Container parameters

| Parameter | Default | Description |
|---|---:|---|
| `--mk_flu_docker_image` | `mk-flu-pipe/mk_flu_tools:local` | Docker image for the main workflow tools. |
| `--medaka_docker_image` | `mk-flu-pipe/medaka_tools:local` | Docker image for Medaka tools. |
| `--irma_docker_image` | `cdcgov/irma:v1.3.2` | Docker image for IRMA. |
| `--mk_flu_singularity_image` | `${projectDir}/containers/sif/mk_flu_tools_local.sif` | Singularity image for main workflow tools. |
| `--medaka_singularity_image` | `${projectDir}/containers/sif/medaka_tools_local.sif` | Singularity image for Medaka tools. |
| `--irma_singularity_image` | `docker://cdcgov/irma:v1.3.2` | Singularity/Apptainer IRMA image source. |
| `--singularity_cache_dir` | `${HOME}/.singularity/cache` | Singularity cache directory. |
| `--container_uid` | `1000` | Container user ID helper value. |
| `--container_gid` | `1000` | Container group ID helper value. |

### 8.10. Database and reference parameters

| Parameter | Default | Description |
|---|---:|---|
| `--human_db_dir` | `${HOME}/mk_flupipe_db/human_genome` | Human genome and host depletion index directory. Docker/Singularity profiles reset this under `${projectDir}/mk_flupipe_db`. |
| `--human_fasta_name` | `GRCh38_no_alt.fna` | Human reference FASTA filename. |
| `--human_index_prefix` | `GRCh38` | Bowtie2 human index prefix. |
| `--human_genome_url` | NCBI GRCh38 URL | URL used to download the human reference FASTA. |
| `--blast_db_dir` | `${HOME}/mk_flupipe_db` | Influenza database/cache directory. Docker/Singularity profiles reset this under `${projectDir}/mk_flupipe_db`. |
| `--blast_db_fasta` | `${params.blast_db_dir}/influenza.fna` | Influenza BLAST FASTA path. |
| `--blast_db_prefix` | `${params.blast_db_dir}/influenza_blast_db` | BLAST database prefix. |
| `--blast_db_timestamp` | `${params.blast_db_dir}/.last_update` | Timestamp file used to decide whether the BLAST database should be refreshed. |
| `--blast_db_url` | NCBI Influenza FTP URL | URL used to download the Influenza BLAST FASTA. |
| `--blast_db_max_days` | `30` | Maximum BLAST database age before refresh. |
| `--nextclade_datasets_dir` | `${params.blast_db_dir}/nextclade_datasets` | Directory for cached Nextclade datasets. |
| `--nextclade_max_days` | `30` | Maximum Nextclade dataset age before refresh. |
| `--canonical_refs_dir` | `${params.blast_db_dir}/canonical_refs` | Directory for canonical references. |
| `--refseq_segments_dir` | `${params.blast_db_dir}/refseq_segments` | Directory for RefSeq segment references and GFF3 files used by full variant calling. |

## 9. Outputs

The workflow writes all results under `--output_dir`.

### 9.1. Top-level output folders

| Path | What it contains |
|---|---|
| `bootstrap/` | Discovered sample tables, run planning files, and metadata used to start the workflow. |
| `qc_reports/` | QC reports and tabular QC summaries used by the dashboard and MultiQC. |
| `preprocessed_reads/` | Reads after `fastp` or `Filtlong`. |
| `depleted_reads/` | Reads after host depletion with `Bowtie2` or `minimap2`. |
| `irma_runs_short/` | Per-sample IRMA short-read run directories. |
| `irma_runs_long/` | Per-sample IRMA long-read run directories. |
| `assembly_final/` | Final consensus FASTA files, segment FASTA files, typing, Nextclade, resistance, H5, coinfection, and assembly QC outputs. |
| `depth_per_position/` | Per-sample depth tables generated from final alignments. |
| `variant_calls/` | Canonical variant calling outputs from `iVar` or Medaka variant workflows. |
| `variant_calls_canonical_long/` | Long-read canonical Medaka outputs. |
| `full_variant_calls/` | Full protein mutation reports and merged protein mutation table. |
| `Surveillance_Outputs/` | Main final delivery folder with dashboard, integrated tables, FASTA exports, GISAID-ready files, and copied final reports. |
| `legacy_bridge/` | Optional outputs only when `--run_legacy_bridge true` is used. |

### 9.2. QC outputs

| Path | What it contains |
|---|---|
| `qc_reports/fastqc_raw/` | Raw FastQC output folders for each sample. |
| `qc_reports/fastp/` | `fastp` HTML and JSON reports for short-read runs. |
| `qc_reports/filtlong/` | `Filtlong` statistics tables for long-read runs. |
| `qc_reports/host_depletion_bowtie2/` | Short-read host depletion logs and statistics. |
| `qc_reports/host_depletion_minimap2/` | Long-read host depletion logs and statistics. |
| `qc_reports/assembly_qc/` | Per-sample assembly QC tables. |
| `qc_reports/samtools_depth/` | Per-sample depth summary tables. |
| `qc_reports/multiqc/` | Full MultiQC report folder. |

### 9.3. Assembly and typing outputs

| Path | What it contains |
|---|---|
| `assembly_final/*.fasta` | Final normalized consensus FASTA files copied from successful IRMA outputs. Degenerate bases are converted to `N`. |
| `assembly_final/irma_status.tsv` | Per-sample IRMA status report, including samples that failed to produce amended consensus sequences. |
| `assembly_final/segments/` | Single-segment FASTA files and merged multi-sample segment FASTA files. |
| `assembly_final/assembly_qc_report.tsv` | Merged assembly QC summary across samples. |
| `assembly_final/depth_summary.tsv` | Merged depth summary across samples. |
| `assembly_final/blast_results/blast_typing_summary.tsv` | BLAST-based type, HA, NA, and hit metadata summary. |
| `assembly_final/nextclade_results/nextclade_summary.tsv` | Nextclade clade, dataset, and QC summary. |
| `assembly_final/coinfection/coinfection_report.tsv` | Coinfection and subtype mixing summary per sample. |

### 9.4. Variant, resistance, and mutation outputs

| Path | What it contains |
|---|---|
| `variant_calls/` | Canonical variant calling outputs from `iVar` or Medaka variant runs. |
| `variant_calls_canonical_long/` | Long-read canonical Medaka outputs used for downstream interpretation. |
| `assembly_final/antiviral_resistance/antiviral_resistance.tsv` | Antiviral resistance calls based on canonical references. |
| `assembly_final/h5_virulence/h5_virulence_markers.tsv` | H5 virulence marker results when H5 is detected. |
| `full_variant_calls/*.fullvarcall` | Per-sample full protein mutation reports. |
| `full_variant_calls/all_samples_protein_mutations.tsv` | Consolidated protein mutation table across all processed samples. |

### 9.5. Final surveillance outputs

`Surveillance_Outputs/` is the main folder for end users.

| Path | What it contains |
|---|---|
| `Surveillance_Outputs/surveillance_report.html` | Interactive final HTML dashboard with Overview, QC Dashboard, Typing, Resistance and H5, Coinfection, Protein Mutations, optional Metadata/Phylogeny, and Downloads tabs. |
| `Surveillance_Outputs/multiqc_report.html` | Copy of the full MultiQC report when available. |
| `Surveillance_Outputs/typing_results.tsv` | Integrated typing table combining BLAST, Nextclade, assembly QC, and hit metadata. |
| `Surveillance_Outputs/coverage_per_segment.tsv` | Per-segment coverage summary. |
| `Surveillance_Outputs/preprocessing_summary.tsv` | `fastp` or `Filtlong` preprocessing summary used by the dashboard. |
| `Surveillance_Outputs/host_depletion_summary.tsv` | Read count and retention summary before and after host depletion. |
| `Surveillance_Outputs/run_summary.tsv` | Compact integrated per-sample run summary. |
| `Surveillance_Outputs/run_summary.json` | JSON version of the integrated run summary. |
| `Surveillance_Outputs/multisample_consensus.fasta` | Multi-sample final consensus FASTA. Segment headers are preserved and degenerate bases are converted to `N`. |
| `Surveillance_Outputs/metadata.csv` | Validated sample metadata exported only when `--metadata_csv` is provided. |
| `Surveillance_Outputs/phylogeny/phylogeny_summary.tsv` | Tree-generation status per HA/NA analysis group when `--run_phylogeny true`. |
| `Surveillance_Outputs/phylogeny/<group>/<group>.json` | Auspice-compatible JSON dataset for each generated HA/NA tree. |
| `Surveillance_Outputs/phylogeny/<group>/tree.nwk` | Refined Newick tree for each generated HA/NA group. |
| `Surveillance_Outputs/coinfection/coinfection_report.tsv` | Final coinfection table copied for dashboard and download access. |
| `Surveillance_Outputs/antiviral_resistance/antiviral_resistance.tsv` | Final antiviral resistance table copied for dashboard and download access. |
| `Surveillance_Outputs/h5_virulence/h5_virulence_markers.tsv` | Final H5 virulence marker table copied for dashboard and download access. |
| `Surveillance_Outputs/full_variant_calls/all_samples_protein_mutations.tsv` | Final consolidated protein mutation table copied for dashboard and download access. |
| `Surveillance_Outputs/README_outputs.txt` | Plain-text explanation of the final output folder. |

### 9.6. GISAID-ready outputs

GISAID outputs are generated only when `--gisaid_location` is provided.

| Path | What it contains |
|---|---|
| `Surveillance_Outputs/GISAID_ready/gisaid_sequences.fasta` | FASTA file intended for GISAID preparation. It uses the same sequence headers and segment separation as `multisample_consensus.fasta`, preserving per-sample/per-segment identifiers. |
| `Surveillance_Outputs/GISAID_ready/gisaid_metadata.csv` | CSV metadata template with isolate identifiers, isolate names, type/subtype, clade, location, host, originating lab, submitting lab, collection date, and submitter fields. |

### 9.7. Optional outputs

| Path | When it appears |
|---|---|
| `legacy_bridge/` | Only when `--run_legacy_bridge true`. |
| `variant_calls_canonical_long/` | Long-read runs with `--run_medaka true`. |
| `qc_reports/filtlong/` | Long-read runs. |
| `qc_reports/fastp/` | Short-read runs. |
| `depleted_reads/bowtie2/` | Short-read runs with `--host_depletion true`. |
| `depleted_reads/minimap2/` | Long-read runs with `--host_depletion true`. |
| `Surveillance_Outputs/GISAID_ready/` | Runs with `--gisaid_location` set. |
| `Surveillance_Outputs/metadata.csv` | Runs with `--metadata_csv` set and valid metadata supplied. |
| `Surveillance_Outputs/phylogeny/` | Runs with `--run_phylogeny true` and valid metadata; contains one group per eligible HA/NA tree. |

## 10. Databases and cache behavior

The workflow automatically creates and updates `mk_flupipe_db/` as needed. This includes:
- human genome and host depletion index;
- Influenza BLAST database;
- Nextclade datasets;
- canonical references;
- RefSeq segment references and GFF3 files for full protein mutation calling;
- antiviral resistance marker databases.

If `mk_flupipe_db/` is deleted, it will be rebuilt on the next run.

## 11. Frequently asked questions

### Do I need to provide every parameter?

No. Parameters are optional unless they are required for your specific run. If a parameter is omitted, the value defined in `nextflow.config` is used.

### Does the pipeline require Conda?

No. The recommended execution strategy is based on Docker or Singularity / Apptainer containers.

### Does IRMA need to be installed on the host system?

No. IRMA runs through:

```text
cdcgov/irma:v1.3.2
```

### Can I run this on WSL?

Yes. The `linux` profile has been validated on WSL.

### Can I run this on native Ubuntu?

Yes. The same `linux` profile is intended for native Ubuntu.

### Why can the `work/` directory become large?

Nextflow stores staged inputs, intermediate files, task logs, and cached task results in `work/`. This is normal and allows `-resume` to reuse completed tasks. The folder can become large during FASTQ preprocessing, host depletion, and assembly.

### What happens if one IRMA sample fails?

The pipeline records the failed sample in the IRMA status output and continues downstream with samples that produced usable consensus FASTA files.

### What happens to degenerate bases?

Final sequences in `assembly_final/`, `Surveillance_Outputs/multisample_consensus.fasta`, and GISAID-ready FASTA exports are normalized so that degenerate bases are converted to `N`.

## 12. Automated checks

The repository includes a lightweight GitHub Actions workflow in `.github/workflows/tests.yml`.

These checks do not run IRMA, BLAST, Nextclade, Medaka, or other containerized analysis steps. They are designed as fast guardrails for:

- validating `nextflow_schema.json`;
- checking Python syntax for sample discovery;
- checking Python syntax and behavior for optional sample metadata validation;
- confirming short-read discovery maps `--seq_type short` to `short_paired`;
- confirming standard Illumina `_S<number>` tokens are not retained in reported sample identifiers;
- confirming long-read discovery produces `seq_type=long`;
- confirming invalid parameter combinations fail before heavy tasks are launched;
- running `nf-test` pipeline and process tests for parameter validation, `DISCOVER_SAMPLES`, and `VALIDATE_METADATA`.

Run the local smoke test from the repository root with:

```bash
python3 tests/check_discover_samples.py
python3 tests/check_metadata.py
```

If `nf-test` is installed, run the local nf-test suite with:

```bash
nf-test test tests/nf-test --ci
```

If the local Nextflow history file causes a duration-related runtime error, run nf-test with an isolated Nextflow home:

```bash
NXF_HOME=/tmp/mkflupipe-nxfhome-test nf-test test tests/nf-test --ci
```

The GitHub workflow installs Nextflow and nf-test, then verifies selected parameter-validation failures, including missing input directories and invalid long-read module/Medaka combinations.

## 13. Citation

If you use MK Flu-Pipe Nextflow in your research, please cite:

> Nascimento, J. (2026). *MK Flu-Pipe Nextflow: A reproducible DSL2 workflow for Influenza short-read and long-read genomic surveillance* (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.20100567

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20100567.svg)](https://doi.org/10.5281/zenodo.20100567)
