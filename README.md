<div align="center">
  <img src="docs/mk_flupipe_nextflow_workflow.svg" alt="MK Flu-Pipe Nextflow workflow" width="1100" />

# MK Flu-Pipe Nextflow

**A reproducible DSL2 Nextflow workflow for Influenza short-read and long-read analysis**

[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-23aa62?style=for-the-badge)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Containers-Docker-2496ED?style=for-the-badge)](https://www.docker.com/)
[![Singularity](https://img.shields.io/badge/Containers-Singularity%20%2F%20Apptainer-1f6feb?style=for-the-badge)](https://sylabs.io/docs/)
[![IRMA](https://img.shields.io/badge/IRMA-v1.3.2-8a2be2?style=for-the-badge)](https://hub.docker.com/r/cdcgov/irma)
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
- [8. Main parameters](#8-main-parameters)
- [9. Outputs](#9-outputs)
- [10. Databases and cache behavior](#10-databases-and-cache-behavior)
- [11. GitHub releases and packages](#11-github-releases-and-packages)
- [12. Frequently asked questions](#12-frequently-asked-questions)

## 1. Overview
`MK Flu-Pipe Nextflow` is the Nextflow DSL2 implementation of the MK Flu-Pipe Influenza workflow. It is designed to be reproducible, modular, and friendly to both first-time and experienced users.

The workflow supports:
- short-read Illumina data;
- long-read ONT data;
- Docker and Singularity / Apptainer execution;
- Linux, native Ubuntu, and WSL environments.

## 2. What the pipeline does
The pipeline covers the core Influenza analysis chain from raw reads to final surveillance outputs.

### Short-read branch
1. Sample discovery and automatic run planning.
2. Raw read QC with `FastQC`.
3. Preprocessing with `fastp`.
4. Host depletion with `Bowtie2`.
5. Assembly with `IRMA` using `FLU` / `FLU-utr` modules.
6. Segment extraction and post-assembly QC.
7. Typing with `BLAST`.
8. Clade assignment with `Nextclade`.
9. Canonical variant calling with `iVar`.
10. Antiviral resistance screening.
11. H5 virulence marker screening when applicable.
12. Full protein mutation calling against `RefSeq NC_* + GFF3` (`Step 10b`).
13. Coinfection / subtype mixing analysis.
14. Final dashboard and surveillance output generation.

### Long-read branch
1. Sample discovery and automatic run planning.
2. Raw read QC with `FastQC`.
3. Preprocessing with `Filtlong`.
4. Host depletion with `minimap2`.
5. Assembly with `IRMA` using `FLU-minion`.
6. Segment extraction and post-assembly QC.
7. Typing with `BLAST`.
8. Clade assignment with `Nextclade`.
9. Canonical variant calling with `Medaka`.
10. Antiviral resistance screening.
11. H5 virulence marker screening when applicable.
12. Full protein mutation calling against `RefSeq NC_* + GFF3` (`Step 10b`).
13. Coinfection / subtype mixing analysis.
14. Final dashboard and surveillance output generation.

## 3. Current implementation status
The following modules are implemented and validated in the current Nextflow version:
- sample discovery;
- execution planning;
- `FastQC`;
- `fastp`;
- `Filtlong`;
- host depletion with `Bowtie2` and `minimap2`;
- `IRMA` short and long branches;
- segment extraction and segment merging;
- assembly QC and `samtools depth` summaries;
- `BLAST` typing;
- `Nextclade`;
- canonical short-read calling with `iVar`;
- canonical long-read calling with `Medaka`;
- antiviral resistance analysis;
- H5 virulence analysis;
- full protein mutation calling (`Step 10b`);
- coinfection analysis;
- MultiQC aggregation;
- interactive HTML dashboard and surveillance output generation.

## 4. Requirements
Recommended environment:
- Linux, Ubuntu, or WSL;
- `Nextflow` installed and working;
- either `Docker` or `Singularity / Apptainer` installed;
- internet access on the first run for database and image downloads;
- at least 8 GB RAM for very small tests;
- more memory and CPUs for full runs.

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

### 5.3. What is versioned in GitHub and what is generated later
The GitHub repository stores the workflow source code, documentation, container build recipes, and helper scripts.

The repository does **not** store:
- prebuilt Docker images;
- prebuilt Singularity `.sif` files;
- downloaded databases under `mk_flupipe_db/`;
- execution outputs;
- `work/` directories or caches.

After `git clone`, the normal flow is:
1. build the local workflow images;
2. run the pipeline;
3. let the pipeline download or rebuild the required databases automatically.

## 6. Container strategy
The workflow currently uses three container groups:
- `irma_tools`
  - resolved automatically to `cdcgov/irma:v1.3.2`;
- `mk_flu_tools`
  - local image with the main workflow tool stack;
- `medaka_tools`
  - local image with Medaka-related tools.

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
- `-profile linux,docker`
- `-profile linux,singularity`

The `wsl` and `ubuntu` profiles remain as compatibility aliases, but `linux` is the recommended profile.

### 7.2. Example: short reads with Docker
```bash
nextflow run main.nf \
  -resume \
  -profile linux,docker \
  --input_dir /path/to/FLU/ \
  --output_dir mk-flupipe_short_results \
  --irma_module FLU-utr \
  --host_depletion true \
  --run_ivar true \
  --run_antiviral true \
  --run_h5_virulence true \
  --run_fullvarcall true
```

### 7.3. Example: long reads with Singularity
```bash
nextflow run main.nf \
  -resume \
  -profile linux,singularity \
  --input_dir /path/to/FLU_long/ \
  --output_dir mk-flupipe_long_results \
  --irma_module FLU-minion \
  --seq_type long \
  --host_depletion true \
  --min_len_long 200 \
  --max_len_long 0 \
  --filtlong_min_mean_q 10 \
  --run_medaka true \
  --run_antiviral true \
  --run_h5_virulence true \
  --run_fullvarcall true
```

### 7.4. Example with resource control
```bash
nextflow run main.nf \
  -resume \
  -profile linux,singularity \
  --input_dir /path/to/FLU/ \
  --output_dir mk-flupipe_results \
  --irma_module FLU-utr \
  --max_cpus 8 \
  --max_memory "24 GB" \
  --queue_size 4
```

## 8. Main parameters
| Parameter | Description |
|---|---|
| `--input_dir` | Folder containing the input FASTQ / FASTQ.GZ files |
| `--output_dir` | Folder where results will be written |
| `--irma_module` | IRMA module such as `FLU-utr` or `FLU-minion` |
| `--seq_type` | `auto`, `short`, or `long` |
| `--host_depletion` | Enables or disables host depletion |
| `--run_ivar` | Enables canonical `iVar` calling for short reads |
| `--run_medaka` | Enables canonical `Medaka` calling for long reads |
| `--run_antiviral` | Enables antiviral resistance analysis |
| `--run_h5_virulence` | Enables H5 virulence analysis |
| `--run_fullvarcall` | Enables `Step 10b` full protein mutation calling |
| `--max_cpus` | Global CPU cap per process |
| `--max_memory` | Global memory cap per process |
| `--queue_size` | Maximum number of concurrent local tasks |

A sample parameter file is available in:
```text
params.example.yml
```

## 9. Outputs
This section summarizes the outputs generated by the pipeline and what each one contains.

### 9.1. Top-level folders created under `--output_dir`
| Path | What it contains |
|---|---|
| `bootstrap/` | Run planning files such as discovered sample tables and run metadata. |
| `qc_reports/` | QC outputs and tabular QC summaries used by the dashboard and MultiQC. |
| `preprocessed_reads/` | Reads after `fastp` or `Filtlong`. |
| `depleted_reads/` | Reads after host depletion with `Bowtie2` or `minimap2`. |
| `irma_runs_short/` | Per-sample IRMA short-read run directories. Present only for short-read runs. |
| `irma_runs_long/` | Per-sample IRMA long-read run directories. Present only for long-read runs. |
| `assembly_final/` | Final consensus assemblies, merged segment FASTA files, typing, clade, resistance, H5, and coinfection outputs. |
| `depth_per_position/` | Per-sample depth tables generated from final alignments. |
| `variant_calls/` | Canonical variant calling outputs from `iVar` or Medaka variant workflows. |
| `variant_calls_canonical_long/` | Long-read canonical Medaka outputs. Present only for long-read runs with `--run_medaka true`. |
| `full_variant_calls/` | Full protein mutation outputs generated in `Step 10b`. |
| `Surveillance_Outputs/` | Final dashboard, final integrated TSV files, run summaries, MultiQC copy, and multi-sample FASTA outputs. |
| `legacy_bridge/` | Optional outputs only when `--run_legacy_bridge true` is used. |

### 9.2. `qc_reports/`
| Path | What it contains |
|---|---|
| `qc_reports/fastqc_raw/` | Raw `FastQC` output folders for each sample. |
| `qc_reports/fastp/` | `fastp` HTML and JSON reports. Present for short-read runs. |
| `qc_reports/filtlong/` | `Filtlong` statistics tables. Present for long-read runs. |
| `qc_reports/host_depletion_bowtie2/` | Host depletion statistics tables for short-read runs. |
| `qc_reports/host_depletion_minimap2/` | Host depletion statistics tables for long-read runs. |
| `qc_reports/assembly_qc/` | Per-sample assembly QC tables generated from extracted segments. |
| `qc_reports/samtools_depth/` | Per-sample depth summary tables generated from `samtools depth`. |
| `qc_reports/multiqc/` | Full `MultiQC` report folder. |

### 9.3. `assembly_final/`
| Path | What it contains |
|---|---|
| `assembly_final/*.fasta` | Final normalized consensus FASTA files copied from IRMA outputs. Degenerate bases are converted to `N`. |
| `assembly_final/segments/` | Single-segment FASTA files and merged multi-sample segment FASTA files. |
| `assembly_final/assembly_qc_report.tsv` | Merged assembly QC summary across samples. |
| `assembly_final/depth_summary.tsv` | Merged depth summary across samples. |
| `assembly_final/blast_results/blast_typing_summary.tsv` | BLAST-based HA / NA typing summary. |
| `assembly_final/nextclade_results/nextclade_summary.tsv` | Nextclade clade, dataset, and QC summary. |
| `assembly_final/antiviral_resistance/antiviral_resistance.tsv` | Antiviral resistance calls based on canonical references. |
| `assembly_final/h5_virulence/h5_virulence_markers.tsv` | H5 virulence marker results when H5 is detected. |
| `assembly_final/coinfection/coinfection_report.tsv` | Coinfection and subtype mixing summary per sample. |

### 9.4. Variant and mutation outputs
| Path | What it contains |
|---|---|
| `variant_calls/` | Canonical variant calling outputs from `iVar` or Medaka variant runs. File types depend on branch and caller. |
| `variant_calls_canonical_long/` | Long-read canonical Medaka outputs used for downstream interpretation. |
| `full_variant_calls/*.fullvarcall` | Per-sample full protein mutation reports from `Step 10b`. |
| `full_variant_calls/all_samples_protein_mutations.tsv` | Consolidated protein mutation table across all processed samples. |

### 9.5. Final dashboard and surveillance outputs
`Surveillance_Outputs/` is the main delivery folder for end users.

| Path | What it contains |
|---|---|
| `Surveillance_Outputs/surveillance_report.html` | Interactive final HTML dashboard with tabs for overview, QC, typing, resistance, coinfection, protein mutations, and downloads. |
| `Surveillance_Outputs/multiqc_report.html` | A copy of the full MultiQC report for direct opening from the final output folder. |
| `Surveillance_Outputs/typing_results.tsv` | Integrated typing table combining BLAST, Nextclade, assembly QC, and hit metadata. |
| `Surveillance_Outputs/preprocessing_summary.tsv` | `fastp` or `Filtlong` summary table used by the dashboard preprocessing section. |
| `Surveillance_Outputs/host_depletion_summary.tsv` | Read count and retention summary before and after host depletion. |
| `Surveillance_Outputs/run_summary.tsv` | Compact integrated run summary per sample. |
| `Surveillance_Outputs/run_summary.json` | JSON version of the integrated run summary. |
| `Surveillance_Outputs/multisample_consensus.fasta` | Multi-sample final consensus FASTA. Degenerate bases are converted to `N`. |
| `Surveillance_Outputs/coinfection/coinfection_report.tsv` | Local copy of the final coinfection summary used by the dashboard. |
| `Surveillance_Outputs/full_variant_calls/all_samples_protein_mutations.tsv` | Local copy of the consolidated protein mutation table used by the dashboard. |
| `Surveillance_Outputs/README_outputs.txt` | Plain-text explanation of the main final outputs. |

### 9.6. Optional outputs
| Path | What it contains |
|---|---|
| `legacy_bridge/` | Outputs created only when the optional legacy Bash bridge is enabled. |
| `variant_calls_canonical_long/` | Only generated for long-read Medaka runs. |
| `qc_reports/filtlong/` | Only generated for long-read runs. |
| `qc_reports/fastp/` | Only generated for short-read runs. |
| `depleted_reads/bowtie2/` and `qc_reports/host_depletion_bowtie2/` | Only generated for short-read runs when host depletion is enabled. |
| `depleted_reads/minimap2/` and `qc_reports/host_depletion_minimap2/` | Only generated for long-read runs when host depletion is enabled. |

## 10. Databases and cache behavior
The workflow automatically recreates and populates `mk_flupipe_db/` as needed. This includes:
- the human genome and host depletion index;
- the Influenza BLAST database;
- Nextclade datasets;
- canonical references;
- `RefSeq NC_*` references and `GFF3` files for `Step 10b`;
- antiviral resistance marker databases.

If `mk_flupipe_db/` is deleted, it will be rebuilt on the next run.

## 11. GitHub releases and packages
The repository already hosts the workflow source code. The remaining GitHub-facing pieces are:
- a versioned **Release**;
- published **container packages**.

### 11.1. Release workflow
A typical first release flow is:
```bash
git add README.md docs/mk_flupipe_nextflow_workflow.svg .github/workflows/publish-ghcr.yml
git commit -m "Docs: refresh README and workflow diagram"
git push origin main

git tag -a v1.0.0 -m "MK Flu-Pipe Nextflow v1.0.0"
git push origin v1.0.0
```

Then open GitHub and create the Release from tag `v1.0.0`.

Recommended release assets or notes:
- workflow version and highlights;
- validated profiles (`linux,docker`, `linux,singularity`);
- validated short-read and long-read support;
- updated dashboard and MultiQC integration;
- any limitations still under active refinement.

### 11.2. Packages via GitHub Container Registry (GHCR)
This repository now includes a GitHub Actions workflow in:
```text
.github/workflows/publish-ghcr.yml
```

It publishes two Docker images to GitHub Container Registry:
- `ghcr.io/<owner>/mk-flupipe-nf-mk-flu-tools`
- `ghcr.io/<owner>/mk-flupipe-nf-medaka-tools`

The workflow can be triggered in two ways:
- automatically when you push a tag such as `v1.0.0`;
- manually from the **Actions** tab using **Run workflow**.

Once the workflow succeeds, those images will appear in the repository's **Packages** section on GitHub.

## 12. Frequently asked questions
### Does the pipeline require Conda?
No. The recommended execution strategy is based on Docker and Singularity / Apptainer.

### Does IRMA need to be installed on the host system?
No. The workflow uses:
```text
cdcgov/irma:v1.3.2
```

### Can I run this on WSL?
Yes. The `linux` profile has been validated on WSL.

### Can I run this on native Ubuntu?
Yes. The same `linux` profile is intended for native Ubuntu.

### Do I have to launch the workflow from inside the repository directory?
No. You may run the workflow by pointing `nextflow run` to the project directory or directly to `main.nf`, as long as you provide valid input and output paths.

### What happens to degenerate bases?
Final sequences in:
- `assembly_final/`;
- `Surveillance_Outputs/multisample_consensus.fasta`;
- optional downstream FASTA exports;

are normalized so that degenerate bases are converted to `N`, matching the original workflow logic.