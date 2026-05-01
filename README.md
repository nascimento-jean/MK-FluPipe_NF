<div align="center">
  <img src="docs/mk_flupipe_nextflow_workflow.svg" alt="MK Flu-Pipe Nextflow workflow" width="1100" />

# MK Flu-Pipe Nextflow

**A reproducible DSL2 Nextflow workflow for Influenza short-read and long-read analysis**  
**Um workflow reprodutível em Nextflow DSL2 para análise de Influenza com short reads e long reads**

[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-23aa62?style=for-the-badge)](https://www.nextflow.io/)  
[![Docker](https://img.shields.io/badge/Containers-Docker-2496ED?style=for-the-badge)](https://www.docker.com/)  
[![Singularity](https://img.shields.io/badge/Containers-Singularity%20%2F%20Apptainer-1f6feb?style=for-the-badge)](https://sylabs.io/docs/)  
[![IRMA](https://img.shields.io/badge/IRMA-v1.3.2-8a2be2?style=for-the-badge)](https://hub.docker.com/r/cdcgov/irma)
</div>

---

## Português (Brasil)

### Índice
- [1. O que é este projeto?](#1-o-que-é-este-projeto)
- [2. O que o pipeline faz?](#2-o-que-o-pipeline-faz)
- [3. Estado atual da implementação](#3-estado-atual-da-implementação)
- [4. Requisitos mínimos](#4-requisitos-mínimos)
- [5. Instalação rápida](#5-instalação-rápida)
- [6. Estratégia de containers](#6-estratégia-de-containers)
- [7. Como executar o pipeline](#7-como-executar-o-pipeline)
- [8. Principais parâmetros](#8-principais-parâmetros)
- [9. Estrutura dos resultados](#9-estrutura-dos-resultados)
- [10. Bancos de dados baixados automaticamente](#10-bancos-de-dados-baixados-automaticamente)
- [11. Dúvidas frequentes](#11-dúvidas-frequentes)
- [12. English version](#12-english-version)

### 1. O que é este projeto?
`MK Flu-Pipe Nextflow` é a migração do pipeline original em Bash/GTK para **Nextflow DSL2**, com foco em:

- reprodutibilidade;
- execução modular;
- suporte a **short reads** e **long reads**;
- uso de **containers** para reduzir dependências manuais;
- documentação clara para usuários iniciantes e avançados.

O pipeline foi organizado para rodar em **Linux**, incluindo:
- **Ubuntu nativo**;
- **WSL (Windows Subsystem for Linux)**.

### 2. O que o pipeline faz?
O fluxo analítico cobre a cadeia principal de análise de Influenza:

1. descoberta das amostras;
2. planejamento automático da execução;
3. controle de qualidade bruto (`FastQC`);
4. pré-processamento:
   - `fastp` para short reads;
   - `Filtlong` para long reads;
5. remoção de hospedeiro:
   - `Bowtie2` para short reads;
   - `minimap2` para long reads;
6. montagem com `IRMA`;
7. extração dos 8 segmentos e QC pós-montagem;
8. tipagem por `BLAST`;
9. classificação por `Nextclade`;
10. variant calling canônico:
   - `iVar` para short reads;
   - `Medaka` para long reads;
11. resistência antiviral;
12. marcadores de virulência H5;
13. variant calling completo contra `RefSeq NC_* + GFF3` (`Step 10b`);
14. análise de coinfecção / subtype mixing;
15. consolidação final de saídas para vigilância e preparação para GISAID.

### 3. Estado atual da implementação
No momento, os seguintes módulos já foram implementados e testados com sucesso:

- descoberta de amostras;
- `FastQC`;
- `fastp`;
- `Filtlong`;
- remoção de hospedeiro com `Bowtie2` e `minimap2`;
- `IRMA` short e long;
- extração de segmentos;
- `assembly_qc` e `samtools depth`;
- tipagem por `BLAST`;
- `Nextclade`;
- `iVar` canônico;
- `Medaka`;
- resistência antiviral;
- virulência H5;
- `Step 10b` (`RefSeq + GFF3`);
- coinfecção;
- relatórios finais para vigilância.

### 4. Requisitos mínimos
Para executar a versão atual do pipeline, recomenda-se:

- Linux / Ubuntu ou WSL;
- `Nextflow` instalado e funcional;
- **Docker** ou **Singularity / Apptainer** instalado;
- acesso à internet na primeira execução para download de imagens e bancos;
- pelo menos 8 GB de RAM para testes pequenos;
- mais memória e CPUs para lotes maiores.

### 5. Instalação rápida
#### 5.1. Clonar ou copiar o projeto
Coloque o projeto em um diretório de trabalho, por exemplo:

```bash
/home/usuario/MK_Flu-Pipe-nextflow
```

#### 5.2. Confirmar o Nextflow
```bash
nextflow -version
```

#### 5.3. Escolher backend de container
Você pode usar uma destas opções:

- **Docker**
- **Singularity / Apptainer**

### 6. Estratégia de containers
O pipeline foi configurado para usar **3 grupos de imagens**:

- `irma_tools`
  - usa automaticamente `cdcgov/irma:v1.3.2`
- `mk_flu_tools`
  - imagem local com o stack principal do pipeline
- `medaka_tools`
  - imagem local para rotinas que dependem do `Medaka`

#### 6.1. Construir imagens Docker locais
No diretório raiz do projeto:

```bash
bash containers/build_docker_images.sh
```

Isso gera localmente:

- `mk-flu-pipe/mk_flu_tools:local`
- `mk-flu-pipe/medaka_tools:local`

#### 6.2. Construir imagens Singularity / Apptainer
Depois que as imagens Docker locais estiverem prontas:

```bash
bash containers/build_singularity_images.sh
```

Isso gera, por padrão:

- `containers/sif/mk_flu_tools_local.sif`
- `containers/sif/medaka_tools_local.sif`

### 7. Como executar o pipeline
#### 7.1. Perfis recomendados
Use sempre o perfil base `linux` junto com o backend desejado:

- `-profile linux,docker`
- `-profile linux,singularity`

Os perfis `wsl` e `ubuntu` continuam disponíveis como aliases de compatibilidade, mas o perfil recomendado é `linux`.

#### 7.2. Exemplo: short reads com Docker
```bash
nextflow run main.nf \
  -resume \
  -profile linux,docker \
  --input_dir /caminho/para/FLU/ \
  --output_dir mk-flupipe_short_results \
  --irma_module FLU-utr \
  --host_depletion true \
  --run_ivar true \
  --run_antiviral true \
  --run_h5_virulence true \
  --run_fullvarcall true
```

#### 7.3. Exemplo: long reads com Singularity
```bash
nextflow run main.nf \
  -resume \
  -profile linux,singularity \
  --input_dir /caminho/para/FLU_long/ \
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

#### 7.4. Exemplo com controle de recursos
```bash
nextflow run main.nf \
  -resume \
  -profile linux,docker \
  --input_dir /caminho/para/FLU/ \
  --output_dir mk-flupipe_results \
  --irma_module FLU-utr \
  --max_cpus 8 \
  --max_memory "24 GB" \
  --queue_size 4
```

### 8. Principais parâmetros
| Parâmetro | Descrição |
|---|---|
| `--input_dir` | Pasta contendo os FASTQ/FASTQ.GZ de entrada |
| `--output_dir` | Pasta onde os resultados serão gravados |
| `--irma_module` | Módulo do IRMA, por exemplo `FLU-utr` ou `FLU-minion` |
| `--seq_type` | `auto`, `long`, `short`, etc. |
| `--host_depletion` | Ativa ou desativa remoção de hospedeiro |
| `--run_ivar` | Ativa variant calling canônico por `iVar` (short) |
| `--run_medaka` | Ativa variant calling por `Medaka` (long) |
| `--run_antiviral` | Ativa análise de resistência antiviral |
| `--run_h5_virulence` | Ativa análise de virulência H5 |
| `--run_fullvarcall` | Ativa o `Step 10b` com `RefSeq + GFF3` |
| `--max_cpus` | Limite global de CPUs por processo |
| `--max_memory` | Limite global de memória por processo |
| `--queue_size` | Número máximo de tarefas locais em paralelo |

### 9. Estrutura dos resultados
Os diretórios principais normalmente incluem:

- `qc_reports/`
- `preprocessed_reads/`
- `depleted_reads/`
- `irma_runs_short/` ou `irma_runs_long/`
- `assembly_final/`
- `variant_calls/`
- `variant_calls_canonical_long/`
- `full_variant_calls/`
- `depth_per_position/`
- `Surveillance_Outputs/`

Arquivos de interesse frequente:

- `assembly_final/blast_results/blast_typing_summary.tsv`
- `assembly_final/nextclade_results/nextclade_summary.tsv`
- `assembly_final/antiviral_resistance/antiviral_resistance.tsv`
- `assembly_final/h5_virulence/h5_virulence_markers.tsv`
- `assembly_final/coinfection/coinfection_report.tsv`
- `full_variant_calls/all_samples_protein_mutations.tsv`
- `Surveillance_Outputs/surveillance_report.html`
- `Surveillance_Outputs/multisample_consensus.fasta`

### 10. Bancos de dados baixados automaticamente
O pipeline recria e popula `mk_flupipe_db/` automaticamente conforme necessário. Entre os recursos baixados/preparados estão:

- genoma humano / índice para host depletion;
- banco de `BLAST` de influenza;
- datasets do `Nextclade`;
- referências canônicas;
- referências `RefSeq NC_*` e arquivos `GFF3`;
- base de marcadores de resistência antiviral.

Se `mk_flupipe_db/` for removido, ele é reconstruído na próxima execução.

### 11. Dúvidas frequentes
#### O pipeline precisa de Conda?
Não. A estratégia recomendada agora é baseada em **Docker** e **Singularity**.

#### O IRMA precisa estar instalado no sistema?
Não necessariamente. O pipeline usa a imagem:

```text
cdcgov/irma:v1.3.2
```

#### Posso rodar no WSL?
Sim. O perfil `linux` foi validado em WSL.

#### Posso rodar no Ubuntu nativo?
Sim. O mesmo perfil `linux` foi pensado para Ubuntu nativo.

#### O que acontece com bases degeneradas?
As sequências finais usadas em:
- `assembly_final/`
- `Surveillance_Outputs/multisample_consensus.fasta`
- `Surveillance_Outputs/GISAID_ready/gisaid_sequences.fasta`

são normalizadas para converter bases degeneradas em `N`, seguindo a lógica do pipeline original.

---

## English Version

### Contents
- [1. What is this project?](#1-what-is-this-project)
- [2. What does the pipeline do?](#2-what-does-the-pipeline-do)
- [3. Current implementation status](#3-current-implementation-status)
- [4. Minimum requirements](#4-minimum-requirements)
- [5. Quick installation](#5-quick-installation)
- [6. Container strategy](#6-container-strategy)
- [7. How to run the pipeline](#7-how-to-run-the-pipeline)
- [8. Main parameters](#8-main-parameters)
- [9. Output structure](#9-output-structure)
- [10. Databases downloaded automatically](#10-databases-downloaded-automatically)
- [11. Frequently asked questions](#11-frequently-asked-questions)

### 1. What is this project?
`MK Flu-Pipe Nextflow` is the migration of the original Bash/GTK workflow to **Nextflow DSL2**, with emphasis on:

- reproducibility;
- modular execution;
- support for **short-read** and **long-read** Influenza data;
- container-based deployment to reduce manual dependencies;
- beginner-friendly and advanced-user documentation.

The workflow is intended to run on **Linux**, including:
- **native Ubuntu**;
- **WSL (Windows Subsystem for Linux)**.

### 2. What does the pipeline do?
The analytical workflow covers the main Influenza analysis chain:

1. sample discovery;
2. automatic execution planning;
3. raw-read QC (`FastQC`);
4. preprocessing:
   - `fastp` for short reads;
   - `Filtlong` for long reads;
5. host depletion:
   - `Bowtie2` for short reads;
   - `minimap2` for long reads;
6. assembly with `IRMA`;
7. segment extraction and post-assembly QC;
8. `BLAST` typing;
9. `Nextclade` classification;
10. canonical variant calling:
   - `iVar` for short reads;
   - `Medaka` for long reads;
11. antiviral resistance;
12. H5 virulence markers;
13. complete variant calling against `RefSeq NC_* + GFF3` (`Step 10b`);
14. coinfection / subtype mixing analysis;
15. final surveillance and GISAID-ready output consolidation.

### 3. Current implementation status
The following modules are already implemented and tested:

- sample discovery;
- `FastQC`;
- `fastp`;
- `Filtlong`;
- host depletion with `Bowtie2` and `minimap2`;
- `IRMA` short and long;
- segment extraction;
- `assembly_qc` and `samtools depth`;
- `BLAST` typing;
- `Nextclade`;
- canonical `iVar` workflow;
- `Medaka`;
- antiviral resistance;
- H5 virulence;
- `Step 10b` (`RefSeq + GFF3`);
- coinfection;
- surveillance outputs.

### 4. Minimum requirements
Recommended environment:

- Linux / Ubuntu or WSL;
- working `Nextflow` installation;
- **Docker** or **Singularity / Apptainer** installed;
- internet access on first run for image and database downloads;
- at least 8 GB RAM for small tests;
- more memory and CPUs for full runs.

### 5. Quick installation
#### 5.1. Clone or copy the project
Place the project in a working directory, for example:

```bash
/home/user/MK_Flu-Pipe-nextflow
```

#### 5.2. Confirm Nextflow
```bash
nextflow -version
```

#### 5.3. Choose a container backend
You may use one of the following:

- **Docker**
- **Singularity / Apptainer**

### 6. Container strategy
The workflow is configured around **3 image groups**:

- `irma_tools`
  - automatically uses `cdcgov/irma:v1.3.2`
- `mk_flu_tools`
  - local image with the main workflow stack
- `medaka_tools`
  - local image for `Medaka`-dependent routines

#### 6.1. Build local Docker images
From the project root:

```bash
bash containers/build_docker_images.sh
```

This creates:

- `mk-flu-pipe/mk_flu_tools:local`
- `mk-flu-pipe/medaka_tools:local`

#### 6.2. Build Singularity / Apptainer images
After the Docker images are available:

```bash
bash containers/build_singularity_images.sh
```

This creates by default:

- `containers/sif/mk_flu_tools_local.sif`
- `containers/sif/medaka_tools_local.sif`

### 7. How to run the pipeline
#### 7.1. Recommended profiles
Use the base `linux` profile together with the desired backend:

- `-profile linux,docker`
- `-profile linux,singularity`

The `wsl` and `ubuntu` profiles remain available as compatibility aliases, but `linux` is the recommended one.

#### 7.2. Example: short reads with Docker
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

#### 7.3. Example: long reads with Singularity
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

#### 7.4. Example with resource control
```bash
nextflow run main.nf \
  -resume \
  -profile linux,docker \
  --input_dir /path/to/FLU/ \
  --output_dir mk-flupipe_results \
  --irma_module FLU-utr \
  --max_cpus 8 \
  --max_memory "24 GB" \
  --queue_size 4
```

### 8. Main parameters
| Parameter | Description |
|---|---|
| `--input_dir` | Folder containing the input FASTQ/FASTQ.GZ files |
| `--output_dir` | Folder where results will be written |
| `--irma_module` | IRMA module, e.g. `FLU-utr` or `FLU-minion` |
| `--seq_type` | `auto`, `long`, `short`, etc. |
| `--host_depletion` | Enables or disables host depletion |
| `--run_ivar` | Enables canonical `iVar` variant calling (short) |
| `--run_medaka` | Enables `Medaka` variant calling (long) |
| `--run_antiviral` | Enables antiviral resistance analysis |
| `--run_h5_virulence` | Enables H5 virulence analysis |
| `--run_fullvarcall` | Enables `Step 10b` using `RefSeq + GFF3` |
| `--max_cpus` | Global CPU cap per process |
| `--max_memory` | Global memory cap per process |
| `--queue_size` | Maximum number of local tasks in parallel |

### 9. Output structure
Main output folders usually include:

- `qc_reports/`
- `preprocessed_reads/`
- `depleted_reads/`
- `irma_runs_short/` or `irma_runs_long/`
- `assembly_final/`
- `variant_calls/`
- `variant_calls_canonical_long/`
- `full_variant_calls/`
- `depth_per_position/`
- `Surveillance_Outputs/`

Common files of interest:

- `assembly_final/blast_results/blast_typing_summary.tsv`
- `assembly_final/nextclade_results/nextclade_summary.tsv`
- `assembly_final/antiviral_resistance/antiviral_resistance.tsv`
- `assembly_final/h5_virulence/h5_virulence_markers.tsv`
- `assembly_final/coinfection/coinfection_report.tsv`
- `full_variant_calls/all_samples_protein_mutations.tsv`
- `Surveillance_Outputs/surveillance_report.html`
- `Surveillance_Outputs/multisample_consensus.fasta`

### 10. Databases downloaded automatically
The workflow recreates and populates `mk_flupipe_db/` automatically as needed. Resources include:

- human genome / index for host depletion;
- influenza `BLAST` database;
- `Nextclade` datasets;
- canonical references;
- `RefSeq NC_*` references and `GFF3` files;
- antiviral resistance marker database.

If `mk_flupipe_db/` is deleted, it is rebuilt on the next run.

### 11. Frequently asked questions
#### Does the pipeline require Conda?
No. The recommended strategy is now **Docker** and **Singularity**.

#### Does IRMA need to be installed on the host system?
Not necessarily. The workflow uses:

```text
cdcgov/irma:v1.3.2
```

#### Can I run this on WSL?
Yes. The `linux` profile has been validated on WSL.

#### Can I run this on native Ubuntu?
Yes. The same `linux` profile is intended for native Ubuntu.

#### What happens to degenerate bases?
Final sequences used in:
- `assembly_final/`
- `Surveillance_Outputs/multisample_consensus.fasta`
- `Surveillance_Outputs/GISAID_ready/gisaid_sequences.fasta`

are normalized so that degenerate bases are converted to `N`, matching the original workflow logic.