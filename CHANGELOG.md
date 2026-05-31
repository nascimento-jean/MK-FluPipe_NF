# Changelog

All notable changes to MK-FluPipe NF are documented in this file.

The project follows semantic versioning whenever possible:

- `MAJOR` versions may introduce breaking workflow or output changes.
- `MINOR` versions add functionality while preserving existing behavior.
- `PATCH` versions fix bugs, improve documentation, or refine infrastructure.

## [Unreleased]

### Added

- Added tiny short-read and long-read `-stub-run` integration tests that traverse the main workflow and verify final dashboard, surveillance, metadata, phylogeny, GISAID, antiviral, H5, iVar, Medaka, and full-variant-call outputs.
- Added a structured functional annotation output for full variant calls (`functional_annotation/functional_annotation.tsv`) derived from RefSeq GFF3 and iVar protein mutation tables, with dashboard and download integration.
- Added optional experimental SnpEff annotation for short-read full variant calls (`--run_snpeff true`), with dashboard and download integration.

### Fixed

- Disabled the short fastp startup watchdog by default (`fastp_startup_timeout = 0`) to avoid killing healthy fastp jobs that take several minutes before producing final output files.
- Updated parameter validation so `--fastp_startup_timeout 0` is accepted as the documented "disabled" value.

## [v0.1.2] - 2026-05-30

### Added

- Added public GHCR container image support through the `ghcr` profile.
- Added Docker and Singularity/Apptainer usage paths for public container images.
- Added automated nf-test coverage for short-read, long-read, Medaka, antiviral, iVar, H5, full-variant-call, and surveillance-output execution paths using stubbed/minimal workflows.
- Added Python smoke tests for sample discovery, metadata validation, phylogeny helpers, light integration checks, and surveillance output generation.
- Added optional Augur-based phylogeny support for HA and NA, including user metadata integration and offline Auspice HTML links in the dashboard.
- Added validated sample metadata support through `--metadata_csv`, with reporting integration in the dashboard.

### Changed

- Updated the README and example parameter files to document public container usage, optional phylogeny inputs, metadata fields, and expanded runtime parameters.
- Improved GitHub Actions coverage for schema validation, profile validation, smoke tests, and nf-test.
- Updated package publication workflow for GHCR container images.
- Improved GISAID sequence output so `gisaid_sequences.fasta` follows the same segment-aware naming structure as `multisample_consensus.fasta`.

### Fixed

- Improved fastp execution robustness by writing uncompressed FASTQ files first and compressing them with external `pigz`, reducing the risk of fastp gzip-writer deadlocks in containerized runs.
- Added retry handling for timeout-killed fastp tasks.
- Improved handling of IRMA failures so failed samples can be recorded while successful samples continue through downstream analyses.
- Fixed dashboard subtype distribution so FluB samples are excluded from FluA subtype charts.

### Infrastructure

- Published container packages for the main MK-FluPipe tool image and Medaka image on GHCR.
- Added package/release documentation for running the pipeline without local container builds.

## [v0.1.1] - 2026-05-09

### Added

- Added Zenodo metadata for citation and DOI-oriented repository publication.
- Added repository license metadata.

### Changed

- Improved local Docker and Singularity execution settings.
- Tuned local execution resources and concurrency configuration for smoother workstation runs.

### Fixed

- Fixed Docker user/group write-permission handling in work directories.
- Improved Singularity bind/cache behavior for local execution.

## [v0.1.0] - 2026-05-03

### Added

- Initial public Nextflow DSL2 implementation of MK Flu-Pipe.
- Added short-read and long-read influenza workflow structure.
- Added core workflow stages for sample discovery, preprocessing, host depletion, IRMA assembly, segment extraction, assembly QC, typing/subtyping, Nextclade, antiviral resistance, H5 virulence screening, coinfection checks, full variant calls, and dashboard/report generation.
- Added initial documentation, workflow diagram, container build scripts, and GitHub release/package scaffolding.

[Unreleased]: https://github.com/nascimento-jean/MK-FluPipe_NF/compare/v0.1.2...HEAD
[v0.1.2]: https://github.com/nascimento-jean/MK-FluPipe_NF/releases/tag/v0.1.2
[v0.1.1]: https://github.com/nascimento-jean/MK-FluPipe_NF/releases/tag/v0.1.1
[v0.1.0]: https://github.com/nascimento-jean/MK-FluPipe_NF/releases/tag/v0.1.0
