# Container Profiles

This project supports containerized execution on Linux via:
- `-profile linux,docker`
- `-profile linux,singularity`

The execution labels are split as follows:
- `irma_tools`: IRMA-only processes (`cdcgov/irma:v1.3.2`)
- `mk_flu_tools`: main pipeline tools (`mk-flu-pipe/mk_flu_tools:local`)
- `medaka_tools`: Medaka-specific tools (`mk-flu-pipe/medaka_tools:local`)

## Build Docker images

```bash
bash containers/build_docker_images.sh
```

## Build Singularity images from Docker images

```bash
bash containers/build_singularity_images.sh
```

## Typical commands

Docker:

```bash
nextflow run main.nf \
  -resume \
  -profile linux,docker \
  --input_dir /path/to/input \
  --output_dir mk_flupipe_out \
  --irma_module FLU-utr
```

Singularity:

```bash
nextflow run main.nf \
  -resume \
  -profile linux,singularity \
  --input_dir /path/to/input \
  --output_dir mk_flupipe_out \
  --irma_module FLU-utr
```