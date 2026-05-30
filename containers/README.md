# Container Profiles

This project supports containerized execution on Linux via:
- `-profile linux,docker`
- `-profile linux,singularity`

The execution labels are split as follows:
- `irma_tools`: IRMA-only processes (`cdcgov/irma:v1.3.2`)
- `mk_flu_tools`: main pipeline tools (`mk-flu-pipe/mk_flu_tools:local`)
- `medaka_tools`: Medaka-specific tools (`mk-flu-pipe/medaka_tools:local`)

Public GHCR images can be used with the extra `ghcr` profile:

```bash
-profile linux,docker,ghcr
```

or, with Singularity/Apptainer pulling from GHCR:

```bash
-profile linux,singularity,ghcr
```

The `ghcr` profile maps the local image names to:

- `ghcr.io/nascimento-jean/mk-flupipe-nf-mk-flu-tools:<tag>`
- `ghcr.io/nascimento-jean/mk-flupipe-nf-medaka-tools:<tag>`

Use `--container_tag` to select a specific published version, for example `--container_tag v0.1.2`.

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

Docker with GHCR images:

```bash
nextflow run main.nf \
  -resume \
  -profile linux,docker,ghcr \
  --container_tag latest \
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

Singularity/Apptainer with GHCR images:

```bash
nextflow run main.nf \
  -resume \
  -profile linux,singularity,ghcr \
  --container_tag latest \
  --input_dir /path/to/input \
  --output_dir mk_flupipe_out \
  --irma_module FLU-utr
```
