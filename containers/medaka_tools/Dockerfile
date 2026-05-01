FROM mambaorg/micromamba:2.3.3

ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends bash ca-certificates procps tini \
    && rm -rf /var/lib/apt/lists/*

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /workspace
USER root