# MutaXplorer images: docker build --build-arg ENV_FILE=envs/<env>.yml -t <name>:1.0.0 .
FROM mambaorg/micromamba:1.5.10
ARG ENV_FILE=envs/mutaxplorer.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER ${ENV_FILE} /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean -a -y
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH
LABEL org.opencontainers.image.title="mutaxplorer" org.opencontainers.image.version="1.0.0"
