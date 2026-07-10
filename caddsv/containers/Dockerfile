FROM condaforge/miniforge3:latest

ARG ENV_FILE
ARG CADD_SV_VERSION=development

LABEL org.opencontainers.image.title="CADD-SV rule environment"
LABEL org.opencontainers.image.description="Containerized dependencies for a CADD-SV Snakemake rule group"
LABEL org.opencontainers.image.source="https://github.com/kircherlab/CADD-SV"
LABEL org.opencontainers.image.version="${CADD_SV_VERSION}"

COPY ${ENV_FILE} /tmp/caddsv-environment.yml

RUN sed -i '/^name:/d; /^prefix:/d' /tmp/caddsv-environment.yml && \
    conda env create --name caddsv --file /tmp/caddsv-environment.yml --solver libmamba && \
    conda clean --all --yes && \
    rm /tmp/caddsv-environment.yml

# Apptainer exec does not rely on Docker entrypoint activation, so expose the
# isolated CADD-SV environment directly. Conda invokes libmamba only as its
# solver; no mamba command is used.
ENV PATH="/opt/conda/envs/caddsv/bin:/opt/conda/bin:${PATH}"
ENV CONDA_PREFIX="/opt/conda/envs/caddsv"
ENV PYTHONNOUSERSITE="1"
