# arraykaryo - container image for HPC (Docker build -> Apptainer/Singularity run).
#
#   docker build -t arraykaryo .
#   docker run --rm -v "$PWD:/data" -w /data arraykaryo inspect cohort.bcf
#
# On HPC via Apptainer (CI publishes ghcr.io/<owner>/arraykaryo with the tags `main` and `sha-<commit>`; there is no `latest`):
#   apptainer build arraykaryo.sif docker://ghcr.io/<owner>/arraykaryo:main
#   apptainer run arraykaryo.sif run --bcf cohort.bcf --out out/
ARG PY=3.12-slim
FROM python:${PY} AS runtime

LABEL org.opencontainers.image.title="arraykaryo" \
      org.opencontainers.image.description="Large chromosomal events in a genotyping-array cohort from its VCF/BCF export: QC, pedigree from the genotypes, sex-chromosome complement, near-chromosomal gains, losses, copy-neutral LOH and ROH with cell fraction, one report" \
      org.opencontainers.image.licenses="MIT"

RUN apt-get update && apt-get install -y --no-install-recommends bcftools tabix && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/arraykaryo
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir . && python -c "import pysam, numpy, matplotlib, cyvcf2, arraykaryo"

RUN useradd --create-home --uid 1000 runner
USER runner
WORKDIR /data
ENV MPLBACKEND=Agg MPLCONFIGDIR=/tmp/mpl

ENTRYPOINT ["arraykaryo"]
CMD ["--help"]
