# syntax=docker/dockerfile:1
###############################################################################
# preAlignmentQC pipeline image
#
# Containerises the tools used by GeneralSteps/preAlignmentQC.sh so that the
# pipeline runs deterministically across dev, staging, and production. The
# bash script itself is mounted from the host or built from this repo at
# /opt/ngsmodule (see usage below).
#
# Build:
#   docker build -f pipelines/docker/preAlignmentQC.Dockerfile \
#                -t ngsmodule/preAlignmentQC:1.0.0 .
#
# Run (the way the Celery worker invokes it):
#   docker run --rm \
#     -v /data/ngsmodule_work:/data/work \
#     -v $(pwd):/opt/ngsmodule:ro \
#     -e NGS_TASK_ID=$task_id \
#     -e NGS_STRICT=1 \
#     ngsmodule/preAlignmentQC:1.0.0 \
#     /opt/ngsmodule/GeneralSteps/preAlignmentQC.sh ...
#
# Reproducibility notes:
#   - Tool versions are pinned via conda (mamba) to the same versions
#     specified in CheckENV.sh.
#   - The image runs as a non-root user (uid 1000) so writes to a mounted
#     work directory match the host user when --user is omitted on Linux.
#   - For production, build with `--platform linux/amd64` and push to a
#     private registry; tag with a content-addressable digest in the
#     PipelineTemplate database row.
###############################################################################

FROM mambaorg/micromamba:1.5.6-jammy

LABEL org.opencontainers.image.title="ngsmodule-preAlignmentQC"
LABEL org.opencontainers.image.source="https://github.com/zhanghao-njmu/NGSmodule"
LABEL org.opencontainers.image.description="Pre-alignment QC: fastp, FastQC, fastq-screen, SortMeRNA, BBTools"
LABEL org.opencontainers.image.licenses="MIT"

USER root

# Install the bioinformatics tools used by GeneralSteps/preAlignmentQC.sh.
# Versions track CheckENV.sh; bump deliberately when promoting an image.
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        fastqc=0.12.1 \
        fastp=0.23.4 \
        fastq-screen=0.15.3 \
        bowtie2=2.5.3 \
        sortmerna=4.3.6 \
        bbmap=39.06 \
        pigz=2.8 \
        samtools=1.19 \
        coreutils \
        gawk \
        gzip \
    && micromamba clean --all --yes

# A small wrapper so subprocess calls don't need to know about micromamba.
ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# Non-root execution by default (the worker passes --user $(id -u) on prod).
RUN useradd -m -u 1000 -s /bin/bash ngs
USER ngs
WORKDIR /home/ngs

# Health check verifies the most critical tool is callable.
HEALTHCHECK --interval=30s --timeout=5s --retries=3 \
    CMD fastp --version >/dev/null 2>&1 || exit 1

# Default to bash; the worker provides the explicit pipeline script path.
CMD ["bash"]
