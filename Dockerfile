# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:1.5.10
LABEL authors="Malte Rühlemann & Eike Wacker" \
      description="Docker image for the MAG Scoring and Refinement Tool (MAGScoT)"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/MAGScoT

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "Rscript", "/opt/MAGScoT/MAGScoT.R"]
CMD ["--help"]
