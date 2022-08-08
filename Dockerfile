# syntax=docker/dockerfile:1
FROM nfcore/base:latest
LABEL authors="Malte Rühlemann & Eike Wacker" \
      description="Docker image for Metagenome-Assembled Genome Scoring Tool (MAGScoT)"
#RUN cd /opt && git clone https://github.com/ikmb/MAGScoT.git
COPY . /opt
WORKDIR /opt

RUN conda env create -f /opt/environment.yml && conda clean -a
ENV PATH /opt/conda/envs/magscot_env/bin::$PATH

CMD ["--help"]
ENTRYPOINT [ "Rscript", "/opt/MAGScoT.R" ]
