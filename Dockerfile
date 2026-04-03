FROM mambaorg/micromamba:1.5.10

USER root

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --chown=$MAMBA_USER:$MAMBA_USER <<'EOF' /tmp/environment.yml
name: strint3
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.7.12
  - samtools=1.21
  - bedtools=2.31.1
  - minimap2=2.28
  - bioframe=0.8.0
  - numpy=1.21.5
  - pysam=0.16.0
  - seaborn=0.12.2
  - matplotlib-base=3.5.3
  - pip
  - pip:
      - pandas==1.3.5
      - editdistance==0.6.1
      - tqdm
      - networkx
      - fast_edit_distance
EOF

RUN micromamba create -y -f /tmp/environment.yml && \
    micromamba clean --all --yes

SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/envs/strint3/bin:$PATH

CMD ["bash"]
