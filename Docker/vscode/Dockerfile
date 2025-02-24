FROM linuxserver/code-server:4.97.2

# Install base utilities
RUN apt-get update && \
    apt-get install -y wget bc tree && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget -O ~/miniconda.sh --quiet \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install mamba -y -n base -c conda-forge

RUN mamba create -y \
-c conda-forge -c bioconda \
-n lr-tools \
"samtools>=1.10" \
minimap2 \
fastqc \
pbmm2 \
parallel \
trgt \
ncurses

# install nanoplot with pip (with conda gives errors)
RUN pip install NanoPlot

RUN mamba create -y -n flair -c conda-forge -c bioconda flair

RUN mamba create -y -n assembly \
    prokka flye busco

# RUN cat /etc/skel/.bashrc >> /etc/bash.bashrc

# intall required software with conda
# COPY environment.yml ./

# RUN conda env create \
#     -n ngs-tools \
#     -f environment.yml

RUN conda init 

COPY settings.json /config/data/User/
COPY keybindings.json /config/data/User/
RUN mkdir -p /config/workdir
