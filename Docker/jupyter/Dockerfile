FROM jupyter/minimal-notebook:9fe5186aba96

USER root
RUN cat /etc/skel/.bashrc >> /etc/bash.bashrc
USER $NB_UID

RUN rm -r /home/jovyan/work

RUN /opt/conda/bin/conda install -y --quiet -c bioconda \
samtools=1.10 \
minimap2=2.17 \
fastqc=0.11.9 \
pbmm2=1.4.0 \
parallel=20170422

RUN parallel --citation <<< 'will cite'

# install nanoplot with pip (with conda gives errors)
RUN pip install NanoPlot

# install flair environment from .yml file
COPY flair_env.yml ./
RUN /opt/conda/bin/conda env create -f flair_env.yml
RUN rm flair_env.yml

# install pacbio environment from .yml file
COPY pacbio_env.yml ./
RUN /opt/conda/bin/conda env create -f pacbio_env.yml
RUN rm pacbio_env.yml

# install assembly environment from .yml file
COPY assembly_env.yml ./
RUN /opt/conda/bin/conda env create -f assembly_env.yml
RUN rm assembly_env.yml
