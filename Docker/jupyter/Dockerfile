FROM jupyter/minimal-notebook:e407f93c8dcc

USER root
RUN cat /etc/skel/.bashrc >> /etc/bash.bashrc
USER $NB_UID

RUN rm -r /home/jovyan/work

RUN /opt/conda/bin/conda config --add channels bioconda
RUN /opt/conda/bin/conda config --add channels conda-forge

RUN /opt/conda/bin/conda install mamba -y -n base -c conda-forge

RUN mamba install -y --quiet \
samtools>=1.10 \
minimap2 \
fastqc \
pbmm2 \
parallel

# install nanoplot with pip (with conda gives errors)
RUN pip install NanoPlot

RUN mamba create -y -n flair -c conda-forge -c bioconda flair

RUN mamba create -y -n pacbio \
    pbccs lima pbmm2 matplotlib numpy \
    pandas seaborn mappy pysam scikit-learn pbcore

RUN mamba create -y -n assembly \
    prokka flye busco
