FROM jupyter/minimal-notebook

USER root
RUN cat /etc/skel/.bashrc >> /etc/bash.bashrc
USER $NB_UID

RUN rm -r /home/jovyan/work

RUN /opt/conda/bin/conda install -y --quiet -c bioconda \
samtools \
minimap2 \
fastqc \
pbmm2 \
parallel

RUN parallel --citation <<< 'will cite'

RUN pip install NanoPlot

COPY flair_env.yml ./

RUN /opt/conda/bin/conda env create -f flair_env.yml

RUN rm flair_env.yml

COPY pacbio_env.yml ./

RUN /opt/conda/bin/conda env create -f pacbio_env.yml

RUN rm pacbio_env.yml
