FROM jupyter/base-notebook

RUN /opt/conda/bin/conda install -y --quiet -c bioconda \
samtools \
minimap2 \
fastqc \
pbmm2

RUN pip install NanoPlot
