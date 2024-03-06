#!/usr/bin/env bash

# download data
cd ~/project
wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/project1.tar.gz
tar -xvf project1.tar.gz
rm project1.tar.gz

cd ~/project/project1

mkdir -p nanoplot

NanoPlot \
--fastq_rich reads/Cell_2.fastq.gz \
--outdir nanoplot/Cell_2

NanoPlot \
--fastq_rich reads/EV_2.fastq.gz \
--outdir nanoplot/EV_2

mkdir -p alignments

for sample in EV_2 Cell_2; do
    minimap2 \
    -a \
    -x splice \
    -t 4 \
    references/Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa \
    reads/"$sample".fastq.gz \
    | samtools sort \
    | samtools view -bh > alignments/"$sample".bam

    ## indexing for IGV
    samtools index alignments/"$sample".bam
done
