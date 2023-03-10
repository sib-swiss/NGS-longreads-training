#!/usr/bin/env bash

mkdir ~/workdir

cd ~/workdir

NanoPlot \
--fastq data/reads/cerebellum-5238-batch2.fastq.gz \
--outdir nanoplot_output

mkdir ~/workdir/alignments

cd ~/workdir

minimap2 \
-a \
-x splice \
-G 500k \
-t 4 \
data/reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa \
data/reads/cerebellum-5238-batch2.fastq.gz \
| samtools sort \
| samtools view -bh > alignments/cerebellum-5238-batch2.bam

## indexing for IGV
samtools index alignments/cerebellum-5238-batch2.bam
