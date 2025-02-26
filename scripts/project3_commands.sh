#!/usr/bin/env bash

# use with conda environment:
# conda activate assembly

# for issues with prokka:
# cpanm Bio::SearchIO::hmmer --force

wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/project3.tar.gz
tar -xvf project3.tar.gz
rm project3.tar.gz

cd ~/project/project3

# bash function for performing the assembly and annotation
assembly_annotation () {
    SAMPLE=$1
    
    mkdir -p "$SAMPLE"_assembly

    # run flye with default options
    flye \
    --pacbio-hifi "$SAMPLE".fastq.gz \
    --threads 4 \
    --out-dir "$SAMPLE"_assembly

    # run BUSCO with auto lineage detection for prokaryotes (--auto-lineage-prok)
    # mode for genome and force overwrite (for re-running)
    busco \
    -i "$SAMPLE"_assembly/assembly.fasta \
    --auto-lineage-prok \
    --out_path "$SAMPLE"_assembly \
    --out BUSCO \
    --mode genome \
    --cpu 4 \
    --force

    # generate the BUSCO plot
    generate_plot.py \
    -wd "$SAMPLE"_assembly/BUSCO

    # run prokka for annotation
    prokka \
    --cpus 4 \
    --outdir "$SAMPLE"_assembly/prokka \
    "$SAMPLE"_assembly/assembly.fasta
}

# run the functions. Here an example for LWH7 on HiFi reads
assembly_annotation sample_1

