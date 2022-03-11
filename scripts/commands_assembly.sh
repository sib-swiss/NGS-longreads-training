#!/usr/bin/env bash

# cpanm Bio::SearchIO::hmmer --force
# for issues with prokka

assembly_annotation () {
    SAMPLE=$1
    READTYPE=$2

    cd /home/jovyan/workdir/groupwork_assembly/"$SAMPLE"
    
    if [ $READTYPE == "HiFi" ]
    then
        MODE="--pacbio-hifi"
    else
        MODE="--pacbio-raw"
    fi
    
    flye \
    "$MODE" "$READTYPE"/reads/"$SAMPLE"_"$READTYPE".fasta.gz \
    --threads 4 \
    --out-dir "$READTYPE"/assembly

    busco \
    -i "$READTYPE"/assembly/assembly.fasta \
    --auto-lineage-prok \
    --out_path "$READTYPE"/ \
    --out BUSCO \
    --mode genome \
    --cpu 4 \
    --force

    generate_plot.py \
    -wd "$READTYPE"/BUSCO

    prokka \
    --cpus 4 \
    --outdir "$READTYPE"/prokka \
    "$READTYPE"/assembly/assembly.fasta
}

# assembly_annotation LWH7 HiFi
assembly_annotation LWH7 CLR
