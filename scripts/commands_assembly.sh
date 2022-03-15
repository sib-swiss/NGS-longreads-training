#!/usr/bin/env bash

# use with conda environment:
# conda activate assembly

# for issues with prokka:
# cpanm Bio::SearchIO::hmmer --force

# bash function for performing the assembly and annotation
assembly_annotation () {
    SAMPLE=$1
    READTYPE=$2

    cd /home/jovyan/workdir/groupwork_assembly/"$SAMPLE"
    
    # if the readtype is HiFi take --pacbio-hifi as option to flye
    # if the readtype is something else (i.e. CLR) take --pacbio-raw as option to flye
    if [ $READTYPE == "HiFi" ]
    then
        MODE="--pacbio-hifi"
    else
        MODE="--pacbio-raw"
    fi
    
    # run flye with default options
    flye \
    "$MODE" "$READTYPE"/reads/"$SAMPLE"_"$READTYPE".fasta.gz \
    --threads 4 \
    --out-dir "$READTYPE"/assembly

    # run BUSCO with auto lineage detection for prokaryotes (--auto-lineage-prok)
    # mode for genome and force overwrite (for re-running)
    busco \
    -i "$READTYPE"/assembly/assembly.fasta \
    --auto-lineage-prok \
    --out_path "$READTYPE"/ \
    --out BUSCO \
    --mode genome \
    --cpu 4 \
    --force

    # generate the BUSCO plot
    generate_plot.py \
    -wd "$READTYPE"/BUSCO

    # run prokka for annotation
    prokka \
    --cpus 4 \
    --outdir "$READTYPE"/prokka \
    "$READTYPE"/assembly/assembly.fasta
}

# run the functions. Here an example for LWH7 on HiFi reads
assembly_annotation LWH7 HiFi
# and CLR reads:
assembly_annotation LWH7 CLR
