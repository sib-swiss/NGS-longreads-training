#!/usr/bin/env bash

#SBATCH --mail-user=geert.vangeest@unibe.ch
#SBATCH --mail-type=fail
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=run_nanoplot
#SBATCH --partition=pibu_el8
#SBATCH --array=1-6
#SBATCH --output=/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/slurm_logs/%x_%A_%a.out
#SBATCH --error=/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/slurm_logs/%x_%A_%a.err

# define the input and output directories
PROJDIR="/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/"
fastq_dir="$PROJDIR/subset_fastq"

cd $fastq_dir
fastq_file=$(ls *.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

out_dir="$PROJDIR"/nanoplot/$(basename $fastq_file .fastq.gz)

mkdir -p $out_dir

NanoPlot --fastq_rich $fastq_dir/$fastq_file \
--outdir $out_dir
