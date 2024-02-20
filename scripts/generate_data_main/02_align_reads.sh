#!/usr/bin/env bash

#SBATCH --mail-user=geert.vangeest@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=12
#SBATCH --job-name=align_reads
#SBATCH --partition=pibu_el8
#SBATCH --array=1-8
#SBATCH --output=/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/slurm_logs/%x_%A_%a.out
#SBATCH --error=/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/slurm_logs/%x_%A_%a.err

# align the RNA-seq ONT reads with minimap2
module add minimap2/2.20-GCCcore-10.3.0
module add SAMtools/1.13-GCC-10.3.0

# define the input and output directories
PROJDIR="/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/"
input_dir="$PROJDIR/raw_data"
output_dir="$PROJDIR/alignments"

mkdir -p $output_dir

# define the input and output files
# take sample_info.csv and use the first column to get the file names. The csv file has headers, so we skip the first line
input_file=$(awk -F ',' 'NR>1 {print $1}' $PROJDIR/sample_info.csv | sed -n ${SLURM_ARRAY_TASK_ID}p)
output_file=$(basename $input_file .fastq.gz).bam

# align the reads
# minimap2 -ax splice -t ${SLURM_CPUS_PER_TASK} \
# /data/references/Homo_sapiens/Ensembl/GRCh38/Sequence/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
# "$input_dir"/"$input_file" \
# | samtools view -bS - \
# | samtools sort -o "$output_dir"/"$output_file"

samtools index "$output_dir"/"$output_file"
