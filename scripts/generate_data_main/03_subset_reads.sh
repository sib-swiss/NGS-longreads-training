#!/usr/bin/env bash

#SBATCH --mail-user=geert.vangeest@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=subset_reads
#SBATCH --partition=pibu_el8
#SBATCH --array=1-8
#SBATCH --output=/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/slurm_logs/%x_%A_%a.out
#SBATCH --error=/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/slurm_logs/%x_%A_%a.err

# align the RNA-seq ONT reads with minimap2
module add minimap2/2.20-GCCcore-10.3.0
module add SAMtools/1.13-GCC-10.3.0

# define the input and output directories
PROJDIR="/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/"
bam_dir="$PROJDIR/alignments"

# define the input and output files
# take sample_info.csv and use the first column to get the file names. The csv file has headers, so we skip the first line
fastq_file=$(awk -F ',' 'NR>1 {print $1}' $PROJDIR/sample_info.csv | sed -n ${SLURM_ARRAY_TASK_ID}p)
input_file=$(basename $fastq_file .fastq.gz).bam
output_file=$(basename $input_file .bam).chr5_6_X.bam

# subset bam files on chromsome 5,6 and X
samtools view -b "$bam_dir"/"$input_file" 5:132802367-144318699 6:42465188-64136196 X:153161980-155636690 > "$bam_dir"/"$output_file"
samtools index "$bam_dir"/"$output_file"
