#!/usr/bin/env bash

#SBATCH --mail-user=geert.vangeest@unibe.ch
#SBATCH --mail-type=fail
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=prepare_ref_files
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/slurm_logs/%x_%A_%a.out
#SBATCH --error=/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/slurm_logs/%x_%A_%a.err

PROJDIR="/data/users/gvangeest/repositories/NGS-longreads-training/generate_data_main/"
cd $PROJDIR

GTF=/data/references/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/build111/Homo_sapiens.GRCh38.111.gtf
FASTA=/data/references/Homo_sapiens/Ensembl/GRCh38/Sequence/Homo_sapiens.GRCh38.dna.primary_assembly.fa

grep ^# $GTF > references/Homo_sapiens.GRCh38.111.chr5.chr6.chrX.gtf

bedtools intersect -wa -a \
$GTF \
-b references/regions.bed \
>> references/Homo_sapiens.GRCh38.111.chr5.chr6.chrX.gtf

echo "5
6
X" > references/chrom_of_interest.txt

apptainer exec \
/data/users/gvangeest/containers/seqtk:1.4--he4a0461_1 \
seqtk subseq $FASTA references/chrom_of_interest.txt \
> references/Homo_sapiens.GRCh38.dna.primary_assembly.chr5.chr6.chrX.fa
