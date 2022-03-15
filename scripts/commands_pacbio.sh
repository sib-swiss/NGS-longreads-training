#!/usr/bin/env bash

cd ~/workdir/groupwork_pacbio/

# generate reference for minimap2
minimap2 \
-x asm20 \
-d reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.mmi \
reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa

mkdir alignments

# loop over all samples and perform alignments
for sample in `seq 1015 1022`
do
  minimap2 \
  -a \
  -x asm20 \
  reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.mmi \
  reads/$sample.fastq.gz \
  | samtools sort \
  | samtools view -bh \
  > alignments/"$sample".bam

  samtools index alignments/"$sample".bam
done

# clone the pacbio application scripts repo
git clone https://github.com/PacificBiosciences/apps-scripts.git

# for both genes generate the reports with makeReports.sh
for gene in gene1 gene2
do
  apps-scripts/RepeatAnalysisTools/makeReports.sh \
  targets/target_"$gene"_hg38.bed \
  reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  expansion_report_"$gene" \
  alignments/????.bam
done
