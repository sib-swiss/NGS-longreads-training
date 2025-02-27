#!/usr/bin/env bash

cd ~/project

wget https://ngs-longreads-training.s3.eu-central-1.amazonaws.com/project2.tar.gz
tar -xvf project2.tar.gz
rm project2.tar.gz

cd ~/project/project2/

mkdir alignments

# loop over all samples and perform alignments
for sample in `seq 1015 1022`
do
  minimap2 \
  -a \
  -x asm20 \
  reference/Homo_sapiens.GRCh38.dna.primary_assembly.chrX.chr4.fa \
  reads/$sample.fastq.gz \
  | samtools sort \
  | samtools view -bh \
  > alignments/"$sample".bam

  samtools index alignments/"$sample".bam
done

mkdir -p trgt_output

for sample in `seq 1015 1022`
do
  trgt genotype --genome reference/Homo_sapiens.GRCh38.dna.primary_assembly.chrX.chr4.fa \
        --repeats targets/targets.bed \
        --reads alignments/"$sample".bam \
        --output-prefix trgt_output/"$sample" \
        --sample-name "$sample"

      samtools sort -o trgt_output/"$sample".spanning.sorted.bam trgt_output/"$sample".spanning.bam
      samtools index trgt_output/"$sample".spanning.sorted.bam
done

for sample in `seq 1015 1022`
do
  for gene in gene1 gene2
  do
    trgt plot --genome reference/Homo_sapiens.GRCh38.dna.primary_assembly.chrX.chr4.fa \
          --repeats targets/targets.bed \
          --vcf trgt_output/"$sample".vcf.gz \
          --spanning-reads trgt_output/"$sample".spanning.sorted.bam \
          --repeat-id "$gene" \
          --image trgt_output/"$gene".$sample.allele.svg
  done
done

# # clone the pacbio application scripts repo
# git clone https://github.com/PacificBiosciences/apps-scripts.git

# # for both genes generate the reports with makeReports.sh
# for gene in gene1 gene2
# do
#   apps-scripts/RepeatAnalysisTools/makeReports.sh \
#   targets/target_"$gene"_hg38.bed \
#   reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#   expansion_report_"$gene" \
#   alignments/????.bam
# done
