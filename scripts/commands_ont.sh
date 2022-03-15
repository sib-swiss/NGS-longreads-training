#!/usr/bin/env bash

cd ~/workdir/groupwork_ont/reads
# get a variable with sample names out of fastq file names
BASENAMES=`ls *.fastq.gz | cut -f 1 -d "."`

cd ..

# generate a minimap2 reference, so it is not rerun every time
minimap2 \
-x splice \
-d reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa.mmi \
reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa

# loop over the sample names
# and perform the alignment (set -x to splice for spliced alignment)
for name in $BASENAMES
do
  minimap2 \
  -a \
  -x splice \
  -G 500k \
  -t 4 \
  reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa.mmi \
  reads/$name.fastq.gz \
  | samtools sort \
  | samtools view -bh > alignments/$name.bam
done

# merge all alignments in order to annotate the splice variants together:
samtools merge alignments/merged.bam alignments/*.bam
samtools index alignments/merged.bam

# clone the FLAIR git repo
git clone https://github.com/BrooksLabUCSC/flair.git

# convert the bam file to a bed file
python flair/bin/bam2Bed12.py \
-i alignments/merged.bam \
> alignments/merged.bed

mkdir flair_output

# correct the bed files based on the gtf and genome
python flair/flair.py correct \
-q alignments/merged.bed \
-g reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa \
-f reference/Homo_sapiens.GRCh38.102.gtf \
-o flair_output/CACNA1C

# generate a comma-seperated list of fastq files
READS=`ls reads/*.fastq.gz | tr "\n" ","`
READS="${READS%?}" #remove last comma

# collapse the individual reads to corrected splice variants
python flair/flair.py collapse \
-g reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa \
-f reference/Homo_sapiens.GRCh38.102.gtf \
-r $READS \
-q flair_output/CACNA1C_all_corrected.bed \
-o flair_output/CACNA1C.collapse

# quantify the splice variants per sample
python flair/flair.py quantify \
-r reads_manifest.tsv \
-i flair_output/CACNA1C.collapse.isoforms.fa

# nice plot for general overview of isoform usage
python flair/bin/plot_isoform_usage.py \
flair_output/CACNA1C.collapse.isoforms.bed \
counts_matrix.tsv \
ENSG00000151067
